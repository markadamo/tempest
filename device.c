/*
 *  device.c
 *
 *  Created by Brendan Faherty on 9/24/10.
 *  Copyright 2010 Dartmouth College. All rights reserved.
 *
 */
#include "tempest.h"

cl_mem cl_iPeakCounts;
cl_mem cl_lPeakIndices;
cl_mem cl_iPeakBins;
cl_mem cl_fPeakInts;
cl_mem cl_fSpectra;
cl_mem cl_fScratch;
cl_mem cl_init_bTheoData;
cl_mem cl_init_fSpectra;
cl_mem cl_init_fScratch;
cl_mem cl_bTheoData;
cl_mem cl_cCandidates;
cl_mem cl_fScores;
cl_mem cl_mPSMs;
cl_mem cl_fNextScores;
cl_mem cl_nlValuesNterm;
cl_mem cl_nlValuesCterm;
cl_mem cl_nlValuesAA;
std::stack<cl_mem> unusedBuffers;

cl_context clContext;
cl_platform_id clPlatformID;
cl_device_id clDeviceID;
cl_program clProgram;
cl_kernel clKernel;
cl_command_queue clCommandQueue;

/*
 * Set the GPU device to use, and check for compatibility.
 */

extern void initialize_device() {
    cl_platform_id platforms[100];
    cl_uint platforms_n = 0;
    int err;

    err = clGetPlatformIDs(100, platforms, &platforms_n);
    check_cl_error(__FILE__, __LINE__, err, "Unable to enumerate OpenCL platforms");

    if (platforms_n == 0) {
        printf("  No OpenCL platforms found!\n");
        printf("  Make sure OpenCL drivers are installed for your hardware.\n");
        tempest_exit(EXIT_FAILURE);
    }
    
    if (config.iPlatform >= platforms_n) {
        printf("  Selected platform ID (%d) is out of range!\n", config.iPlatform);
        printf("  Available OpenCL platform IDs range from 0 to %d.\n", platforms_n-1);
        tempest_exit(EXIT_FAILURE);
    }

    clPlatformID = platforms[config.iPlatform];
    
    cl_device_id devices[100];
    cl_uint devices_n = 0;
    err = clGetDeviceIDs(clPlatformID, CL_DEVICE_TYPE_ALL, 100, devices, &devices_n);
    check_cl_error(__FILE__, __LINE__, err, "Unable to enumerate OpenCL devices");

    if (config.iDevice >= devices_n) {
        printf("  Selected device ID (%d) is out of range!\n", config.iDevice);
        printf("  Available OpenCL device IDs on platform %d range from 0 to %d.\n", config.iPlatform, devices_n-1);
        tempest_exit(EXIT_FAILURE);
    }

    clDeviceID = devices[config.iDevice];
    
    char buffer[10240];
    clGetDeviceInfo(clDeviceID, CL_DEVICE_NAME, sizeof(buffer), buffer, NULL);
    printf(" » Using %s\n", buffer);

    //create context
    clContext = clCreateContext(NULL, 1, devices+config.iDevice, NULL, NULL, &err);
    check_cl_error(__FILE__, __LINE__, err, "Error creating context");

    /* Create a command queue */
    clCommandQueue = clCreateCommandQueue(clContext, clDeviceID, CL_QUEUE_PROFILING_ENABLE, &err);
    if(err < 0) {
        perror("Couldn't create a command queue");
        exit(1);   
    }
    
}


/*
 *  Prepare the device for scoring.
 *
 *   - create events
 *   - copy constant memory
 *   - determine configuration
 *   - setup global device memory (allocate and copy)
 *   - prebuild observed MS/MS spectra
 *   - free host memory
 */

extern void setup_device() {
    size_t zMemFree, zMemTotal, zMemRequired;
    size_t zMemSpectrum, zMemPeaksInfo, zMemPeaks, zMemSpectrumb, zMemCandidates, zMemResults;
    size_t zMemSetup_noprebuild, zMemSetup_parallel, zMemSetup_sequential;
    size_t zMemScore_prebuilt, zMemScore_peaks, zMemScore_peaks_shared;
    int err;
    
    err  = clGetDeviceInfo(clDeviceID, CL_DEVICE_MAX_MEM_ALLOC_SIZE,  sizeof(cl_ulong), &(cl_info.lMaxMemAllocSize),  NULL);
    err |= clGetDeviceInfo(clDeviceID, CL_DEVICE_ADDRESS_BITS,        sizeof(cl_uint),  &(cl_info.iAddressBits),      NULL);
    err |= clGetDeviceInfo(clDeviceID, CL_DEVICE_HOST_UNIFIED_MEMORY, sizeof(cl_bool),  &(cl_info.bUnifiedMemory),      NULL); 
    err |= clGetDeviceInfo(clDeviceID, CL_DEVICE_GLOBAL_MEM_SIZE,     sizeof(cl_ulong), &(cl_info.lGlobalMemSize),    NULL);  
    err |= clGetDeviceInfo(clDeviceID, CL_DEVICE_MAX_COMPUTE_UNITS,   sizeof(cl_uint),  &(cl_info.iMaxComputeUnits),  NULL);
    err |= clGetDeviceInfo(clDeviceID, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(cl_ulong), &(cl_info.lMaxWorkGroupSize), NULL);
    check_cl_error(__FILE__, __LINE__, err, "Unable to get OpenCL device attributes");

    //limit global memory size to addressable memory (e.g. NVidia may list more than is addressable)
    if (cl_info.iAddressBits == 32)
        if (cl_info.lGlobalMemSize > pow(2, 32))
            cl_info.lGlobalMemSize = pow(2, 32);
    
    MAX_BLOCKDIM = cl_info.lMaxWorkGroupSize;
    //DEFAULT_BLOCKDIM_SCORE = cl_info.lMaxWorkGroupSize > DEFAULT_CANDIDATE_BUFFER_SIZE ? DEFAULT_CANDIDATE_BUFFER_SIZE : cl_info.lMaxWorkGroupSize;
    //BLOCKDIM_BUILD = cl_info.lMaxWorkGroupSize;
    BLOCKDIM_TRANSFORM = cl_info.lMaxWorkGroupSize;
    BLOCKDIM_REDUCE = cl_info.lMaxWorkGroupSize > DEFAULT_CANDIDATE_BUFFER_SIZE/2 ? DEFAULT_CANDIDATE_BUFFER_SIZE/2 : cl_info.lMaxWorkGroupSize;
    zMemFree = cl_info.lGlobalMemSize;

    /* Build program */
    clProgram = build_program(clContext, clDeviceID, KERNEL_FILE);

    // Create kernels
    create_kernels();
    
    // Determine Configuration
    
    // peaks
    size_t size_iPeakBins = tempest.lNumMS2Peaks * sizeof(cl_int);
    size_t size_fPeakInts = tempest.lNumMS2Peaks * sizeof(cl_float);
    cl_iPeakBins = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, size_iPeakBins, &(host_iPeakBins[0]), &err);
    check_cl_error(__FILE__, __LINE__, err, "Unable to allocate device memory for peak bins.");
    cl_fPeakInts = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, size_fPeakInts, &(host_fPeakInts[0]), &err);
    check_cl_error(__FILE__, __LINE__, err, "Unable to allocate device memory for peak intensities.");
    
    // cleanup host
    std::vector<int>().swap(host_iPeakBins);
    std::vector<float>().swap(host_fPeakInts);

    //cudaMalloc((void**) &gpu_fSpectra, tempest.iNumMS2Bins * sizeof(float));
    //cl_fSpectra = clCreateBuffer(clContext, CL_MEM_READ_WRITE, tempest.iNumMS2Bins * sizeof(float), NULL, &err);
    float * init_fSpectra = (float *) calloc(tempest.iNumMS2Bins, sizeof(float));
    size_t size_init_fSpectra = tempest.iNumMS2Bins * sizeof(cl_float);
    cl_init_fSpectra = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, size_init_fSpectra, init_fSpectra, &err);
    free(init_fSpectra);
    
    // candidate and results
    mObj * init_mPSMs = (mObj *) calloc(tempest.iNumSpectra * params.iNumOutputPSMs, sizeof(mObj)); 
    float * init_fNextScores = (float *) calloc(tempest.iNumSpectra, sizeof(float));
    size_t size_cCandidates = sizeof(cObj) * config.iCandidateBufferSize;
    size_t size_fScores = sizeof(cl_float)  * config.iCandidateBufferSize;
    size_t size_mPSMs = sizeof(mObj)  * tempest.iNumSpectra * params.iNumOutputPSMs;
    size_t size_fNextScores = sizeof(float) * tempest.iNumSpectra;
    cl_cCandidates = clCreateBuffer(clContext, CL_MEM_READ_ONLY, size_cCandidates, NULL, &err);
    cl_fScores = clCreateBuffer(clContext, CL_MEM_READ_WRITE, size_fScores, NULL, &err);  
    cl_mPSMs = clCreateBuffer(clContext, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, size_mPSMs , init_mPSMs, &err);
    cl_fNextScores = clCreateBuffer(clContext, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, size_fNextScores, init_fNextScores, &err);
    //MEA: need to block free until previous clCreateBuffer commands complete?
    free(init_mPSMs);
    free(init_fNextScores);
    check_cl_error(__FILE__, __LINE__, err, "Unable to allocate device memory for candidates and results.");

    //determine how many spectra can be kept in device memory at a time
    size_t availMemSpectra = cl_info.lGlobalMemSize
        - size_iPeakBins
        - size_fPeakInts
        - size_init_fSpectra
        - size_cCandidates
        - size_fScores
        - size_mPSMs
        - size_fNextScores;
    long maxCachedSpectra = availMemSpectra / (tempest.iNumMS2Bins*sizeof(cl_float));
    //maxCachedSpectra = 1;
    if (maxCachedSpectra > tempest.iNumSpectra)
        maxCachedSpectra = tempest.iNumSpectra;
    
    printf(" » Allocating %ld bytes of device memory for %ld cached spectra.\n", maxCachedSpectra*tempest.iNumMS2Bins*sizeof(cl_float), maxCachedSpectra);
    for (int i=0; i<maxCachedSpectra; i++) {
        cl_mem newBuffer = clCreateBuffer(clContext, CL_MEM_READ_WRITE, tempest.iNumMS2Bins*sizeof(cl_float), NULL, &err);
        check_cl_error(__FILE__, __LINE__, err, "Unable to allocate spectrum memory on device.");
        unusedBuffers.push(newBuffer);
    }
}

/*
 * Cleanup device memory and events. Prints any timing requests.
 */

extern void cleanup_device()
{
    int err;
    err  = clReleaseMemObject(cl_lPeakIndices);
    err |= clReleaseMemObject(cl_iPeakCounts);
    err |= clReleaseMemObject(cl_iPeakBins);
    err |= clReleaseMemObject(cl_fPeakInts);
    err |= clReleaseMemObject(cl_fSpectra);
    err |= clReleaseMemObject(cl_fScratch);
    err |= clReleaseMemObject(cl_cCandidates);
    err |= clReleaseMemObject(cl_fScores);
    err |= clReleaseMemObject(cl_mPSMs);
    err |= clReleaseMemObject(cl_fNextScores);
    err |= clReleaseMemObject(cl_bTheoData);
    check_cl_error(__FILE__, __LINE__, err, "Unable to free memory object(s)");
}

/* Create program from a file and compile it */
cl_program build_program(cl_context ctx, cl_device_id dev, const char* filename) {

    cl_program program;
    FILE *program_handle;
    char *program_buffer, *program_log;
    size_t program_size, log_size;
    char options[1024];
    int err;
   

    /* Read program file and place content into buffer */
    program_handle = fopen(filename, "r");
    if(program_handle == NULL) {
        perror("Couldn't find the program file");
        exit(1);
    }
    fseek(program_handle, 0, SEEK_END);
    program_size = ftell(program_handle);
    rewind(program_handle);
    program_buffer = (char*)malloc(program_size + 1);
    program_buffer[program_size] = '\0';
    fread(program_buffer, sizeof(char), program_size, program_handle);
    fclose(program_handle);

    /* Create program from file */
    program = clCreateProgramWithSource(ctx, 1, 
                                        (const char**)&program_buffer, &program_size, &err);
    if(err < 0) {
        perror("Couldn't create the program");
        exit(1);
    }
    free(program_buffer);

    /* Define known constants as literals to build into the program */
    float fMassProteinNterm = (float) params.dProteinNtermMass;
    float fMassPeptideNterm = (float) params.dPeptideNtermMass;
    float fMassNtermMod = (float) params.dVariableNtermMassDiff;
    sprintf(options,                   "-DMAX_PEPTIDE_LENGTH=%d ",     MAX_PEPTIDE_LENGTH);
    sprintf(options + strlen(options), "-DBLOCKDIM_REDUCE=%d ",        BLOCKDIM_REDUCE);
    sprintf(options + strlen(options), "-DBLOCKDIM_SCORE=%d " ,        DEFAULT_BLOCKDIM_SCORE);
    sprintf(options + strlen(options), "-DGPU_MASS_PROTEIN_NTERM=%f ", fMassProteinNterm);
    sprintf(options + strlen(options), "-DGPU_MASS_PEPTIDE_NTERM=%f ", fMassPeptideNterm);
    sprintf(options + strlen(options), "-DGPU_MASS_NTERM_MOD=%f ",     fMassNtermMod);
    sprintf(options + strlen(options), "-DGPU_TOLERANCE=%f ",          params.fFragmentTolerance);
    sprintf(options + strlen(options), "-DGPU_NUM_OUTPUT_PSMS=%d ",    params.iNumOutputPSMs);
    sprintf(options + strlen(options), "-DGPU_TRACK_DUPLICATES=%d ",   params.bTrackDuplicates);
    sprintf(options + strlen(options), "-DGPU_FLANKING=%d ",           params.bTheoreticalFlanking);
    sprintf(options + strlen(options), "-DNUM_AA_NL=%d ",              params.numAANL);
    sprintf(options + strlen(options), "-DNUM_NTERM_NL=%d ",           params.numNtermNL);
    sprintf(options + strlen(options), "-DNUM_CTERM_NL=%d ",           params.numCtermNL);
    sprintf(options + strlen(options), "-DGPU_FIX_DELTA_SCORE=%d ",    params.bFixDeltaScore);
    sprintf(options + strlen(options), "-DGPU_XWIDTH=%d ",             tempest.iCrossCorrelationWidth);
    sprintf(options + strlen(options), "-DGPU_NUM_BINS=%d ",           tempest.iNumMS2Bins);
    //sprintf(options + strlen(options), "-cl-nv-maxrregcount=16");
    sprintf(options + strlen(options), "-cl-single-precision-constant");
   
    /* Build program */
    err = clBuildProgram(program, 0, NULL, options, NULL, NULL);
    if(PROFILE || err < 0) {

        /* Find size of log and print to std output */
        clGetProgramBuildInfo(program, dev, CL_PROGRAM_BUILD_LOG, 
                              0, NULL, &log_size);
        program_log = (char*) malloc(log_size + 1);
        program_log[log_size] = '\0';
        clGetProgramBuildInfo(program, dev, CL_PROGRAM_BUILD_LOG, 
                              log_size + 1, program_log, NULL);
        printf("%s\n", program_log);
        free(program_log);
        if (err < 0)
            exit(1);
    }

    return program;
}

