#include "tempest.hpp"

/*
 * Launches a kernel to perform parallel memset on a block of device memory
 */
int Device::cl_memset(cl_mem buffer, int c, unsigned long n, cl_event* evt) {
    int err;
    size_t global_dim = Tempest::mround(n, memset_size);
    err  = clSetKernelArg(__cl_memset, 0, sizeof(cl_mem), &buffer);
    err |= clSetKernelArg(__cl_memset, 1, sizeof(int), &c);
    err |= clSetKernelArg(__cl_memset, 2, sizeof(long), &n);
    err |= clEnqueueNDRangeKernel(clCommandQueue, __cl_memset, 1, NULL, &global_dim, &memset_size, 0, NULL, evt);
    return err;
}

void Device::setup_constant_memory() {
    float fMassAA[256];
    int err;

    // convert the double masses to floats
    for (int iResidue=0; iResidue < 256; iResidue++) {
        fMassAA[iResidue] = (float) Tempest::dMassAA[iResidue];
    }

    cl_mem GPU_MASS_AA      = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, 256*sizeof(float), &fMassAA, &err);
    cl_mem UN_MOD_AA        = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, 256*sizeof(char),  &Tempest::unModAA, &err);
    cl_mem NTERM_MOD_MASSES = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, 5*sizeof(float),  &Tempest::ntermModMasses, &err);
    cl_mem CTERM_MOD_MASSES = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, 5*sizeof(float),  &Tempest::ctermModMasses, &err);
    
    //neutral losses
    //transfer to device
    if (Tempest::params.numNtermNL)
        cl_nlValuesNterm = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(nlValue)*Tempest::nlValuesNterm.size(), &Tempest::nlValuesNterm[0], &err);
    if (Tempest::params.numCtermNL)
        cl_nlValuesCterm = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(nlValue)*Tempest::nlValuesCterm.size(), &Tempest::nlValuesCterm[0], &err);
    if (Tempest::params.numAANL)
        cl_nlValuesAA = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(nlValue)*Tempest::nlValuesAA.size(), &Tempest::nlValuesAA[0], &err);
    Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to copy device constant memory");
    
    //pre-set these constant arrays as arguments to scoring kernels (re-used in every launch)
    err  = clSetKernelArg(__gpu_build, 1, sizeof(int),    &Tempest::tempest.iNumMS2Bins);
    err |= clSetKernelArg(__gpu_build, 3, sizeof(cl_mem), &cl_iPeakBins);
    err |= clSetKernelArg(__gpu_build, 4, sizeof(cl_mem), &cl_fPeakInts);
    Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to set constant args (__gpu_build)");

    err = clSetKernelArg(__gpu_transform, 1, sizeof(int),  &Tempest::tempest.iNumMS2Bins);
    Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to set constant args (__gpu_transform)");

    err  = clSetKernelArg(__gpu_score, 2, sizeof(cl_mem), &cl_cCandidates);
    err |= clSetKernelArg(__gpu_score, 3, sizeof(cl_mem), &cl_fScores);
    err |= clSetKernelArg(__gpu_score, 6, sizeof(int),    &Tempest::tempest.iNumMS2Bins);
    err |= clSetKernelArg(__gpu_score, 7, sizeof(cl_mem), &GPU_MASS_AA);
    err |= clSetKernelArg(__gpu_score, 8, sizeof(cl_mem), &NTERM_MOD_MASSES);
    err |= clSetKernelArg(__gpu_score, 9, sizeof(cl_mem), &CTERM_MOD_MASSES);
    err |= clSetKernelArg(__gpu_score, 10, sizeof(cl_mem), Tempest::params.numNtermNL ? &cl_nlValuesNterm : NULL);
    err |= clSetKernelArg(__gpu_score, 11, sizeof(cl_mem), Tempest::params.numCtermNL ? &cl_nlValuesCterm : NULL);
    err |= clSetKernelArg(__gpu_score, 12, sizeof(cl_mem), Tempest::params.numAANL ? &cl_nlValuesAA : NULL);
    Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to set constant args (__gpu_score)");

    size_t zMemShared = this->score_reduction_size * (sizeof(int) + sizeof(float));
    err  = clSetKernelArg(__gpu_score_reduction, 0, sizeof(int), &(this->candidateBufferSize));
    err |= clSetKernelArg(__gpu_score_reduction, 1, sizeof(cl_mem), &cl_cCandidates);
    err |= clSetKernelArg(__gpu_score_reduction, 2, sizeof(cl_mem), &cl_fScores);
    err |= clSetKernelArg(__gpu_score_reduction, 3, sizeof(cl_mem), &cl_mPSMs);
    err |= clSetKernelArg(__gpu_score_reduction, 6, zMemShared, NULL);
    Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to set constant args (__gpu_score_reduction)");
    
}

void Device::create_kernels() {
    int err;
    size_t paramBuffer;
    
    this->__gpu_build = clCreateKernel(clProgram, "gpu_build", &err);
    Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to create kernel");
    //clGetKernelWorkGroupInfo(__gpu_build, this->clDeviceID, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &paramBuffer, NULL);
    //printf("gpu_build: WORK_GROUP_SIZE=%ld\n", paramBuffer);
    //clGetKernelWorkGroupInfo(__gpu_build, this->clDeviceID, CL_KERNEL_LOCAL_MEM_SIZE, sizeof(size_t), &paramBuffer, NULL);
    //printf("gpu_build: LOCAL_MEM_SIZE=%ld\n", paramBuffer);
    clGetKernelWorkGroupInfo(__gpu_build, this->clDeviceID, CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, sizeof(size_t), &(this->build_size), NULL);
    //printf("gpu_build: PREFERRED_WORK_GROUP_SIZE_MULTIPLE=%ld\n", this->build_size);
    //clGetKernelWorkGroupInfo(__gpu_build, this->clDeviceID, CL_KERNEL_PRIVATE_MEM_SIZE, sizeof(size_t), &paramBuffer, NULL);
    //printf("gpu_build: PRIVATE_MEM_SIZE=%ld\n", paramBuffer);
    
    this->__gpu_transform = clCreateKernel(clProgram, "gpu_transform", &err);
    Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to create kernel");
    // clGetKernelWorkGroupInfo(__gpu_transform, this->clDeviceID, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &paramBuffer, NULL);
    // printf("gpu_transform: WORK_GROUP_SIZE=%ld\n", paramBuffer);
    // clGetKernelWorkGroupInfo(__gpu_transform, this->clDeviceID, CL_KERNEL_LOCAL_MEM_SIZE, sizeof(size_t), &paramBuffer, NULL);
    // printf("gpu_transform: LOCAL_MEM_SIZE=%ld\n", paramBuffer);
    clGetKernelWorkGroupInfo(__gpu_transform, this->clDeviceID, CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, sizeof(size_t), &(this->transform_size), NULL);
    //printf("gpu_transform: PREFERRED_WORK_GROUP_SIZE_MULTIPLE=%ld\n", paramBuffer);
    //clGetKernelWorkGroupInfo(__gpu_transform, this->clDeviceID, CL_KERNEL_PRIVATE_MEM_SIZE, sizeof(size_t), &paramBuffer, NULL);
    //printf("gpu_transform: PRIVATE_MEM_SIZE=%ld\n", paramBuffer);
    
    this->__gpu_score = clCreateKernel(clProgram, "gpu_score", &err);
    Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to create kernel");
    // clGetKernelWorkGroupInfo(__gpu_score, this->clDeviceID, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &paramBuffer, NULL);
    // printf("gpu_score: WORK_GROUP_SIZE=%ld\n", paramBuffer);
    // clGetKernelWorkGroupInfo(__gpu_score, this->clDeviceID, CL_KERNEL_LOCAL_MEM_SIZE, sizeof(size_t), &paramBuffer, NULL);
    // printf("gpu_score: LOCAL_MEM_SIZE=%ld\n", paramBuffer);
    clGetKernelWorkGroupInfo(__gpu_score, this->clDeviceID, CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, sizeof(size_t), &(this->score_size), NULL);
    //this->score_size = 1;
    // printf("gpu_score: PREFERRED_WORK_GROUP_SIZE_MULTIPLE=%ld\n", paramBuffer);
    // clGetKernelWorkGroupInfo(__gpu_score, this->clDeviceID, CL_KERNEL_PRIVATE_MEM_SIZE, sizeof(size_t), &paramBuffer, NULL);
    // printf("gpu_score: PRIVATE_MEM_SIZE=%ld\n", paramBuffer);
    
    this->__gpu_score_reduction = clCreateKernel(clProgram, "gpu_score_reduction", &err);
    Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to create kernel");
    clGetKernelWorkGroupInfo(__gpu_score_reduction, this->clDeviceID, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &(this->score_reduction_size_max), NULL);
    //printf("gpu_score_reduction: WORK_GROUP_SIZE=%ld\n", this->score_reduction_size_max);
    clGetKernelWorkGroupInfo(__gpu_score_reduction, this->clDeviceID, CL_KERNEL_LOCAL_MEM_SIZE, sizeof(size_t), &(this->score_reduction_size_local), NULL);
    //printf("gpu_score_reduction: LOCAL_MEM_SIZE=%ld\n", this->score_reduction_size_local);
    clGetKernelWorkGroupInfo(__gpu_score_reduction, this->clDeviceID, CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, sizeof(size_t), &this->score_reduction_multiple, NULL);
    //printf("gpu_score_reduction: PREFERRED_WORK_GROUP_SIZE_MULTIPLE=%ld\n", this->score_reduction_multiple);
    clGetKernelWorkGroupInfo(__gpu_score_reduction, this->clDeviceID, CL_KERNEL_PRIVATE_MEM_SIZE, sizeof(size_t), &paramBuffer, NULL);
    //printf("gpu_score_reduction: PRIVATE_MEM_SIZE=%ld\n", paramBuffer);
     
    this->__cl_memset = clCreateKernel(clProgram, "cl_memset", &err);
    Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to create kernel");
    // clGetKernelWorkGroupInfo(__cl_memset, this->clDeviceID, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &paramBuffer, NULL);
    // printf("cl_memset: WORK_GROUP_SIZE=%ld\n", paramBuffer);
    // clGetKernelWorkGroupInfo(__cl_memset, this->clDeviceID, CL_KERNEL_LOCAL_MEM_SIZE, sizeof(size_t), &paramBuffer, NULL);
    // printf("cl_memset: LOCAL_MEM_SIZE=%ld\n", paramBuffer);
    clGetKernelWorkGroupInfo(__cl_memset, this->clDeviceID, CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, sizeof(size_t), &(this->memset_size), NULL);
    // printf("cl_memset: PREFERRED_WORK_GROUP_SIZE_MULTIPLE=%ld\n", paramBuffer);
    // clGetKernelWorkGroupInfo(__cl_memset, this->clDeviceID, CL_KERNEL_PRIVATE_MEM_SIZE, sizeof(size_t), &paramBuffer, NULL);
    // printf("cl_memset: PRIVATE_MEM_SIZE=%ld\n", paramBuffer);
    
}

Device::Device(int platformID, int deviceInd) {
    cl_platform_id platforms[100];
    cl_uint platforms_n = 0;
    int err;

    this->deviceID = Tempest::config.iDevices[deviceInd];
    this->platformID = platformID;
    this->deviceInd = deviceInd;

    err = clGetPlatformIDs(100, platforms, &platforms_n);
    Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to enumerate OpenCL platforms");

    if (platforms_n == 0) {
        printf("  No OpenCL platforms found!\n");
        printf("  Make sure OpenCL drivers are installed for your hardware.\n");
        Tempest::tempest_exit(EXIT_FAILURE);
    }
    
    if (platformID >= platforms_n) {
        printf("  Selected platform ID (%d) is out of range!\n", platformID);
        printf("  Available OpenCL platform IDs range from 0 to %d.\n", platforms_n-1);
        Tempest::tempest_exit(EXIT_FAILURE);
    }

    clPlatformID = platforms[platformID];
    
    cl_device_id devices[100];
    cl_uint devices_n = 0;
    err = clGetDeviceIDs(clPlatformID, CL_DEVICE_TYPE_ALL, 100, devices, &devices_n);
    Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to enumerate OpenCL devices");

    if (deviceID >= devices_n) {
        printf("  Selected device ID (%d) is out of range!\n", deviceID);
        printf("  Available OpenCL device IDs on platform %d range from 0 to %d.\n", platformID, devices_n-1);
        Tempest::tempest_exit(EXIT_FAILURE);
    }

    clDeviceID = devices[deviceID];
    
    char buffer[10240];
    clGetDeviceInfo(clDeviceID, CL_DEVICE_NAME, sizeof(buffer), buffer, NULL);
    printf(" » (%d:%d) Using %s\n", platformID, deviceID, buffer);

    //create context
    clContext = clCreateContext(NULL, 1, devices+deviceID, NULL, NULL, &err);
    Tempest::check_cl_error(__FILE__, __LINE__, err, "Error creating context");

    /* Create a command queue */
    clCommandQueue = clCreateCommandQueue(clContext, clDeviceID, PROFILE ? CL_QUEUE_PROFILING_ENABLE : 0, &err);
    if(err < 0) {
        perror("Couldn't create a command queue");
        exit(1);   
    }
    
    err  = clGetDeviceInfo(clDeviceID, CL_DEVICE_MAX_MEM_ALLOC_SIZE,  sizeof(cl_ulong), &(lMaxMemAllocSize),  NULL);
    err |= clGetDeviceInfo(clDeviceID, CL_DEVICE_ADDRESS_BITS,        sizeof(cl_uint),  &(iAddressBits),      NULL);
    err |= clGetDeviceInfo(clDeviceID, CL_DEVICE_HOST_UNIFIED_MEMORY, sizeof(cl_bool),  &(bUnifiedMemory),    NULL); 
    err |= clGetDeviceInfo(clDeviceID, CL_DEVICE_GLOBAL_MEM_SIZE,     sizeof(cl_ulong), &(lGlobalMemSize),    NULL);
    err |= clGetDeviceInfo(clDeviceID, CL_DEVICE_LOCAL_MEM_SIZE,      sizeof(cl_ulong), &(lLocalMemSize),    NULL);
    err |= clGetDeviceInfo(clDeviceID, CL_DEVICE_MAX_COMPUTE_UNITS,   sizeof(cl_uint),  &(iMaxComputeUnits),  NULL);
    err |= clGetDeviceInfo(clDeviceID, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(cl_ulong), &(lMaxWorkGroupSize), NULL);
    Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to get OpenCL device attributes");

    //limit global memory size to addressable memory (e.g. NVidia may list more than is addressable)
    if (iAddressBits == 32)
        if (lGlobalMemSize > pow(2, 32))
            lGlobalMemSize = pow(2, 32);

    FILE *program_handle;
    char *program_buffer, *program_log;
    size_t program_size, log_size;
    char options[1024];

    /* Read program file and place content into buffer */
    program_handle = fopen(KERNEL_FILE, "r");
    if(program_handle == NULL) {
        perror("Couldn't find the .cl program file");
        exit(1);
    }
    fseek(program_handle, 0, SEEK_END);
    program_size = ftell(program_handle);
    rewind(program_handle);
    program_buffer = (char*)malloc(program_size + 1);
    program_buffer[program_size] = '\0';
    int bytesRead = fread(program_buffer, sizeof(char), program_size, program_handle);
    fclose(program_handle);

    /* Create program from file */
    clProgram = clCreateProgramWithSource(clContext, 1, (const char**)&program_buffer, &program_size, &err);
    
    //clProgram = clCreateProgramWithSource(clContext, 1, &kernels_cl_const, &kernels_cl_len, &err);
    if(err < 0) {
        perror("Couldn't create the .cl program");
        exit(1);
    }
    free(program_buffer);

    /* Define known constants as literals to build into the program */
    float fMassProteinNterm = (float) Tempest::params.dProteinNtermMass;
    float fMassPeptideNterm = (float) Tempest::params.dPeptideNtermMass;
    sprintf(options,                   "-DPRIMARY_INTENSITY=%d ",      PRIMARY_INTENSITY);
    sprintf(options + strlen(options), "-DPROTON_MASS=%f ",            PROTON_MASS);
    sprintf(options + strlen(options), "-DH_MASS=%f ",                 H_MASS);
    sprintf(options + strlen(options), "-DC_MASS=%f ",                 C_MASS);
    sprintf(options + strlen(options), "-DO_MASS=%f ",                 O_MASS);
    sprintf(options + strlen(options), "-DN_MASS=%f ",                 N_MASS);
    sprintf(options + strlen(options), "-DOH_MASS=%f ",                OH_MASS);
    sprintf(options + strlen(options), "-DNH3_MASS=%f ",               NH3_MASS);
    sprintf(options + strlen(options), "-DH2O_MASS=%f ",               H2O_MASS);
    sprintf(options + strlen(options), "-DCO_MASS=%f ",                CO_MASS);
    sprintf(options + strlen(options), "-DGPU_MASS_PROTEIN_NTERM=%f ", fMassProteinNterm);
    sprintf(options + strlen(options), "-DGPU_MASS_PEPTIDE_NTERM=%f ", fMassPeptideNterm);
    sprintf(options + strlen(options), "-DMAX_PEPTIDE_LENGTH=%d ",     MAX_PEPTIDE_LENGTH);
    sprintf(options + strlen(options), "-DNUM_OUTPUT_PSMS=%d ",        Tempest::params.numInternalPSMs);
    sprintf(options + strlen(options), "-DUSE_A_IONS=%d ",             Tempest::params.useAIons);
    sprintf(options + strlen(options), "-DUSE_B_IONS=%d ",             Tempest::params.useBIons);
    sprintf(options + strlen(options), "-DUSE_C_IONS=%d ",             Tempest::params.useCIons);
    sprintf(options + strlen(options), "-DUSE_X_IONS=%d ",             Tempest::params.useXIons);
    sprintf(options + strlen(options), "-DUSE_Y_IONS=%d ",             Tempest::params.useYIons);
    sprintf(options + strlen(options), "-DUSE_Z_IONS=%d ",             Tempest::params.useZIons);
    sprintf(options + strlen(options), "-DGPU_TOLERANCE=%f ",          Tempest::params.fFragmentTolerance);
    sprintf(options + strlen(options), "-DFRAGMENT_BIN_OFFSET=%f ",    Tempest::params.fragmentBinOffset);
    sprintf(options + strlen(options), "-DGPU_NUM_OUTPUT_PSMS=%d ",    Tempest::params.numInternalPSMs);
    sprintf(options + strlen(options), "-DFLANKING_INTENSITY=%f ",     Tempest::params.flankingIntensity);
    sprintf(options + strlen(options), "-DNUM_AA_NL=%d ",              Tempest::params.numAANL);
    sprintf(options + strlen(options), "-DNUM_NTERM_NL=%d ",           Tempest::params.numNtermNL);
    sprintf(options + strlen(options), "-DNUM_CTERM_NL=%d ",           Tempest::params.numCtermNL);
    sprintf(options + strlen(options), "-DGPU_FIX_DELTA_SCORE=%d ",    Tempest::params.bFixDeltaScore && (Tempest::params.iNumMods > 0));
    sprintf(options + strlen(options), "-DXCORR_TRANSFORM_WIDTH=%d ",  Tempest::params.xcorrTransformWidth);
    //sprintf(options + strlen(options), "-DGPU_NUM_BINS=%d ",           Tempest::tempest.iNumMS2Bins);
    //sprintf(options + strlen(options), "-cl-nv-maxrregcount=16");
    //sprintf(options + strlen(options), "-g ");
    sprintf(options + strlen(options), "-cl-single-precision-constant");
   
    /* Build program */
    err = clBuildProgram(clProgram, 0, NULL, options, NULL, NULL);
    if(PROFILE || err < 0) {

        /* Find size of log and print to std output */
        clGetProgramBuildInfo(clProgram, clDeviceID, CL_PROGRAM_BUILD_LOG, 
                              0, NULL, &log_size);
        program_log = (char*) malloc(log_size + 1);
        program_log[log_size] = '\0';
        clGetProgramBuildInfo(clProgram, clDeviceID, CL_PROGRAM_BUILD_LOG, 
                              log_size + 1, program_log, NULL);
        printf("%s\n", program_log);
        free(program_log);
        if (err < 0)
            exit(1);
    }
    
    create_kernels();
}

void Device::setup(unsigned int minScan, unsigned int maxScan) {
    int err;
    this->minScan = minScan;
    this->maxScan = maxScan;
    
    // Determine Configuration
    this->candidateBufferSize = this->score_reduction_multiple;
    this->score_reduction_size = this->score_reduction_multiple;

    size_t hostMem = sizeof(mObj) * Tempest::tempest.iNumSpectra * Tempest::params.numInternalPSMs
      + sizeof(eObj) * Tempest::tempest.iNumSpectra
      + sizeof(cl_mem) * Tempest::tempest.iNumSpectra
      + sizeof(std::vector<int>) + sizeof(int)*Tempest::host_iPeakBins.size()
      + sizeof(std::vector<float>) + sizeof(float)*Tempest::host_fPeakInts.size()
      + sizeof(int)*Tempest::tempest.iNumSpectra
      + sizeof(long)*Tempest::tempest.iNumSpectra;
    for (int candidateBufferSize=this->score_reduction_multiple; hostMem + candidateBufferSize*Tempest::tempest.iNumSpectra*sizeof(cObj) < Tempest::config.maxHostMem; candidateBufferSize += this->score_reduction_multiple) {
        for (int scoreReductionSize = 1;
             scoreReductionSize <= candidateBufferSize
                 && scoreReductionSize <= this->score_reduction_size_max
                 && scoreReductionSize*(sizeof(int)+sizeof(float)) + this->score_reduction_size_local <= this->lLocalMemSize;
             scoreReductionSize *= 2) {
            if (scoreReductionSize%(this->score_reduction_multiple) == 0 && candidateBufferSize%scoreReductionSize == 0)
                if (candidateBufferSize * scoreReductionSize > this->candidateBufferSize * this->score_reduction_size) {
                    this->candidateBufferSize = candidateBufferSize;
                    this->score_reduction_size = scoreReductionSize;
                }	    
        }
    }
    //printf("Candidate buffer size=%ld\n", this->candidateBufferSize);
    //printf("Reduction size: %ld\n", this->score_reduction_size);


    for (int i=minScan+deviceInd; i<maxScan; i+=Tempest::config.iDevices.size()) {
        eObj* e = Tempest::eScans[i];
        e->candidateBuffer = (cObj*)malloc(this->candidateBufferSize * sizeof(cObj));
        e->candidateBufferSize = this->candidateBufferSize;
        e->clEventSent = clCreateUserEvent(clContext, NULL);
        clSetUserEventStatus(e->clEventSent, 0);
        e->device = this;
    }

    // peaks
    size_t size_iPeakBins = Tempest::tempest.lNumMS2Peaks * sizeof(cl_int);
    size_t size_fPeakInts = Tempest::tempest.lNumMS2Peaks * sizeof(cl_float);
    cl_iPeakBins = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, size_iPeakBins, &(Tempest::host_iPeakBins[0]), &err);
    Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to allocate device memory for peak bins.");
    cl_fPeakInts = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, size_fPeakInts, &(Tempest::host_fPeakInts[0]), &err);
    Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to allocate device memory for peak intensities.");
    
    // cleanup host
    //std::vector<int>().swap(Tempest::host_iPeakBins);
    //std::vector<float>().swap(Tempest::host_fPeakInts);

    //cudaMalloc((void**) &gpu_fSpectra, Tempest::tempest.iNumMS2Bins * sizeof(float));
    //cl_fSpectra = clCreateBuffer(clContext, CL_MEM_READ_WRITE, Tempest::tempest.iNumMS2Bins * sizeof(float), NULL, &err);
    float * init_fSpectra = (float *) calloc(Tempest::tempest.iNumMS2Bins, sizeof(float));
    size_t size_init_fSpectra = Tempest::tempest.iNumMS2Bins * sizeof(cl_float);
    cl_init_fSpectra = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, size_init_fSpectra, init_fSpectra, &err);
    free(init_fSpectra);
    
    // candidate and results
    mObj * init_mPSMs = (mObj *) calloc(Tempest::tempest.iNumSpectra * Tempest::params.numInternalPSMs, sizeof(mObj)); 
    //float * init_fNextScores = (float *) calloc(Tempest::tempest.iNumSpectra, sizeof(float));
    size_t size_cCandidates = sizeof(cObj) * this->candidateBufferSize;
    size_t size_fScores = sizeof(cl_float)  * this->candidateBufferSize;
    size_t size_mPSMs = sizeof(mObj)  * Tempest::tempest.iNumSpectra * Tempest::params.numInternalPSMs;
    //size_t size_fNextScores = sizeof(float) * Tempest::tempest.iNumSpectra;
    cl_cCandidates = clCreateBuffer(clContext, CL_MEM_READ_ONLY, size_cCandidates, NULL, &err);
    cl_fScores = clCreateBuffer(clContext, CL_MEM_READ_WRITE, size_fScores, NULL, &err);  
    cl_mPSMs = clCreateBuffer(clContext, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, size_mPSMs , init_mPSMs, &err);
    //cl_fNextScores = clCreateBuffer(clContext, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, size_fNextScores, init_fNextScores, &err);
    //MEA: need to block free until previous clCreateBuffer commands complete?
    free(init_mPSMs);
    //free(init_fNextScores);
    Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to allocate device memory for candidates and results.");

    //determine how many spectra can be kept in device memory at a time
    size_t availMemSpectra = lGlobalMemSize
        - size_iPeakBins
        - size_fPeakInts
        - size_init_fSpectra
        - size_cCandidates
        - size_fScores
        - size_mPSMs;
    if (availMemSpectra > Tempest::config.maxDeviceMem)
        availMemSpectra = Tempest::config.maxDeviceMem;
    long maxCachedSpectra = availMemSpectra / (Tempest::tempest.iNumMS2Bins*sizeof(cl_float));
    if (maxCachedSpectra > (long)ceil(float(Tempest::tempest.iNumSpectra)/Tempest::devices.size()))
        maxCachedSpectra = (long)ceil(float(Tempest::tempest.iNumSpectra)/Tempest::devices.size());
    if (maxCachedSpectra <= 0)
        maxCachedSpectra = 1;
    
    printf(" » (%d:%d) Allocating %.2f MB of device memory for %ld cached %s.\n", platformID, deviceID, (float)maxCachedSpectra*Tempest::tempest.iNumMS2Bins*sizeof(cl_float)/MB, maxCachedSpectra, maxCachedSpectra==1 ? "spectrum" : "spectra");
    for (int i=0; i<maxCachedSpectra; i++) {
        cl_mem newBuffer = clCreateBuffer(clContext, CL_MEM_READ_WRITE, Tempest::tempest.iNumMS2Bins*sizeof(cl_float), NULL, &err);
        Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to allocate spectrum memory on device.");
        unusedBuffers.push(newBuffer);
    }
    
    setup_constant_memory();
    
    //initialize profiling variables
    scoreEvent = clCreateUserEvent(clContext, NULL);
    reduceEvent = clCreateUserEvent(clContext, NULL);
    buildEvent = clCreateUserEvent(clContext, NULL);
    memsetEvent = clCreateUserEvent(clContext, NULL);
    transformEvent = clCreateUserEvent(clContext, NULL);
    totalScoreTime = 0;
    totalReduceTime = 0;
    totalBuildTime = 0;
    totalTransformTime = 0;
    totalMemsetTime = 0;
    totalSendTime = 0;
    buildLaunches = 0;
    scoreKernelLaunches = 0;
    lastBuildIndex = -1;

}

void Device::scoreCandidates(eObj *e) {
    //e->iNumBufferedCandidates = 0;
    //return;
    //MEA: static?
    static cObj* p;
    //static size_t iNumBlocks;
    static size_t stGlobalDim;
    static size_t globalTransDim = Tempest::mround(Tempest::tempest.iNumMS2Bins, this->transform_size);
    static float fElapsedTime;
    long lSpectrumOffset = e->lIndex*Tempest::tempest.iNumMS2Bins;
    long lScratchOffset = (long)Tempest::tempest.iCrossCorrelationWidth;
    long lNoOffset = 0;
    int err;
    cl_ulong start;
    cl_ulong end;
    
    err = clEnqueueWriteBuffer(clCommandQueue, cl_cCandidates, CL_FALSE, 0, sizeof(cObj) * e->iNumBufferedCandidates, e->candidateBuffer, 0, NULL, &(e->clEventSent));
    Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to copy candidate data from host to GPU");
	
    stGlobalDim = Tempest::mround(Tempest::host_iPeakCounts[e->lIndex], this->build_size);
    cl_mem spectrumBuffer;

    std::map<long,cl_mem>::iterator s2bElem = spectrum2buffer.find(e->lIndex);
    if (s2bElem == spectrum2buffer.end()) { //spectrum not cached
        if (!unusedBuffers.empty()) {
            spectrumBuffer = unusedBuffers.top();
            unusedBuffers.pop();
        }
        else {
            spectrumBuffer = spectrum2buffer.begin()->second;
            spectrum2buffer.erase(spectrum2buffer.begin());
        }
        spectrum2buffer[e->lIndex] = spectrumBuffer;

        //initialize buffer
        err = clEnqueueCopyBuffer(clCommandQueue, cl_init_fSpectra, spectrumBuffer, 0, 0, Tempest::tempest.iNumMS2Bins*sizeof(cl_float), 0, NULL, PROFILE ? &memsetEvent : NULL);
        //Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to clear spectrum memory");
        if (err != 0) {
            //memory cap reached. Stop filling new buffers.
            unusedBuffers = std::stack<cl_mem>();
            spectrumBuffer = spectrum2buffer.begin()->second;
            spectrum2buffer.erase(spectrum2buffer.begin());
            spectrum2buffer[e->lIndex] = spectrumBuffer;
            err = clEnqueueCopyBuffer(clCommandQueue, cl_init_fSpectra, spectrumBuffer, 0, 0, Tempest::tempest.iNumMS2Bins*sizeof(cl_float), 0, NULL, PROFILE ? &memsetEvent : NULL);
            Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to clear spectrum memory");
        }
        if (PROFILE) {
            clFinish(clCommandQueue);
            clGetEventProfilingInfo(memsetEvent, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, NULL);
            clGetEventProfilingInfo(memsetEvent, CL_PROFILING_COMMAND_END,   sizeof(cl_ulong), &end,   NULL);
            totalMemsetTime += (end-start);
            clReleaseEvent(memsetEvent);
        }
        
        // build
        err  = clSetKernelArg(__gpu_build, 0, sizeof(cl_mem), &spectrumBuffer);
        err |= clSetKernelArg(__gpu_build, 2, sizeof(int), &(Tempest::host_iPeakCounts[e->lIndex]));
        err |= clSetKernelArg(__gpu_build, 5, sizeof(long), &(Tempest::host_lPeakIndices[e->lIndex]));
        err |= clEnqueueNDRangeKernel(clCommandQueue, __gpu_build, 1, NULL, &stGlobalDim, &(this->build_size), 0, NULL, PROFILE ? &buildEvent : NULL);
        Tempest::check_cl_error(__FILE__, __LINE__, err, "Could not build spectrum (gpu_build kernel)");
        if (PROFILE) {
            clFinish(clCommandQueue);
            clGetEventProfilingInfo(buildEvent, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, NULL);
            clGetEventProfilingInfo(buildEvent, CL_PROFILING_COMMAND_END,   sizeof(cl_ulong), &end,   NULL);
            totalBuildTime += (end-start);
            buildLaunches += 1;
            clReleaseEvent(buildEvent);
        }

        // transform
        if (Tempest::params.xcorrTransformWidth) {
            //size_t localDim = CROSS_CORRELATION_WINDOW * 2;
            //size_t globalDim = localDim * Tempest::tempest.iNumMS2Bins;
            size_t globalDim = Tempest::mround(Tempest::tempest.iNumMS2Bins, this->transform_size);
            err  = clSetKernelArg(__gpu_transform, 0, sizeof(cl_mem), &spectrumBuffer);
            err |= clEnqueueNDRangeKernel(clCommandQueue, __gpu_transform, 1, NULL, &globalDim, &(this->transform_size), 0, NULL, PROFILE ? & transformEvent : NULL);
            Tempest::check_cl_error(__FILE__, __LINE__, err, "Could not transform spectrum (gpu_transform kernel)");
            if (PROFILE) {
                clFinish(clCommandQueue);
                clGetEventProfilingInfo(transformEvent, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, NULL);
                clGetEventProfilingInfo(transformEvent, CL_PROFILING_COMMAND_END,   sizeof(cl_ulong), &end,   NULL);
                totalTransformTime += (end-start);
                clReleaseEvent(transformEvent);
            }
        }
    }
    else {
        //move spectrum entry to end of map by reinserting
        spectrumBuffer = s2bElem->second;
        spectrum2buffer.erase(s2bElem);
        spectrum2buffer[e->lIndex] = spectrumBuffer;
    }
        
    // score
    err  = clSetKernelArg(__gpu_score, 0, sizeof(int), &(e->iPrecursorCharge));
    err |= clSetKernelArg(__gpu_score, 1, sizeof(int), &(e->iNumBufferedCandidates));
    err |= clSetKernelArg(__gpu_score, 4, sizeof(cl_mem), &spectrumBuffer);
    err |= clSetKernelArg(__gpu_score, 5, sizeof(long), &lNoOffset);
    err |= clEnqueueNDRangeKernel(clCommandQueue, __gpu_score, 1, NULL, &(this->candidateBufferSize), &(this->score_size), 0, NULL, PROFILE ? &scoreEvent : NULL);
    Tempest::check_cl_error(__FILE__, __LINE__, err, "Could not score candidates (gpu_score kernel)");
    if (PROFILE) {
        clFinish(clCommandQueue);
        clGetEventProfilingInfo(scoreEvent, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, NULL);
        clGetEventProfilingInfo(scoreEvent, CL_PROFILING_COMMAND_END,   sizeof(cl_ulong), &end,   NULL);
        totalScoreTime += (end-start);
        clReleaseEvent(scoreEvent);
        scoreKernelLaunches++;
    }
    
    // Process Scores
	
    // TODO what if buffer size is less than 512?
    long lPSMsOffset = e->lIndex * Tempest::params.numInternalPSMs;
    err |= clSetKernelArg(__gpu_score_reduction, 4, sizeof(long), &lPSMsOffset);
    err |= clSetKernelArg(__gpu_score_reduction, 5, sizeof(long), &(e->lIndex));
    err |= clEnqueueNDRangeKernel(clCommandQueue, __gpu_score_reduction, 1, NULL, &(this->score_reduction_size), &(this->score_reduction_size), 0, NULL, PROFILE ? &reduceEvent : NULL);
    //err |= clEnqueueTask(clCommandQueue, __gpu_score_reduction, 0, NULL, PROFILE ? &reduceEvent : NULL);
    Tempest::check_cl_error(__FILE__, __LINE__, err, "Could not process scores (gpu_score_reduction kernel)");
    if (PROFILE) {
        clFinish(clCommandQueue);
        clGetEventProfilingInfo(reduceEvent, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, NULL);
        clGetEventProfilingInfo(reduceEvent, CL_PROFILING_COMMAND_END,   sizeof(cl_ulong), &end,   NULL);
        totalReduceTime += (end-start);
        clReleaseEvent(reduceEvent);
    }
	
    // Stats
    Tempest::gpu_info.iNumScoringKernels += 1;

    // reset buffer
    e->iNumBufferedCandidates = 0;
}

void Device::finish() {
    clFinish(clCommandQueue);
}

int Device::get_mPSMs(mObj* destination) {
    mObj* temp_mPSMs = (mObj*)malloc(sizeof(mObj)*Tempest::tempest.iNumSpectra*Tempest::params.numInternalPSMs);
    int err;
    err = clEnqueueReadBuffer(clCommandQueue, cl_mPSMs, CL_TRUE, 0, sizeof(mObj)*Tempest::tempest.iNumSpectra*Tempest::params.numInternalPSMs, temp_mPSMs, 0, NULL, NULL);
    if (err != 0)
        return err;
    //read interleaved results from this device (possibility of multiple devices)
    for (int i=this->minScan+this->deviceInd; i<this->maxScan; i+=Tempest::config.iDevices.size())
        for (int j=0; j<Tempest::params.numInternalPSMs; j++)
            destination[i*Tempest::params.numInternalPSMs + j] = temp_mPSMs[i*Tempest::params.numInternalPSMs + j];
    free(temp_mPSMs);
    return err;       
}

// int Device::get_fNextScores(float* destination) {
//     return clEnqueueReadBuffer(clCommandQueue, cl_fNextScores, CL_TRUE, 0, sizeof(float)*Tempest::tempest.iNumSpectra, destination, 0, NULL, NULL);
// }

void Device::printProfilingData() {
    printf("Total memset time:       %fs\n", (float)totalMemsetTime/1000000000);
    printf("Total build time:        %fs\n", (float)totalBuildTime/1000000000);
    printf("Total transform time:    %fs\n", (float)totalTransformTime/1000000000);
    printf("Total send time:         %fs\n", (float)totalSendTime/1000000000);
    printf("Total score time:        %fs\n", (float)totalScoreTime/1000000000);
    printf("Total reduce time:       %fs\n", (float)totalReduceTime/1000000000);
    printf("Build launches:          %ld\n", buildLaunches);
    printf("Scoring kernel launches: %ld\n", scoreKernelLaunches);
}

/*
 * Cleanup device memory and events. Prints any timing requests.

 extern void cleanup_device()
 {
 int err;
 err  = clReleaseMemObject(cl_lPeakIndices);
 err |= clReleaseMemObject(cl_iPeakCounts);
 err |= clReleaseMemObject(cl_iPeakBins);
 err |= clReleaseMemObject(cl_fPeakInts);
 err |= clReleaseMemObject(cl_fSpectra);
 err |= clReleaseMemObject(cl_cCandidates);
 err |= clReleaseMemObject(cl_fScores);
 err |= clReleaseMemObject(cl_mPSMs);
 err |= clReleaseMemObject(cl_fNextScores);
 Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to free memory object(s)");
 }
*/
