#include "tempest.hpp"
#include "kernels.h"

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
        fMassAA[iResidue] = (float) Tempest::params.dMassAA[iResidue];
    }

    cl_mem CL_MASS_AA      = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, 256*sizeof(float), &fMassAA, &err);
    cl_mem UN_MOD_AA        = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, 256*sizeof(char),  &Tempest::params.unModAA, &err);
    cl_mem NTERM_MOD_MASSES = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, 5*sizeof(float),  &Tempest::params.ntermModMasses, &err);
    cl_mem CTERM_MOD_MASSES = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, 5*sizeof(float),  &Tempest::params.ctermModMasses, &err);
    
    //neutral losses
    //transfer to device
    if (Tempest::params.numNtermNL)
        cl_nlValuesNterm = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(nlValue)*Tempest::params.nlValuesNterm.size(), &Tempest::params.nlValuesNterm[0], &err);
    if (Tempest::params.numCtermNL)
        cl_nlValuesCterm = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(nlValue)*Tempest::params.nlValuesCterm.size(), &Tempest::params.nlValuesCterm[0], &err);
    if (Tempest::params.numAANL)
        cl_nlValuesAA = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(nlValue)*Tempest::params.nlValuesAA.size(), &Tempest::params.nlValuesAA[0], &err);
    Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to copy device constant memory");
    
    //pre-set these constant arrays as arguments to scoring kernels (re-used in every launch)
    err  = clSetKernelArg(__cl_build, 2, sizeof(cl_mem), &cl_iPeakBins);
    err |= clSetKernelArg(__cl_build, 3, sizeof(cl_mem), &cl_fPeakInts);
    Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to set constant args (__cl_build)");

    err = clSetKernelArg(__cl_transform, 1, sizeof(int),  &Tempest::data.iNumMS2Bins);
    Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to set constant args (__cl_transform)");

    err  = clSetKernelArg(__cl_score, 2, sizeof(cl_mem), &cl_cCandidates);
    err |= clSetKernelArg(__cl_score, 3, sizeof(cl_mem), &cl_fScores);
    err |= clSetKernelArg(__cl_score, 6, sizeof(int),    &Tempest::data.iNumMS2Bins);
    err |= clSetKernelArg(__cl_score, 7, sizeof(cl_mem), &CL_MASS_AA);
    err |= clSetKernelArg(__cl_score, 8, sizeof(cl_mem), &NTERM_MOD_MASSES);
    err |= clSetKernelArg(__cl_score, 9, sizeof(cl_mem), &CTERM_MOD_MASSES);
    err |= clSetKernelArg(__cl_score, 10, sizeof(cl_mem), Tempest::params.numNtermNL ? &cl_nlValuesNterm : NULL);
    err |= clSetKernelArg(__cl_score, 11, sizeof(cl_mem), Tempest::params.numCtermNL ? &cl_nlValuesCterm : NULL);
    err |= clSetKernelArg(__cl_score, 12, sizeof(cl_mem), Tempest::params.numAANL ? &cl_nlValuesAA : NULL);
    Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to set constant args (__cl_score)");

    size_t matchesLocalSize = this->reduce_scores_size * sizeof(int);
    size_t scoresLocalSize = this->reduce_scores_size * sizeof(float);
    err  = clSetKernelArg(__cl_reduce_scores, 0, sizeof(int), &(this->candidateBufferSize));
    err |= clSetKernelArg(__cl_reduce_scores, 1, sizeof(cl_mem), &cl_cCandidates);
    err |= clSetKernelArg(__cl_reduce_scores, 2, sizeof(cl_mem), &cl_fScores);
    err |= clSetKernelArg(__cl_reduce_scores, 3, sizeof(cl_mem), &cl_mPSMs);
    if (Tempest::config.parallelReduce) {
        //err |= clSetKernelArg(__cl_reduce_scores, 5, zMemShared, NULL);
        err |= clSetKernelArg(__cl_reduce_scores, 5, this->reduce_scores_size * sizeof(int), NULL);
        err |= clSetKernelArg(__cl_reduce_scores, 6, this->reduce_scores_size * sizeof(float), NULL);
    }
    Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to set constant args (__cl_reduce_scores)");
    
}

void Device::create_kernels() {
    int err;
    size_t paramBuffer;
    
    this->__cl_build = clCreateKernel(clProgram, "cl_build", &err);
    Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to create kernel");
    //clGetKernelWorkGroupInfo(__cl_build, this->clDeviceID, CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, sizeof(size_t), &(this->build_size), NULL);
    clGetKernelWorkGroupInfo(__cl_build, this->clDeviceID, CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, sizeof(size_t), &(this->build_size), NULL);
    if (this->build_size < Tempest::config.minWorkSize)
        this->build_size = Tempest::config.minWorkSize;
    
    this->__cl_transform = clCreateKernel(clProgram, "cl_transform", &err);
    Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to create kernel");
    clGetKernelWorkGroupInfo(__cl_transform, this->clDeviceID, CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, sizeof(size_t), &(this->transform_size), NULL);
    if (this->transform_size < Tempest::config.minWorkSize)
        this->transform_size = Tempest::config.minWorkSize;   
    
    this->__cl_score = clCreateKernel(clProgram, "cl_score", &err);
    Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to create kernel");
    clGetKernelWorkGroupInfo(__cl_score, this->clDeviceID, CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, sizeof(size_t), &(this->score_size), NULL);
    if (this->score_size < Tempest::config.minWorkSize)
      this->score_size = Tempest::config.minWorkSize;  
    
    this->__cl_reduce_scores = clCreateKernel(clProgram, Tempest::config.parallelReduce ? "cl_reduce_scores_parallel" : "cl_reduce_scores_sequential", &err);
    Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to create kernel");
    clGetKernelWorkGroupInfo(__cl_reduce_scores, this->clDeviceID, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &(this->reduce_scores_size_max), NULL);
    clGetKernelWorkGroupInfo(__cl_reduce_scores, this->clDeviceID, CL_KERNEL_LOCAL_MEM_SIZE, sizeof(size_t), &(this->reduce_scores_size_local), NULL);
    clGetKernelWorkGroupInfo(__cl_reduce_scores, this->clDeviceID, CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, sizeof(size_t), &this->reduce_scores_multiple, NULL);
     
    // this->__cl_memset = clCreateKernel(clProgram, "cl_memset", &err);
    // Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to create kernel");
    // clGetKernelWorkGroupInfo(__cl_memset, this->clDeviceID, CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, sizeof(size_t), &(this->memset_size), NULL);
    // printf("cl_memset: PREFERRED_WORK_GROUP_SIZE_MULTIPLE=%ld\n", paramBuffer);
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
    clCommandQueue = clCreateCommandQueue(clContext, clDeviceID, Tempest::config.profile ? CL_QUEUE_PROFILING_ENABLE : 0, &err);
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
    /*
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
    */

    // Create program from file
    //clProgram = clCreateProgramWithSource(clContext, 1, (const char**)&program_buffer, &program_size, &err);

    // Create program from embedded source
    const char* kernels_ptr = (const char*)&(kernels_cl[0]);
    size_t kernels_size = (size_t)kernels_cl_len;
    clProgram = clCreateProgramWithSource(clContext, 1, &kernels_ptr, &kernels_size, &err);
    if(err < 0) {
        perror("Couldn't create the kernels program");
        exit(1);
    }
    //free(program_buffer);

    // Define known constants as literals to build into the program
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
    sprintf(options + strlen(options), "-DMASS_PROTEIN_NTERM=%f ",     fMassProteinNterm);
    sprintf(options + strlen(options), "-DMASS_PEPTIDE_NTERM=%f ",     fMassPeptideNterm);
    sprintf(options + strlen(options), "-DMAX_PEPTIDE_LENGTH=%d ",     MAX_PEPTIDE_LENGTH);
    sprintf(options + strlen(options), "-DNUM_INTERNAL_PSMS=%d ",      Tempest::params.numInternalPSMs);
    sprintf(options + strlen(options), "-DUSE_A_IONS=%d ",             Tempest::params.useAIons);
    sprintf(options + strlen(options), "-DUSE_B_IONS=%d ",             Tempest::params.useBIons);
    sprintf(options + strlen(options), "-DUSE_C_IONS=%d ",             Tempest::params.useCIons);
    sprintf(options + strlen(options), "-DUSE_X_IONS=%d ",             Tempest::params.useXIons);
    sprintf(options + strlen(options), "-DUSE_Y_IONS=%d ",             Tempest::params.useYIons);
    sprintf(options + strlen(options), "-DUSE_Z_IONS=%d ",             Tempest::params.useZIons);
    sprintf(options + strlen(options), "-DFRAGMENT_TOLERANCE=%f ",     Tempest::params.fFragmentTolerance);
    sprintf(options + strlen(options), "-DFRAGMENT_BIN_OFFSET=%f ",    Tempest::params.fragmentBinOffset);
    sprintf(options + strlen(options), "-DNUM_AA_NL=%d ",              Tempest::params.numAANL);
    sprintf(options + strlen(options), "-DNUM_NTERM_NL=%d ",           Tempest::params.numNtermNL);
    sprintf(options + strlen(options), "-DNUM_CTERM_NL=%d ",           Tempest::params.numCtermNL);
    sprintf(options + strlen(options), "-DXCORR_TRANSFORM_WIDTH=%d ",  Tempest::params.xcorrTransformWidth);
    sprintf(options + strlen(options), "-cl-single-precision-constant");
   
    /* Build program */
    err = clBuildProgram(clProgram, 0, NULL, options, NULL, NULL);
    if(Tempest::config.profile || err < 0) {

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
    this->candidateBufferSize = this->reduce_scores_multiple;
    this->reduce_scores_size = this->reduce_scores_multiple;

    size_t hostMem = sizeof(mObj) * Tempest::data.iNumSpectra * Tempest::params.numInternalPSMs
      + sizeof(eObj) * Tempest::data.iNumSpectra
      + sizeof(cl_mem) * Tempest::data.iNumSpectra
      + sizeof(std::vector<int>) + sizeof(int)*Tempest::data.host_iPeakBins.size()
      + sizeof(std::vector<float>) + sizeof(float)*Tempest::data.host_fPeakInts.size()
      + sizeof(int)*Tempest::data.iNumSpectra
      + sizeof(long)*Tempest::data.iNumSpectra;
    for (int candidateBufferSize=this->reduce_scores_multiple; hostMem + candidateBufferSize*Tempest::data.iNumSpectra*sizeof(cObj) < Tempest::config.maxHostMem; candidateBufferSize += this->reduce_scores_multiple) {
        for (int reduceScoresSize = 1;
             reduceScoresSize <= candidateBufferSize
                 && reduceScoresSize <= this->reduce_scores_size_max
                 && reduceScoresSize*(sizeof(int) + sizeof(float)) + this->reduce_scores_size_local <= this->lLocalMemSize;
             reduceScoresSize *= 2) {
            if (reduceScoresSize%(this->reduce_scores_multiple) == 0 && candidateBufferSize%reduceScoresSize == 0)
                if (candidateBufferSize * reduceScoresSize > this->candidateBufferSize * this->reduce_scores_size) {
                    this->candidateBufferSize = candidateBufferSize;
                    this->reduce_scores_size = reduceScoresSize;
                }	    
        }
    }
    if (Tempest::config.profile) {
        printf("cl_build: local_work_size=%ld\n", this->build_size);
        printf("cl_transform: local_work_size=%ld\n", this->transform_size);
        printf("cl_score: local_work_size=%ld\n", this->score_size);
        printf("candidate buffer size=%ld\n", this->candidateBufferSize);
        printf("cl_reduce_scores: local_work_size=%ld\n", this->reduce_scores_size);
    }

    for (int i=minScan+deviceInd; i<maxScan; i+=Tempest::config.iDevices.size()) {
        eObj* e = Tempest::data.eScans[i];
        e->candidateBuffer = (cObj*)malloc(this->candidateBufferSize * sizeof(cObj));
        e->candidateBufferSize = this->candidateBufferSize;
        e->clEventSent = clCreateUserEvent(clContext, NULL);
        clSetUserEventStatus(e->clEventSent, 0);
        e->device = this;
    }

    // peaks
    size_t size_iPeakBins = Tempest::data.lNumMS2Peaks * sizeof(cl_int);
    size_t size_fPeakInts = Tempest::data.lNumMS2Peaks * sizeof(cl_float);
    cl_iPeakBins = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, size_iPeakBins, &(Tempest::data.host_iPeakBins[0]), &err);
    Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to allocate device memory for peak bins.");
    cl_fPeakInts = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, size_fPeakInts, &(Tempest::data.host_fPeakInts[0]), &err);
    Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to allocate device memory for peak intensities.");
    
    // cleanup host
    //std::vector<int>().swap(Tempest::data.host_iPeakBins);
    //std::vector<float>().swap(Tempest::data.host_fPeakInts);

    //cudaMalloc((void**) &cl_fSpectra, Tempest::data.iNumMS2Bins * sizeof(float));
    //cl_fSpectra = clCreateBuffer(clContext, CL_MEM_READ_WRITE, Tempest::data.iNumMS2Bins * sizeof(float), NULL, &err);
    float * init_fSpectra = (float *) calloc(Tempest::data.iNumMS2Bins, sizeof(float));
    size_t size_init_fSpectra = Tempest::data.iNumMS2Bins * sizeof(cl_float);
    cl_init_fSpectra = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, size_init_fSpectra, init_fSpectra, &err);
    free(init_fSpectra);
    
    // candidate and results
    mObj * init_mPSMs = (mObj *) calloc(Tempest::data.iNumSpectra * Tempest::params.numInternalPSMs, sizeof(mObj));
    // for (int i=0; i<Tempest::data.iNumSpectra * Tempest::params.numInternalPSMs; i++)
    //     init_mPSMs[i].fScore = MIN_SCORE;
    //float * init_fNextScores = (float *) calloc(Tempest::data.iNumSpectra, sizeof(float));
    size_t size_cCandidates = sizeof(cObj) * this->candidateBufferSize;
    size_t size_fScores = sizeof(cl_float)  * this->candidateBufferSize;
    size_t size_mPSMs = sizeof(mObj)  * Tempest::data.iNumSpectra * Tempest::params.numInternalPSMs;
    //size_t size_fNextScores = sizeof(float) * Tempest::data.iNumSpectra;
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
    long maxCachedSpectra = availMemSpectra / (Tempest::data.iNumMS2Bins*sizeof(cl_float));
    if (maxCachedSpectra > (long)ceil(float(Tempest::data.iNumSpectra)/Tempest::devices.size()))
        maxCachedSpectra = (long)ceil(float(Tempest::data.iNumSpectra)/Tempest::devices.size());
    if (maxCachedSpectra <= 0)
        maxCachedSpectra = 1;
    
    printf(" » (%d:%d) Allocating %.2f MB of device memory for %ld cached %s.\n", platformID, deviceID, (float)maxCachedSpectra*Tempest::data.iNumMS2Bins*sizeof(cl_float)/MB, maxCachedSpectra, maxCachedSpectra==1 ? "spectrum" : "spectra");
    for (int i=0; i<maxCachedSpectra; i++) {
        cl_mem newBuffer = clCreateBuffer(clContext, CL_MEM_READ_WRITE, Tempest::data.iNumMS2Bins*sizeof(cl_float), NULL, &err);
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
    static size_t globalTransDim = Tempest::mround(Tempest::data.iNumMS2Bins, this->transform_size);
    static float fElapsedTime;
    long lSpectrumOffset = e->lIndex*Tempest::data.iNumMS2Bins;
    long lScratchOffset = (long)Tempest::data.iCrossCorrelationWidth;
    long lNoOffset = 0;
    int err;
    cl_ulong start;
    cl_ulong end;
    
    err = clEnqueueWriteBuffer(clCommandQueue, cl_cCandidates, CL_FALSE, 0, sizeof(cObj) * e->iNumBufferedCandidates, e->candidateBuffer, 0, NULL, &(e->clEventSent));
    Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to copy candidate data from host to GPU");
	
    stGlobalDim = Tempest::mround(Tempest::data.host_iPeakCounts[e->lIndex], this->build_size);
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
        err = clEnqueueCopyBuffer(clCommandQueue, cl_init_fSpectra, spectrumBuffer, 0, 0, Tempest::data.iNumMS2Bins*sizeof(cl_float), 0, NULL, Tempest::config.profile ? &memsetEvent : NULL);
        //Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to clear spectrum memory");
        if (err != 0) {
            //memory cap reached. Stop filling new buffers.
            unusedBuffers = std::stack<cl_mem>();
            spectrumBuffer = spectrum2buffer.begin()->second;
            spectrum2buffer.erase(spectrum2buffer.begin());
            spectrum2buffer[e->lIndex] = spectrumBuffer;
            err = clEnqueueCopyBuffer(clCommandQueue, cl_init_fSpectra, spectrumBuffer, 0, 0, Tempest::data.iNumMS2Bins*sizeof(cl_float), 0, NULL, Tempest::config.profile ? &memsetEvent : NULL);
            Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to clear spectrum memory");
        }
        if (Tempest::config.profile) {
            clFinish(clCommandQueue);
            clGetEventProfilingInfo(memsetEvent, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, NULL);
            clGetEventProfilingInfo(memsetEvent, CL_PROFILING_COMMAND_END,   sizeof(cl_ulong), &end,   NULL);
            totalMemsetTime += (end-start);
            clReleaseEvent(memsetEvent);
        }
        
        // build
        err  = clSetKernelArg(__cl_build, 0, sizeof(cl_mem), &spectrumBuffer);
        err |= clSetKernelArg(__cl_build, 1, sizeof(int), &(Tempest::data.host_iPeakCounts[e->lIndex]));
        err |= clSetKernelArg(__cl_build, 4, sizeof(long), &(Tempest::data.host_lPeakIndices[e->lIndex]));
        err |= clEnqueueNDRangeKernel(clCommandQueue, __cl_build, 1, NULL, &stGlobalDim, &(this->build_size), 0, NULL, Tempest::config.profile ? &buildEvent : NULL);
        Tempest::check_cl_error(__FILE__, __LINE__, err, "Could not build spectrum (cl_build kernel)");
        if (Tempest::config.profile) {
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
            //size_t globalDim = localDim * Tempest::data.iNumMS2Bins;
            size_t globalDim = Tempest::mround(Tempest::data.iNumMS2Bins, this->transform_size);
            err  = clSetKernelArg(__cl_transform, 0, sizeof(cl_mem), &spectrumBuffer);
            err |= clEnqueueNDRangeKernel(clCommandQueue, __cl_transform, 1, NULL, &globalDim, &(this->transform_size), 0, NULL, Tempest::config.profile ? & transformEvent : NULL);
            Tempest::check_cl_error(__FILE__, __LINE__, err, "Could not transform spectrum (cl_transform kernel)");
            if (Tempest::config.profile) {
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
    err  = clSetKernelArg(__cl_score, 0, sizeof(int), &(e->iPrecursorCharge));
    err |= clSetKernelArg(__cl_score, 1, sizeof(int), &(e->iNumBufferedCandidates));
    err |= clSetKernelArg(__cl_score, 4, sizeof(cl_mem), &spectrumBuffer);
    err |= clSetKernelArg(__cl_score, 5, sizeof(long), &lNoOffset);
    err |= clEnqueueNDRangeKernel(clCommandQueue, __cl_score, 1, NULL, &(this->candidateBufferSize), &(this->score_size), 0, NULL, Tempest::config.profile ? &scoreEvent : NULL);
    Tempest::check_cl_error(__FILE__, __LINE__, err, "Could not score candidates (cl_score kernel)");
    if (Tempest::config.profile) {
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
    err |= clSetKernelArg(__cl_reduce_scores, 4, sizeof(long), &lPSMsOffset);
    if (Tempest::config.parallelReduce)
        err |= clEnqueueNDRangeKernel(clCommandQueue, __cl_reduce_scores, 1, NULL, &(this->reduce_scores_size), &(this->reduce_scores_size), 0, NULL, Tempest::config.profile ? &reduceEvent : NULL);
    else
        err |= clEnqueueTask(clCommandQueue, __cl_reduce_scores, 0, NULL, Tempest::config.profile ? &reduceEvent : NULL);
    Tempest::check_cl_error(__FILE__, __LINE__, err, "Could not process scores (cl_reduce_scores kernel)");
    if (Tempest::config.profile) {
        clFinish(clCommandQueue);
        clGetEventProfilingInfo(reduceEvent, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, NULL);
        clGetEventProfilingInfo(reduceEvent, CL_PROFILING_COMMAND_END,   sizeof(cl_ulong), &end,   NULL);
        totalReduceTime += (end-start);
        clReleaseEvent(reduceEvent);
    }

    // reset buffer
    e->iNumBufferedCandidates = 0;
}

void Device::finish() {
    clFinish(clCommandQueue);
}

int Device::get_mPSMs(mObj* destination) {
    mObj* temp_mPSMs = (mObj*)malloc(sizeof(mObj)*Tempest::data.iNumSpectra*Tempest::params.numInternalPSMs);
    int err;
    err = clEnqueueReadBuffer(clCommandQueue, cl_mPSMs, CL_TRUE, 0, sizeof(mObj)*Tempest::data.iNumSpectra*Tempest::params.numInternalPSMs, temp_mPSMs, 0, NULL, NULL);
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
//     return clEnqueueReadBuffer(clCommandQueue, cl_fNextScores, CL_TRUE, 0, sizeof(float)*Tempest::data.iNumSpectra, destination, 0, NULL, NULL);
// }

void Device::printProfilingData() {
    printf("(%d:%d)\n", this->platformID, this->deviceID);
    printf("Total memset time:       %fs\n", (float)totalMemsetTime/1000000000);
    printf("Total build time:        %fs\n", (float)totalBuildTime/1000000000);
    printf("Total transform time:    %fs\n", (float)totalTransformTime/1000000000);
    //printf("Total send time:         %fs\n", (float)totalSendTime/1000000000);
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
