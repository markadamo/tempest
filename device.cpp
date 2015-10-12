#include "tempest.hpp"

/*
 * Launches a kernel to perform parallel memset on a block of device memory
 */
int Device::cl_memset(cl_mem buffer, int c, unsigned long n, cl_event* evt) {
    int err;
    size_t global_dim = mround(n, MAX_BLOCKDIM);
    err  = clSetKernelArg(__cl_memset, 0, sizeof(cl_mem), &buffer);
    err |= clSetKernelArg(__cl_memset, 1, sizeof(int), &c);
    err |= clSetKernelArg(__cl_memset, 2, sizeof(long), &n);
    err |= clEnqueueNDRangeKernel(clCommandQueue, __cl_memset, 1, NULL, &global_dim, &MAX_BLOCKDIM, 0, NULL, evt);
    return err;
}

void Device::setup_constant_memory() {
    int iResidue;
    float fMassAA['z'-'A'+1];
    int err;

    // convert the double masses to floats
    for (iResidue=0; iResidue < 'z'-'A'+1; iResidue++) {
        fMassAA[iResidue] = (float) dMassAA[iResidue];
    }

    GPU_MASS_AA            = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float) * ('z'-'A'+1), &fMassAA,                        &err);
    GPU_MASS_PROTON        = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float) * MAX_CHARGE,  &fMassProton,                    &err);
    GPU_MASS_NH3           = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float) * MAX_CHARGE,  &fMassNH3,                       &err);
    GPU_MASS_H2O           = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float) * MAX_CHARGE,  &fMassH2O,                       &err);
    GPU_MASS_CO            = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float) * MAX_CHARGE,  &fMassCO,                        &err);

    //neutral losses
    //transfer to device
    if (params.numNtermNL)
        cl_nlValuesNterm = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(nlValue)*nlValuesNterm.size(), &nlValuesNterm[0], &err);
    if (params.numCtermNL)
        cl_nlValuesCterm = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(nlValue)*nlValuesCterm.size(), &nlValuesCterm[0], &err);
    if (params.numAANL)
        cl_nlValuesAA = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(nlValue)*nlValuesAA.size(), &nlValuesAA[0], &err);
    check_cl_error(__FILE__, __LINE__, err, "Unable to copy device constant memory");
    
    //pre-set these constant arrays as arguments to scoring kernels (re-used in every launch)
    err  = clSetKernelArg(__gpu_build, 3, sizeof(cl_mem), &cl_iPeakBins);
    err |= clSetKernelArg(__gpu_build, 4, sizeof(cl_mem), &cl_fPeakInts);
    check_cl_error(__FILE__, __LINE__, err, "Unable to set constant args (__gpu_build)");

    err  = clSetKernelArg(__gpu_score, 2, sizeof(cl_mem), &cl_cCandidates);
    err |= clSetKernelArg(__gpu_score, 3, sizeof(cl_mem), &cl_fScores);
    err |= clSetKernelArg(__gpu_score, 6, sizeof(cl_mem), &GPU_MASS_AA);
    err |= clSetKernelArg(__gpu_score, 7, sizeof(cl_mem), &GPU_MASS_PROTON);
    err |= clSetKernelArg(__gpu_score, 8, sizeof(cl_mem), params.numNtermNL ? &cl_nlValuesNterm : NULL);
    err |= clSetKernelArg(__gpu_score, 9, sizeof(cl_mem), params.numCtermNL ? &cl_nlValuesCterm : NULL);
    err |= clSetKernelArg(__gpu_score, 10, sizeof(cl_mem), params.numAANL ? &cl_nlValuesAA : NULL);
    check_cl_error(__FILE__, __LINE__, err, "Unable to set constant args (__gpu_score)");
    
}

void Device::create_kernels() {
    int err;
    __gpu_build = clCreateKernel(clProgram, "gpu_build", &err);
    check_cl_error(__FILE__, __LINE__, err, "Unable to create kernel");
    __gpu_transform = clCreateKernel(clProgram, "gpu_transform", &err);
    check_cl_error(__FILE__, __LINE__, err, "Unable to create kernel");
    __gpu_transform_local = clCreateKernel(clProgram, "gpu_transform_local", &err);
    check_cl_error(__FILE__, __LINE__, err, "Unable to create kernel");
    __gpu_score = clCreateKernel(clProgram, "gpu_score", &err);
    check_cl_error(__FILE__, __LINE__, err, "Unable to create kernel");
    __gpu_score_reduction = clCreateKernel(clProgram, "gpu_score_reduction", &err);
    check_cl_error(__FILE__, __LINE__, err, "Unable to create kernel");
    __cl_memset = clCreateKernel(clProgram, "cl_memset", &err);
    check_cl_error(__FILE__, __LINE__, err, "Unable to create kernel");
}

Device::Device(int platformID, int deviceInd, const char* filename, unsigned int minScan, unsigned int maxScan) {
    cl_platform_id platforms[100];
    cl_uint platforms_n = 0;
    int deviceID = config.iDevices[deviceInd];
    int err;

    err = clGetPlatformIDs(100, platforms, &platforms_n);
    check_cl_error(__FILE__, __LINE__, err, "Unable to enumerate OpenCL platforms");

    if (platforms_n == 0) {
        printf("  No OpenCL platforms found!\n");
        printf("  Make sure OpenCL drivers are installed for your hardware.\n");
        tempest_exit(EXIT_FAILURE);
    }
    
    if (platformID >= platforms_n) {
        printf("  Selected platform ID (%d) is out of range!\n", platformID);
        printf("  Available OpenCL platform IDs range from 0 to %d.\n", platforms_n-1);
        tempest_exit(EXIT_FAILURE);
    }

    clPlatformID = platforms[platformID];
    
    cl_device_id devices[100];
    cl_uint devices_n = 0;
    err = clGetDeviceIDs(clPlatformID, CL_DEVICE_TYPE_ALL, 100, devices, &devices_n);
    check_cl_error(__FILE__, __LINE__, err, "Unable to enumerate OpenCL devices");

    if (deviceID >= devices_n) {
        printf("  Selected device ID (%d) is out of range!\n", deviceID);
        printf("  Available OpenCL device IDs on platform %d range from 0 to %d.\n", platformID, devices_n-1);
        tempest_exit(EXIT_FAILURE);
    }

    clDeviceID = devices[deviceID];
    
    char buffer[10240];
    clGetDeviceInfo(clDeviceID, CL_DEVICE_NAME, sizeof(buffer), buffer, NULL);
    printf(" » (%d:%d) Using %s\n", platformID, deviceID, buffer);

    //create context
    clContext = clCreateContext(NULL, 1, devices+deviceID, NULL, NULL, &err);
    check_cl_error(__FILE__, __LINE__, err, "Error creating context");

    /* Create a command queue */
    clCommandQueue = clCreateCommandQueue(clContext, clDeviceID, CL_QUEUE_PROFILING_ENABLE, &err);
    if(err < 0) {
        perror("Couldn't create a command queue");
        exit(1);   
    }
    
    err  = clGetDeviceInfo(clDeviceID, CL_DEVICE_MAX_MEM_ALLOC_SIZE,  sizeof(cl_ulong), &(lMaxMemAllocSize),  NULL);
    err |= clGetDeviceInfo(clDeviceID, CL_DEVICE_ADDRESS_BITS,        sizeof(cl_uint),  &(iAddressBits),      NULL);
    err |= clGetDeviceInfo(clDeviceID, CL_DEVICE_HOST_UNIFIED_MEMORY, sizeof(cl_bool),  &(bUnifiedMemory),      NULL); 
    err |= clGetDeviceInfo(clDeviceID, CL_DEVICE_GLOBAL_MEM_SIZE,     sizeof(cl_ulong), &(lGlobalMemSize),    NULL);  
    err |= clGetDeviceInfo(clDeviceID, CL_DEVICE_MAX_COMPUTE_UNITS,   sizeof(cl_uint),  &(iMaxComputeUnits),  NULL);
    err |= clGetDeviceInfo(clDeviceID, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(cl_ulong), &(lMaxWorkGroupSize), NULL);
    check_cl_error(__FILE__, __LINE__, err, "Unable to get OpenCL device attributes");

    //limit global memory size to addressable memory (e.g. NVidia may list more than is addressable)
    if (iAddressBits == 32)
        if (lGlobalMemSize > pow(2, 32))
            lGlobalMemSize = pow(2, 32);
    
    // Determine Configuration
    // peaks
    size_t size_iPeakBins = tempest.lNumMS2Peaks * sizeof(cl_int);
    size_t size_fPeakInts = tempest.lNumMS2Peaks * sizeof(cl_float);
    cl_iPeakBins = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, size_iPeakBins, &(host_iPeakBins[0]), &err);
    check_cl_error(__FILE__, __LINE__, err, "Unable to allocate device memory for peak bins.");
    cl_fPeakInts = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, size_fPeakInts, &(host_fPeakInts[0]), &err);
    check_cl_error(__FILE__, __LINE__, err, "Unable to allocate device memory for peak intensities.");
    
    // cleanup host
    //std::vector<int>().swap(host_iPeakBins);
    //std::vector<float>().swap(host_fPeakInts);

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
    size_t availMemSpectra = lGlobalMemSize
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
    
    printf(" » (%d:%d) Allocating %ld bytes of device memory for %ld cached spectra.\n", platformID, deviceID, maxCachedSpectra*tempest.iNumMS2Bins*sizeof(cl_float), maxCachedSpectra);
    for (int i=0; i<maxCachedSpectra; i++) {
        cl_mem newBuffer = clCreateBuffer(clContext, CL_MEM_READ_WRITE, tempest.iNumMS2Bins*sizeof(cl_float), NULL, &err);
        check_cl_error(__FILE__, __LINE__, err, "Unable to allocate spectrum memory on device.");
        unusedBuffers.push(newBuffer);
    }

    FILE *program_handle;
    char *program_buffer, *program_log;
    size_t program_size, log_size;
    char options[1024];

    /* Read program file and place content into buffer */
    program_handle = fopen(filename, "r");
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
    clProgram = clCreateProgramWithSource(clContext, 1, 
                                        (const char**)&program_buffer, &program_size, &err);
    if(err < 0) {
        perror("Couldn't create the .cl program");
        exit(1);
    }
    free(program_buffer);

    /* Define known constants as literals to build into the program */
    float fMassProteinNterm = (float) params.dProteinNtermMass;
    float fMassPeptideNterm = (float) params.dPeptideNtermMass;
    float fMassNtermMod = (float) params.dVariableNtermMassDiff;
    sprintf(options,                   "-DMAX_PEPTIDE_LENGTH=%d ",     MAX_PEPTIDE_LENGTH);
    sprintf(options + strlen(options), "-DBLOCKDIM_REDUCE=%lu ",        BLOCKDIM_REDUCE);
    sprintf(options + strlen(options), "-DBLOCKDIM_SCORE=%lu " ,        DEFAULT_BLOCKDIM_SCORE);
    sprintf(options + strlen(options), "-DGPU_MASS_PROTEIN_NTERM=%f ", fMassProteinNterm);
    sprintf(options + strlen(options), "-DGPU_MASS_PEPTIDE_NTERM=%f ", fMassPeptideNterm);
    sprintf(options + strlen(options), "-DGPU_MASS_NTERM_MOD=%f ",     fMassNtermMod);
    sprintf(options + strlen(options), "-DGPU_TOLERANCE=%f ",          params.fFragmentTolerance);
    sprintf(options + strlen(options), "-DGPU_NUM_OUTPUT_PSMS=%d ",    params.iNumOutputPSMs);
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
    setup_constant_memory();

    for (int i=minScan+deviceInd; i<maxScan; i+=config.iDevices.size()) {
        eObj* e = eScans[i];
        e->clEventSent = clCreateUserEvent(clContext, NULL);
        clSetUserEventStatus(e->clEventSent, 0);
        e->device = this;
    }
    
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
    static size_t globalTransDim = mround(tempest.iNumMS2Bins, BLOCKDIM_TRANSFORM);
    static size_t zMemShared;
    static float fElapsedTime;
    long lSpectrumOffset = e->lIndex*tempest.iNumMS2Bins;
    long lScratchOffset = (long)tempest.iCrossCorrelationWidth;
    long lNoOffset = 0;
    int err;
    cl_ulong start;
    cl_ulong end;
    
    err = clEnqueueWriteBuffer(clCommandQueue, cl_cCandidates, CL_FALSE, 0, sizeof(cObj) * e->iNumBufferedCandidates, e->pCandidateBufferFill, 0, NULL, &(e->clEventSent));
    check_cl_error(__FILE__, __LINE__, err, "Unable to copy candidate data from host to GPU");
	
    stGlobalDim = mround(host_iPeakCounts[e->lIndex], BLOCKDIM_BUILD);
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
        err = clEnqueueCopyBuffer(clCommandQueue, cl_init_fSpectra, spectrumBuffer, 0, 0, tempest.iNumMS2Bins*sizeof(cl_float), 0, NULL, PROFILE ? &memsetEvent : NULL);
        //check_cl_error(__FILE__, __LINE__, err, "Unable to clear spectrum memory");
        if (err != 0) {
            //memory cap reached. Stop filling new buffers.
            unusedBuffers = std::stack<cl_mem>();
            spectrumBuffer = spectrum2buffer.begin()->second;
            spectrum2buffer.erase(spectrum2buffer.begin());
            spectrum2buffer[e->lIndex] = spectrumBuffer;
            err = clEnqueueCopyBuffer(clCommandQueue, cl_init_fSpectra, spectrumBuffer, 0, 0, tempest.iNumMS2Bins*sizeof(cl_float), 0, NULL, PROFILE ? &memsetEvent : NULL);
            check_cl_error(__FILE__, __LINE__, err, "Unable to clear spectrum memory");
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
        err |= clSetKernelArg(__gpu_build, 1, sizeof(long), &lNoOffset);
        err |= clSetKernelArg(__gpu_build, 2, sizeof(int), &(host_iPeakCounts[e->lIndex]));
        err |= clSetKernelArg(__gpu_build, 5, sizeof(long), &(host_lPeakIndices[e->lIndex]));
        err |= clEnqueueNDRangeKernel(clCommandQueue, __gpu_build, 1, NULL, &stGlobalDim, &BLOCKDIM_BUILD, 0, NULL, PROFILE ? &buildEvent : NULL);
        check_cl_error(__FILE__, __LINE__, err, "Could not build spectrum (gpu_build kernel)");
        if (PROFILE) {
            clFinish(clCommandQueue);
            clGetEventProfilingInfo(buildEvent, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, NULL);
            clGetEventProfilingInfo(buildEvent, CL_PROFILING_COMMAND_END,   sizeof(cl_ulong), &end,   NULL);
            totalBuildTime += (end-start);
            buildLaunches += 1;
            clReleaseEvent(buildEvent);
        }

        // transform
        if (params.bCrossCorrelation) {
            size_t localDim = CROSS_CORRELATION_WINDOW * 2;
            size_t globalDim = localDim * tempest.iNumMS2Bins;
            err  = clSetKernelArg(__gpu_transform, 0, sizeof(cl_mem), &spectrumBuffer);
            err |= clSetKernelArg(__gpu_transform, 1, sizeof(long), &lNoOffset);
            err |= clEnqueueNDRangeKernel(clCommandQueue, __gpu_transform, 1, NULL, &globalDim, &localDim, 0, NULL, PROFILE ? & transformEvent : NULL);
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
    err |= clEnqueueNDRangeKernel(clCommandQueue, __gpu_score, 1, NULL, &DEFAULT_CANDIDATE_BUFFER_SIZE, &DEFAULT_BLOCKDIM_SCORE, 0, NULL, PROFILE ? &scoreEvent : NULL);
    check_cl_error(__FILE__, __LINE__, err, "Could not score candidates (gpu_score kernel)");
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
    zMemShared = BLOCKDIM_REDUCE * (sizeof(int) + sizeof(float)) + sizeof(int);
    long lPSMsOffset = e->lIndex * params.iNumOutputPSMs;
    err |= clSetKernelArg(__gpu_score_reduction, 0, sizeof(int), &(config.iCandidateBufferSize));
    err |= clSetKernelArg(__gpu_score_reduction, 1, sizeof(cl_mem), &cl_cCandidates);
    err |= clSetKernelArg(__gpu_score_reduction, 2, sizeof(cl_mem), &cl_fScores);
    err |= clSetKernelArg(__gpu_score_reduction, 3, sizeof(cl_mem), &cl_mPSMs);
    err |= clSetKernelArg(__gpu_score_reduction, 4, sizeof(long), &lPSMsOffset);
    err |= clSetKernelArg(__gpu_score_reduction, 5, sizeof(cl_mem), &cl_fNextScores);
    err |= clSetKernelArg(__gpu_score_reduction, 6, sizeof(long), &(e->lIndex));
    err |= clSetKernelArg(__gpu_score_reduction, 7, zMemShared, NULL);
    err |= clEnqueueNDRangeKernel(clCommandQueue, __gpu_score_reduction, 1, NULL, &BLOCKDIM_REDUCE, &BLOCKDIM_REDUCE, 0, NULL, PROFILE ? &reduceEvent : NULL);
    check_cl_error(__FILE__, __LINE__, err, "Could not process scores (gpu_score_reduction kernel)");
    if (PROFILE) {
        clFinish(clCommandQueue);
        clGetEventProfilingInfo(reduceEvent, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, NULL);
        clGetEventProfilingInfo(reduceEvent, CL_PROFILING_COMMAND_END,   sizeof(cl_ulong), &end,   NULL);
        totalReduceTime += (end-start);
        clReleaseEvent(reduceEvent);
    }
	
    // Stats
    gpu_info.iNumScoringKernels += 1;

    // reset buffer
    e->iNumBufferedCandidates = 0;
}

void Device::finish() {
    clFinish(clCommandQueue);
}

int Device::get_mPSMs(mObj* destination) {
    return clEnqueueReadBuffer(clCommandQueue, cl_mPSMs, CL_TRUE, 0, sizeof(mObj)*tempest.iNumSpectra*params.iNumOutputPSMs, destination, 0, NULL, NULL);
}

int Device::get_fNextScores(float* destination) {
    return clEnqueueReadBuffer(clCommandQueue, cl_fNextScores, CL_TRUE, 0, sizeof(float)*tempest.iNumSpectra, destination, 0, NULL, NULL);
}

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
 err |= clReleaseMemObject(cl_fScratch);
 err |= clReleaseMemObject(cl_cCandidates);
 err |= clReleaseMemObject(cl_fScores);
 err |= clReleaseMemObject(cl_mPSMs);
 err |= clReleaseMemObject(cl_fNextScores);
 err |= clReleaseMemObject(cl_bTheoData);
 check_cl_error(__FILE__, __LINE__, err, "Unable to free memory object(s)");
 }
*/
