#include "tempest.h"

#ifdef MAC
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

cl_mem GPU_MASS_AA;
cl_mem GPU_MASS_PROTON;
cl_mem GPU_MASS_NH3;
cl_mem GPU_MASS_H2O;
cl_mem GPU_MASS_CO;

// global kernels
cl_kernel __gpu_build;
cl_kernel __gpu_transform;
cl_kernel __gpu_transform_local;
cl_kernel __gpu_score;
cl_kernel __gpu_score_reduction;
cl_kernel __cl_memset;

std::map<long, cl_mem> spectrum2buffer;

cl_event scoreEvent = clCreateUserEvent(clContext, NULL);
cl_event reduceEvent = clCreateUserEvent(clContext, NULL);
cl_event buildEvent = clCreateUserEvent(clContext, NULL);
cl_event memsetEvent = clCreateUserEvent(clContext, NULL);
cl_event transformEvent = clCreateUserEvent(clContext, NULL);
long totalScoreTime = 0;
long totalReduceTime = 0;
long totalBuildTime = 0;
long totalTransformTime = 0;
long totalMemsetTime = 0;
long totalSendTime = 0;
long buildLaunches = 0;
long scoreKernelLaunches = 0;
long lastBuildIndex = -1;

//=================================================================================================
// Entry Functions
//=================================================================================================

/*
 * Launches a kernel to perform parallel memset on a block of device memory
 */

int cl_memset(cl_mem buffer, int c, unsigned long n, cl_event* evt) {
    int err;
    size_t global_dim = mround(n, MAX_BLOCKDIM);
    err  = clSetKernelArg(__cl_memset, 0, sizeof(cl_mem), &buffer);
    err |= clSetKernelArg(__cl_memset, 1, sizeof(int), &c);
    err |= clSetKernelArg(__cl_memset, 2, sizeof(long), &n);
    err |= clEnqueueNDRangeKernel(clCommandQueue, __cl_memset, 1, NULL, &global_dim, &MAX_BLOCKDIM, 0, NULL, evt);
    return err;
}

/*
 * Copies masses and other constant values into device constant memory;
 */

void setup_constant_memory() {
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
    err |= clSetKernelArg(__gpu_score, 8, sizeof(cl_mem), &cl_nlValuesNterm);
    err |= clSetKernelArg(__gpu_score, 9, sizeof(cl_mem), &cl_nlValuesCterm);
    err |= clSetKernelArg(__gpu_score, 10, sizeof(cl_mem), &cl_nlValuesAA);
    check_cl_error(__FILE__, __LINE__, err, "Unable to set constant args (__gpu_score)");
    
}

void create_kernels() {
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

/*
 * Launches transfers and kernels for gpu scoring
 */

void gpu_score_candidates(eObj *e)
{
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
