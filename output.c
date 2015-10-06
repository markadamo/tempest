
#include "tempest.h"

#ifdef MAC
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#define header "scan,precursor mass,charge,candidates,rank,sequence,ref,occurrences,sequence_mass,ppm,xcorr,dcn\n"

char* standard_peptide(mObj);

/*
 * Write PSMs to the output file(s)
 */

extern void write_psms()
{
    char    sOutfile[STRING_SIZE];
    eObj    *e;
    mObj    *mPSMs;
    mObj    *mScanPSMs;
    float   *fNextScores;
    FILE    *fp;
    int     i;
    int     iRank;
    float   fPrevScore;
    int     iBin;
    int     err;
    
    fp = 0;
    
    // allocate host memory
    if (0 == (mPSMs = (mObj*) malloc(sizeof(mObj)*tempest.iNumSpectra*params.iNumOutputPSMs))) {
        fprintf(stderr, "\nERROR\tUnable to allocate host memory for results\n");
        tempest_exit(EXIT_FAILURE);
    }

    // allocate host memory
    if (0 == (fNextScores = (float*) malloc(sizeof(float)*tempest.iNumSpectra))) {
        fprintf(stderr, "\nERROR\tUnable to allocate host memory for results\n");
        tempest_exit(EXIT_FAILURE);
    }

    // transfer results
    err = clEnqueueReadBuffer(clCommandQueue, cl_mPSMs, CL_TRUE, 0, sizeof(mObj)*tempest.iNumSpectra*params.iNumOutputPSMs, mPSMs, 0, NULL, NULL);
    //cudaMemcpy(mPSMs, gpu_mPSMs, sizeof(mObj)*tempest.iNumSpectra*params.iNumOutputPSMs, cudaMemcpyDeviceToHost);
    err |= clEnqueueReadBuffer(clCommandQueue, cl_fNextScores, CL_TRUE, 0, sizeof(float)*tempest.iNumSpectra, fNextScores, 0, NULL, NULL);
    //cudaMemcpy(fNextScores, gpu_fNextScores, sizeof(float)*tempest.iNumSpectra, cudaMemcpyDeviceToHost);
    check_cl_error(__FILE__, __LINE__, err, "Unable to transfer results from the device");
    
    // setup single output file
    if (!args.bSplitOutput) {
        //filename
        strcpy(sOutfile, args.sOut);
        strcat(sOutfile, ".csv\0");

        if (args.iPrintLevel) {
            printf(" » Writing results to %s... ", sOutfile);
        }

        // open
        if (0 == (fp = (FILE *) fopen(sOutfile, "w"))) {
            fprintf(stderr, "\nERROR\tUnable to open %s (%s)\n", sOutfile, strerror(errno));
            tempest_exit(EXIT_FAILURE);
        }

        // header
        fprintf(fp,header);
    }

    // loop through MS/MS scans
    for (iBin=0;iBin<tempest.iNumPrecursorBins;iBin++) {
        for (e = eScanIndex[iBin]; e; e = e->next) {
            // skip if no results
            if (0 == e->iNumCandidates) {
                fprintf(stderr, "WARNING\tno results for %s\n", e->sName);
                continue;
            }

            // open split output file
            if (args.bSplitOutput) {
                //filename
                if (!(sscanf(e->sName, "%s.dta", sOutfile) == 1)) {
                    strcpy(sOutfile, args.sOut);
                    strcat(sOutfile, ".");
                    strcat(sOutfile, e->sName);
                }
                
                //extension
                strcat(sOutfile,".csv\0");

                //open
                if (0 == (fp = (FILE *) fopen(sOutfile, "w"))) {
                    fprintf(stderr, "\nERROR\tUnable to open %s (%s)\n", args.sOut, strerror(errno));
                    tempest_exit(EXIT_FAILURE);
                }

                //header
                fprintf(fp, header);

            }
            
            // write matches to file
            mScanPSMs = &mPSMs[e->lIndex*params.iNumOutputPSMs];
            iRank=1;
            fPrevScore = FLT_MAX;
            for (i=0; i < params.iNumOutputPSMs; i++) {
                if (mScanPSMs[i].fScore == 0.0f) {
                    fprintf(stderr, "WARNING\tno results for %s\n", e->sName);
                    continue;
                }

                fprintf(fp, "%s,", e->sName);
                fprintf(fp, "%.5f,", e->dPrecursorMass);
                fprintf(fp, "%d,", e->iPrecursorCharge);
                fprintf(fp, "%d,", e->iNumCandidates);
                fprintf(fp, "%d,", iRank);
                fprintf(fp, "%c.%s.%c,", mScanPSMs[i].cBefore, standard_peptide(mScanPSMs[i]), mScanPSMs[i].cAfter);
                fprintf(fp, "%s,", &sProteinReferences[mScanPSMs[i].iProtein * MAX_LENGTH_REFERENCE]);
                fprintf(fp, "%d,", mScanPSMs[i].iNumOccurrences);
                fprintf(fp, "%.5f,", mScanPSMs[i].fPeptideMass);
                fprintf(fp, "%.4f,", 1000000.0f * (e->dPrecursorMass - mScanPSMs[i].fPeptideMass) / mScanPSMs[i].fPeptideMass);
                fprintf(fp, "%.4f,", mScanPSMs[i].fScore);
                if (i==0) fprintf(fp, "%.4f", (mScanPSMs[0].fScore - fNextScores[e->lIndex]) / mScanPSMs[0].fScore);
                else fprintf(fp, "n/a");
                fprintf(fp, "\n");

                // update rank
                if(mScanPSMs[i].fScore != fPrevScore) {
                    iRank += 1;
                    fPrevScore = mScanPSMs[i].fScore;
                }
            }

            // close split output file
            if (args.bSplitOutput) fclose(fp);
        }
    }
    
    // close single output file
    if (!args.bSplitOutput) {
        if (args.iPrintLevel) printf("done\n");
        if (fp) fclose(fp);
    }
    
    // cleanup
    //free(mPSMs);
    //free(fNextScores);
}

extern void write_log()
{
    int i;
    char sOutfile[STRING_SIZE];
    char sDate[STRING_SIZE];
    FILE *log;
    time_t tTime=0;
    struct tm *tDate;

    //filename
    strcpy(sOutfile, args.sOut);
    strcat(sOutfile, ".log\0");
    
    // open log file
    if (0 == (log = (FILE *) fopen(sOutfile, "w"))) {
        fprintf(stderr, "\nERROR\tUnable to open %s for writing (%s)\n", sOutfile, strerror(errno));
        tempest_exit(EXIT_FAILURE);
    }

    //Tempest info
    fprintf(log, "Program: Tempest v%s\n", VERSION_STRING);
    fprintf(log, "Authors: %s\n", AUTHORS_STRING);
    fprintf(log, "License: %s\n", LICENSE_STRING);
    fprintf(log, "\n");

    //date
    tTime = time(0);
    tDate  = localtime(&tTime);
    strftime(sDate, STRING_SIZE, "%D %T", tDate);
    fprintf(log, "Date: %s\n", sDate);

    //params
    fprintf(log, "Spectra:             %s\n", args.sSpectra);
    fprintf(log, "Fasta:               %s\n", args.sFasta);
    fprintf(log, "Digest Sites:        %s\n", params.sDigestSites);
    fprintf(log, "Non-Digest Sites:    %s\n", params.sDigestNoSites);
    fprintf(log, "Digest Offset:       %d\n", params.iDigestOffset);
    switch (params.iDigestSpecificity) {
    case 0: fprintf(log, "Digest Specificity:  None\n"); break;
    case 1: fprintf(log, "Digest Specificity:  Partial\n"); break;
    case 2: fprintf(log, "Digest Specificity:  Full\n"); break;
    }
    fprintf(log, "Missed Cleavages:    %d\n", params.iDigestMaxMissedCleavages);
    fprintf(log, "Digest Length:       %d-%d\n", params.iMinPeptideLength, params.iMaxPeptideLength);

    for (i=0; i<params.iNumMods; i++) {
        if (!params.tMods[i].cSymbol) {
            fprintf(log, "Fixed Modification:  %c %.4f Da\n", params.tMods[i].cAminoAcid, params.tMods[i].dMassDiff);
        }
    }

    if (params.dPeptideNtermMass > 0.0) fprintf(log, "Fixed Modification:  peptide-nterm %f Da\n", params.dPeptideNtermMass);
    if (params.dPeptideCtermMass > 0.0) fprintf(log, "Fixed Modification:  peptide-cterm %f Da\n", params.dPeptideCtermMass);
    if (params.dProteinNtermMass > 0.0) fprintf(log, "Fixed Modification:  protein-nterm %f Da\n", params.dProteinNtermMass);
    if (params.dProteinCtermMass > 0.0) fprintf(log, "Fixed Modification:  protein-cterm %f Da\n", params.dProteinCtermMass);

    for (i=0; i<params.iNumMods; i++) {
        if (params.tMods[i].cSymbol) {
            fprintf(log, "Modification:        %c%c %.4f Da\n", params.tMods[i].cAminoAcid, params.tMods[i].cSymbol, params.tMods[i].dMassDiff);
        }
    }

    if (params.cVariableNtermSymbol) {
        fprintf(log, "Modification:        nterm %c %.4f Da\n", params.cVariableNtermSymbol, params.dVariableNtermMassDiff);
    }

    if (params.cVariableCtermSymbol) {
        fprintf(log, "Modification:        cterm %c %.4f Da\n", params.cVariableCtermSymbol, params.dVariableCtermMassDiff);
    }

    fprintf(log, "Precursor Tolerance: %.4f %s\n", params.fPrecursorTolerance, params.bPrecursorTolerancePPM ? "ppm" : "Da");
    fprintf(log, "Fragment Tolerance:  %.4f %s\n", params.fFragmentTolerance, params.bFragmentTolerancePPM ? "ppm" : "mz");
    fprintf(log, "Remove Precursors:   %d\n", params.bRemovePrecursors);

    for (i=0; i<params.iNumRemoveMzRanges; i++) {
        fprintf(log, "Remove Peaks Range:  %.4f-%.4f m/z\n", params.fRemoveMzRanges[2*i], params.fRemoveMzRanges[2*i+1]);
    }

    fprintf(log, "Similarity Score:    %s\n", params.bCrossCorrelation ? "xcorr" : "dotproduct");
    fprintf(log, "Duplicate Peak Mode: %s\n", params.bTrackDuplicates ? "max" : "sum");
    fprintf(log, "Flanking Peaks:      %d\n", params.bTheoreticalFlanking);
    //fprintf(log, "Neutral Loss Peaks:  %d\n", params.bTheoreticalNL);
    //fprintf(log, "Phos NL Peaks:       %d\n", params.bTheoreticalPhosNL);
    fprintf(log, "Max Output PSMs:     %d\n", params.iNumOutputPSMs);
    fprintf(log, "Fix Delta Score:     %d\n", params.bFixDeltaScore);
    fprintf(log, "\n");

    //config
    if (config.bForceNoGPU) {
        fprintf(log, "CPU-only mode\n");
    }

    else {
        fprintf(log, "GPU Device: %d\n", config.iDevice);
        fprintf(log, "Setup Mode: %s\n", get_setup_mode(config.eSetup));
        fprintf(log, "Score Mode: %s\n", get_score_mode(config.eScore));
        fprintf(log, " Buffer:    %d\n", config.iCandidateBufferSize);
        fprintf(log, " Config:    %d blocks x %d threads per block\n", config.iScoreNumBlocks, config.iScoreBlockDim);
    }
    fprintf(log, "\n");

    //summary
    fprintf(log, "Spectra:  %d\n", tempest.iNumSpectra);
    fprintf(log, "Masses:   %.5f - %.5f Da\n", tempest.fMinPrecursorMass, tempest.fMaxPrecursorMass);
    fprintf(log, "Proteins: %d\n", tempest.iNumProteins);
    fprintf(log, "Peptides: %ld\n", tempest.lNumPeptides);
    fprintf(log, "PSMs:     %ld\n", tempest.lNumPSMs);
    fprintf(log, "Runtime:  %lus\n", (unsigned long int) tTime-tempest.tStart);
    
}

char* standard_peptide(mObj psm) {
    static int i;
    static char* c;
    static char *sModPeptide=0;
    
    if (sModPeptide == 0) sModPeptide = (char*) malloc((MAX_PEPTIDE_LENGTH + params.iModificationsMax + 1) * sizeof(char)); 
    for (c=psm.sPeptide, i=0; *c; c++) {
        if (isupper(*c)) {
            sModPeptide[i++] = *c;
        }
        else {
            sModPeptide[i++] = toupper(*c);
            sModPeptide[i++] = cModSites[aa2i(toupper(*c))];
        }

        if (c==psm.sPeptide && psm.bNtermMod) {
            sModPeptide[i++] = params.cVariableNtermSymbol;
        }
    }

    if (psm.bCtermMod) {
        sModPeptide[i++] = params.cVariableCtermSymbol;
    }

    sModPeptide[i] = '\0';

    return sModPeptide;
}
