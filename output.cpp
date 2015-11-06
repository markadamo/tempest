#include "tempest.hpp"

#ifdef MAC
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#define header "scan,precursor mass,charge,candidates,rank,sequence,ref,occurrences,sequence_mass,ppm,xcorr,dcn,dcn_backbone\n"

char* standard_peptide(mObj);

/*
 * Write PSMs to the output file(s)
 */

extern void Tempest::write_psms()
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
    if (0 == (mPSMs = (mObj*) malloc(sizeof(mObj)*Tempest::tempest.iNumSpectra*Tempest::params.numInternalPSMs))) {
        fprintf(stderr, "\nERROR\tUnable to allocate host memory for results\n");
        Tempest::tempest_exit(EXIT_FAILURE);
    }

    // allocate host memory
    if (0 == (fNextScores = (float*) malloc(sizeof(float)*Tempest::tempest.iNumSpectra))) {
        fprintf(stderr, "\nERROR\tUnable to allocate host memory for results\n");
        Tempest::tempest_exit(EXIT_FAILURE);
    }

    // transfer results
    err = 0;
    for (int i=0; i<Tempest::config.iDevices.size(); i++)
        err |= Tempest::devices[i]->get_mPSMs(mPSMs);
    //cudaMemcpy(mPSMs, gpu_mPSMs, sizeof(mObj)*Tempest::tempest.iNumSpectra*Tempest::params.iNumOutputPSMs, cudaMemcpyDeviceToHost);
    //err |= Tempest::devices[0]->get_fNextScores(fNextScores);
    //cudaMemcpy(fNextScores, gpu_fNextScores, sizeof(float)*Tempest::tempest.iNumSpectra, cudaMemcpyDeviceToHost);
    Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to transfer results from the device");
    
    // setup single output file
    if (!Tempest::args.bSplitOutput) {
        //filename
        strcpy(sOutfile, Tempest::args.sOut);
        strcat(sOutfile, ".csv\0");

        if (Tempest::args.iPrintLevel) {
            printf(" » Writing results to %s... ", sOutfile);
        }

        // open
        if (0 == (fp = (FILE *) fopen(sOutfile, "w"))) {
            fprintf(stderr, "\nERROR\tUnable to open %s (%s)\n", sOutfile, strerror(errno));
            Tempest::tempest_exit(EXIT_FAILURE);
        }

        // header
        fprintf(fp,header);
    }

    // loop through MS/MS scans
    for (iBin=0;iBin<Tempest::tempest.iNumPrecursorBins;iBin++) {
        for (e = Tempest::eScanIndex[iBin]; e; e = e->next) {
            // skip if no results
            if (0 == e->iNumCandidates) {
                fprintf(stderr, "WARNING\tno results for %s\n", e->sName);
                continue;
            }

            // open split output file
            if (Tempest::args.bSplitOutput) {
                //filename
                if (!(sscanf(e->sName, "%s.dta", sOutfile) == 1)) {
                    strcpy(sOutfile, Tempest::args.sOut);
                    strcat(sOutfile, ".");
                    strcat(sOutfile, e->sName);
                }
                
                //extension
                strcat(sOutfile,".csv\0");

                //open
                if (0 == (fp = (FILE *) fopen(sOutfile, "w"))) {
                    fprintf(stderr, "\nERROR\tUnable to open %s (%s)\n", Tempest::args.sOut, strerror(errno));
                    Tempest::tempest_exit(EXIT_FAILURE);
                }

                //header
                fprintf(fp, header);

            }
            
            // write matches to file
            mScanPSMs = &(mPSMs[e->lIndex*Tempest::params.numInternalPSMs]);
            iRank=0;
            fPrevScore = FLT_MAX;
            for (i=0; i < Tempest::params.numOutputPSMs; i++) {
                if (mScanPSMs[i].fScore <= 0.0f) {
                    if (i==0)
                        fprintf(stderr, "WARNING\tno results for %s\n", e->sName);
                    break;
                }

                // update rank
                //if(mScanPSMs[i].fScore != fPrevScore) {
                //    iRank += 1;
                //    fPrevScore = mScanPSMs[i].fScore;
                // }
                
                fprintf(fp, "%s,", e->sName);
                fprintf(fp, "%.5f,", e->dPrecursorMass);
                fprintf(fp, "%d,", e->iPrecursorCharge);
                fprintf(fp, "%d,", e->iNumCandidates);
                fprintf(fp, "%d,", i+1);
                fprintf(fp, "%c.%s.%c,", mScanPSMs[i].cBefore, standard_peptide(mScanPSMs[i]), mScanPSMs[i].cAfter);
                fprintf(fp, "%s,", &Tempest::sProteinReferences[mScanPSMs[i].iProtein * MAX_LENGTH_REFERENCE]);
                fprintf(fp, "%d,", mScanPSMs[i].iNumOccurrences);
                fprintf(fp, "%.5f,", mScanPSMs[i].fPeptideMass);
                fprintf(fp, "%.4f,", 1000000.0f * (e->dPrecursorMass - mScanPSMs[i].fPeptideMass) / mScanPSMs[i].fPeptideMass);
                fprintf(fp, "%.4f,", mScanPSMs[i].fScore);
                if (i < Tempest::params.numInternalPSMs-1)
                    fprintf(fp, "%.4f,", (mScanPSMs[i].fScore - mScanPSMs[i+1].fScore) / mScanPSMs[i].fScore);
                else
                    fprintf(fp, "%.4f,", 1.0);
                if (i < Tempest::params.numInternalPSMs-1) {
                    float nextScore_backbone = 0;
                    for (int j=i+1; j<Tempest::params.numInternalPSMs; j++) {
                        if (!Tempest::backboneMatch(mScanPSMs[i], mScanPSMs[j])) {
                            nextScore_backbone = mScanPSMs[j].fScore;
                            break;
                        }
                    }
                    fprintf(fp, "%.4f", (mScanPSMs[i].fScore - nextScore_backbone) / mScanPSMs[i].fScore);
                }
                else
                    fprintf(fp, "%.4f", 1.0);
                fprintf(fp, "\n");
            }

            // close split output file
            if (Tempest::args.bSplitOutput) fclose(fp);
        }
    }
    
    // close single output file
    if (!Tempest::args.bSplitOutput) {
        if (Tempest::args.iPrintLevel) printf("done\n");
        if (fp) fclose(fp);
    }
    
    // cleanup
    //free(mPSMs);
    //free(fNextScores);
}

extern void Tempest::write_log()
{
    int i;
    char sOutfile[STRING_SIZE];
    char sDate[STRING_SIZE];
    FILE *log;
    time_t tTime=0;
    struct tm *tDate;

    //filename
    strcpy(sOutfile, Tempest::args.sOut);
    strcat(sOutfile, ".log\0");
    
    // open log file
    if (0 == (log = (FILE *) fopen(sOutfile, "w"))) {
        fprintf(stderr, "\nERROR\tUnable to open %s for writing (%s)\n", sOutfile, strerror(errno));
        Tempest::tempest_exit(EXIT_FAILURE);
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
    fprintf(log, "Spectra:             %s\n", Tempest::args.sSpectra);
    fprintf(log, "Fasta:               %s\n", Tempest::args.sFasta);
    fprintf(log, "Digest Sites:        %s\n", Tempest::params.sDigestSites);
    fprintf(log, "Non-Digest Sites:    %s\n", Tempest::params.sDigestNoSites);
    fprintf(log, "Digest Offset:       %d\n", Tempest::params.iDigestOffset);
    switch (Tempest::params.iDigestSpecificity) {
    case 0: fprintf(log, "Digest Specificity:  None\n"); break;
    case 1: fprintf(log, "Digest Specificity:  Partial\n"); break;
    case 2: fprintf(log, "Digest Specificity:  Full\n"); break;
    }
    fprintf(log, "Missed Cleavages:    %d\n", Tempest::params.iDigestMaxMissedCleavages);
    fprintf(log, "Digest Length:       %d-%d\n", Tempest::params.iMinPeptideLength, Tempest::params.iMaxPeptideLength);

    for (i=0; i<Tempest::params.iNumMods; i++) {
        if (!Tempest::params.tMods[i].cSymbol) {
            fprintf(log, "Fixed Modification:  %c %.4f Da\n", Tempest::params.tMods[i].cAminoAcid, Tempest::params.tMods[i].dMassDiff);
        }
    }

    if (Tempest::params.dPeptideNtermMass > 0.0) fprintf(log, "Fixed Modification:  peptide-nterm %f Da\n", Tempest::params.dPeptideNtermMass);
    if (Tempest::params.dPeptideCtermMass > 0.0) fprintf(log, "Fixed Modification:  peptide-cterm %f Da\n", Tempest::params.dPeptideCtermMass);
    if (Tempest::params.dProteinNtermMass > 0.0) fprintf(log, "Fixed Modification:  protein-nterm %f Da\n", Tempest::params.dProteinNtermMass);
    if (Tempest::params.dProteinCtermMass > 0.0) fprintf(log, "Fixed Modification:  protein-cterm %f Da\n", Tempest::params.dProteinCtermMass);

    for (i=0; i<Tempest::params.iNumMods; i++) {
        if (Tempest::params.tMods[i].cSymbol) {
            fprintf(log, "Modification:        %c%c %.4f Da\n", Tempest::params.tMods[i].cAminoAcid, Tempest::params.tMods[i].cSymbol, Tempest::params.tMods[i].dMassDiff);
        }
    }

    for (i=0; i<Tempest::numNtermModSites; i++) {
        fprintf(log, "Modification:        nterm %c %.4f Da\n", Tempest::ntermModSymbols[i], Tempest::ntermModMasses[i]);
    }

    for (i=0; i<Tempest::numCtermModSites; i++) {
        fprintf(log, "Modification:        cterm %c %.4f Da\n", Tempest::ctermModSymbols[i], Tempest::ctermModMasses[i]);
    }

    fprintf(log, "Precursor Tolerance: %.4f %s\n", Tempest::params.fPrecursorTolerance, Tempest::params.bPrecursorTolerancePPM ? "ppm" : "Da");
    fprintf(log, "Fragment Tolerance:  %.4f %s\n", Tempest::params.fFragmentTolerance, Tempest::params.bFragmentTolerancePPM ? "ppm" : "mz");
    fprintf(log, "Remove Precursors:   %d\n", Tempest::params.bRemovePrecursors);

    for (i=0; i<Tempest::params.iNumRemoveMzRanges; i++) {
        fprintf(log, "Remove Peaks Range:  %.4f-%.4f m/z\n", Tempest::params.fRemoveMzRanges[2*i], Tempest::params.fRemoveMzRanges[2*i+1]);
    }

    fprintf(log, "xcorr_transform_window: %d\n", Tempest::params.xcorrTransformWidth);
    fprintf(log, "flanking_intensity:     %f\n", Tempest::params.flankingIntensity);
    fprintf(log, "num_internal_psms:      %d\n", Tempest::params.numInternalPSMs);
    fprintf(log, "num_output_psms:        %d\n", Tempest::params.numOutputPSMs);
    fprintf(log, "backbone_delta_score:   %d\n", Tempest::params.bFixDeltaScore);
    fprintf(log, "\n");

    
    //fprintf(log, " Buffer:    %lu\n", Tempest::config.iCandidateBufferSize);
    //fprintf(log, " Config:    %lu blocks x %lu threads per block\n", Tempest::config.iScoreNumBlocks, Tempest::config.iScoreBlockDim);
        
    fprintf(log, "\n");

    //summary
    fprintf(log, "Spectra:  %d\n", Tempest::tempest.iNumSpectra);
    fprintf(log, "Masses:   %.5f - %.5f Da\n", Tempest::tempest.fMinPrecursorMass, Tempest::tempest.fMaxPrecursorMass);
    fprintf(log, "Proteins: %d\n", Tempest::tempest.iNumProteins);
    fprintf(log, "Peptides: %ld\n", Tempest::tempest.lNumPeptides);
    fprintf(log, "PSMs:     %ld\n", Tempest::tempest.lNumPSMs);
    fprintf(log, "Runtime:  %lus\n", (unsigned long int) tTime-Tempest::tempest.tStart);
    
}

char* standard_peptide(mObj psm) {
    static int i;
    static unsigned char* c;
    static char *sModPeptide=0;
    
    if (sModPeptide == 0) sModPeptide = (char*) malloc((MAX_PEPTIDE_LENGTH + Tempest::params.iModificationsMax + 1) * sizeof(char)); 
    for (c=psm.sPeptide, i=0; *c; c++) {
        if (isupper(*c)) {
            sModPeptide[i++] = *c;
        }
        else {
            sModPeptide[i++] = Tempest::unModAA[*c];
            sModPeptide[i++] = Tempest::cModSites[Tempest::unModAA[*c]][Tempest::getModInd(*c)];
        }

        if (c==psm.sPeptide && psm.ntermMod) {
            sModPeptide[i++] = Tempest::ntermModSymbols[psm.ntermMod-1];
        }
    }

    if (psm.ctermMod) {
        sModPeptide[i++] = Tempest::ctermModSymbols[psm.ctermMod-1];
    }

    sModPeptide[i] = '\0';

    return sModPeptide;
}
