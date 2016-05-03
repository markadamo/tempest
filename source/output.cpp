#include "tempest.hpp"

#ifdef MAC
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#define header "Scan\tPrecursorMass\tCharge\tCandidates\tRank\tSequence\tRef\tDecoy\tOccurrences\tSequenceMass\tppm\tScore\tDeltaScore\tDeltaScore_backbone\n"

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
    if (0 == (mPSMs = (mObj*) malloc(sizeof(mObj)*Tempest::data.iNumSpectra*Tempest::params.numInternalPSMs))) {
        fprintf(stderr, "\nERROR\tUnable to allocate host memory for results\n");
        Tempest::tempest_exit(EXIT_FAILURE);
    }

    // allocate host memory
    if (0 == (fNextScores = (float*) malloc(sizeof(float)*Tempest::data.iNumSpectra))) {
        fprintf(stderr, "\nERROR\tUnable to allocate host memory for results\n");
        Tempest::tempest_exit(EXIT_FAILURE);
    }

    // transfer results
    err = 0;
    for (int i=0; i<Tempest::config.iDevices.size(); i++)
        err |= Tempest::devices[i]->get_mPSMs(mPSMs);
    //cudaMemcpy(mPSMs, gpu_mPSMs, sizeof(mObj)*Tempest::data.iNumSpectra*Tempest::params.iNumOutputPSMs, cudaMemcpyDeviceToHost);
    //err |= Tempest::devices[0]->get_fNextScores(fNextScores);
    //cudaMemcpy(fNextScores, gpu_fNextScores, sizeof(float)*Tempest::data.iNumSpectra, cudaMemcpyDeviceToHost);
    Tempest::check_cl_error(__FILE__, __LINE__, err, "Unable to transfer results from the device");
    
    // setup single output file
    //if (!Tempest::args.bSplitOutput) {
        //filename
        strcpy(sOutfile, Tempest::args.sOut);
        strcat(sOutfile, ".tsv\0");

        if (Tempest::args.iPrintLevel) {
            printf(" » Writing results to %s...\n", sOutfile);
        }

        // open
        if (0 == (fp = (FILE *) fopen(sOutfile, "w"))) {
            fprintf(stderr, "\nERROR\tUnable to open %s (%s)\n", sOutfile, strerror(errno));
            Tempest::tempest_exit(EXIT_FAILURE);
        }

        // header
        fprintf(fp,header);
        //}

    // loop through MS/MS scans
    for (iBin=0;iBin<Tempest::data.iNumPrecursorBins;iBin++) {
        for (e = Tempest::data.eScanIndex[iBin]; e; e = e->next) {
            // skip if no results
            if (0 == e->iNumCandidates) {
                fprintf(stderr, "WARNING\tno results for %s\n", e->sName);
                continue;
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
                
                fprintf(fp, "%s\t", e->sName);
                fprintf(fp, "%.5f\t", e->dPrecursorMass);
                fprintf(fp, "%d\t", e->iPrecursorCharge);
                fprintf(fp, "%d\t", e->iNumCandidates);
                fprintf(fp, "%d\t", i+1);
                fprintf(fp, "%c.%s.%c\t", mScanPSMs[i].cBefore, standard_peptide(mScanPSMs[i]), mScanPSMs[i].cAfter);
                fprintf(fp, "%s\t", &Tempest::data.sProteinReferences[mScanPSMs[i].iProtein * MAX_LENGTH_REFERENCE]);
                fprintf(fp, "%d\t", mScanPSMs[i].decoy);
                fprintf(fp, "%d\t", mScanPSMs[i].iNumOccurrences);
                fprintf(fp, "%.5f\t", mScanPSMs[i].fPeptideMass);
                fprintf(fp, "%.4f\t", 1000000.0f * (e->dPrecursorMass - mScanPSMs[i].fPeptideMass) / mScanPSMs[i].fPeptideMass);
                fprintf(fp, "%.4f\t", mScanPSMs[i].fScore);
                if (i < Tempest::params.numInternalPSMs-1)
                    fprintf(fp, "%.4f\t", (mScanPSMs[i].fScore - mScanPSMs[i+1].fScore) / mScanPSMs[i].fScore);
                else
                    fprintf(fp, "%.4f\t", 1.0);
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
        }
    }
    
    // cleanup
    //free(mPSMs);
    //free(fNextScores);
}

extern void Tempest::write_log(int argc, char** argv)
{
    int i;
    char sOutfile[STRING_SIZE];
    char sDate[STRING_SIZE];
    char sLine[STRING_SIZE];
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
    //fprintf(log, "License: %s\n", LICENSE_STRING);
    fprintf(log, "\n");

    //date
    tTime = time(0);
    tDate  = localtime(&tTime);
    strftime(sDate, STRING_SIZE, "%D %T", tDate);
    fprintf(log, "Date: %s\n", sDate);
    fprintf(log, "\n");  

    //summary
    fprintf(log, "Spectra:  %d\n", Tempest::data.iNumSpectra);
    fprintf(log, "Masses:   %f - %f Da\n", Tempest::data.fMinPrecursorMass, Tempest::data.fMaxPrecursorMass);
    fprintf(log, "Proteins: %d\n", Tempest::data.iNumProteins);
    fprintf(log, "Peptides: %ld\n", Tempest::data.lNumPeptides);
    fprintf(log, "PSMs:     %ld\n", Tempest::data.lNumPSMs);
    fprintf(log, "Runtime:  %lus\n", (unsigned long int) tTime-Tempest::data.tStart);
        //modifications
    // fixed
    fprintf(log, "Modifications:\n");
    for (i=0; i<Tempest::params.iNumMods; i++) {
        if (!Tempest::params.tMods[i].cSymbol) {
            fprintf(log, "%c\t+ %f\t= %f Da\n", Tempest::params.tMods[i].cAminoAcid, Tempest::params.tMods[i].dMassDiff, Tempest::params.dMassAA[Tempest::params.tMods[i].cAminoAcid]);
        }
    }
    if (Tempest::params.dPeptideNtermMass > 0.0) fprintf(log, "peptide-nterm\t+ %f Da\n", Tempest::params.dPeptideNtermMass);
    if (Tempest::params.dPeptideCtermMass > 0.0) fprintf(log, "peptide-cterm\t+ %f Da\n", Tempest::params.dPeptideCtermMass);
    if (Tempest::params.dProteinNtermMass > 0.0) fprintf(log, "protein-nterm\t+ %f Da\n", Tempest::params.dProteinNtermMass);
    if (Tempest::params.dProteinCtermMass > 0.0) fprintf(log, "protein-cterm\t+ %f Da\n", Tempest::params.dProteinCtermMass);
    // variable
    for (i=0; i<Tempest::params.iNumMods; i++) {
        if (Tempest::params.tMods[i].cSymbol) {
            fprintf(log, "%c%c\t+ %f\t= %f Da\n", Tempest::params.tMods[i].cAminoAcid, Tempest::params.tMods[i].cSymbol, Tempest::params.tMods[i].dMassDiff, Tempest::params.dMassAA[Tempest::params.tMods[i].cAminoAcid]+Tempest::params.tMods[i].dMassDiff);
        }
    }
    for (i=0; i<Tempest::params.numNtermModSites; i++) {
        fprintf(log, "%s-nterm%c\t+ %f Da\n", Tempest::params.ntermModProtein[i] ? "protein" : "peptide", Tempest::params.ntermModSymbols[i], Tempest::params.ntermModMasses[i]);
    }
    for (i=0; i<Tempest::params.numCtermModSites; i++) {
        fprintf(log, "%s-cterm%c\t+ %f Da\n", Tempest::params.ctermModProtein[i] ? "protein" : "peptide", Tempest::params.ctermModSymbols[i], Tempest::params.ctermModMasses[i]);
    }
    fprintf(log, "\n");
    fprintf(log, "%s\n\n", std::string(80,'=').c_str());

    //command line arguments
    fprintf(log, "%s\n", std::string(80,'-').c_str());
    fprintf(log, "Command line\n");
    fprintf(log, "%s\n", std::string(80,'-').c_str());
    for (int i=0; i<argc; i++)
        fprintf(log, "%s ", argv[i]);
    fprintf(log, "\n");
    fprintf(log, "%s\n\n", std::string(80,'=').c_str());
    
    //method file
    fprintf(log, "%s\n", std::string(80,'-').c_str());
    fprintf(log, "Method file: %s\n", Tempest::args.sParams);
    fprintf(log, "%s\n", std::string(80,'-').c_str());
    FILE* fp = (FILE*) fopen(Tempest::args.sParams, "r");
    while (fgets(sLine, sizeof(sLine), fp))
        fprintf(log, "%s", sLine);
    fprintf(log, "%s\n\n", std::string(80,'=').c_str());

    //config file
    fprintf(log, "%s\n", std::string(80,'-').c_str());
    fprintf(log, "Config file: %s\n", Tempest::args.sConfig);
    fprintf(log, "%s\n", std::string(80,'-').c_str());
    fp = (FILE*) fopen(Tempest::args.sConfig, "r");
    while (fgets(sLine, sizeof(sLine), fp))
        fprintf(log, "%s", sLine);
    fprintf(log, "%s\n\n", std::string(80,'=').c_str());   
}

char* standard_peptide(mObj psm) {
    static int i;
    static unsigned char* c;
    static char *sModPeptide=0;
    
    //if (sModPeptide == 0) sModPeptide = (char*) malloc((MAX_PEPTIDE_LENGTH + Tempest::params.iModificationsMax + 1) * sizeof(char));
    if (sModPeptide == 0) sModPeptide = (char*) malloc((MAX_PEPTIDE_LENGTH*2 + 2 + 6) * sizeof(char));

    i=0;

    //nterm mod symbol to the left of sequence
    sModPeptide[i++] = 'n';
    if (psm.ntermMod) {
        sModPeptide[i++] = Tempest::params.ntermModSymbols[psm.ntermMod-1];
    }
    
    for (c=psm.sPeptide; *c; c++) {
        if (isupper(*c)) {
            sModPeptide[i++] = *c;
        }
        else {
            //AA mod symbol to the right of AA
            sModPeptide[i++] = Tempest::params.unModAA[*c];
            sModPeptide[i++] = Tempest::params.cModSites[Tempest::params.unModAA[*c]][Tempest::getModInd(*c)];
        }
    }

    //cterm mod symbol to the right of sequence
    sModPeptide[i++] = 'c';
    if (psm.ctermMod) {
        sModPeptide[i++] = Tempest::params.ctermModSymbols[psm.ctermMod-1];
    }

    sModPeptide[i] = '\0';

    return sModPeptide;
}
