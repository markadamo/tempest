#include "tempest.hpp"

#define LENGTH_INPUT_BUFFER 1024

//functions
int count_references(FILE*);
void parse_reference(char*, char*);
void digest_protein(int, int, char*);
void gen_candidates(cObj, bool, bool, bool, bool, int*, int, int);
void store_candidate(cObj cCandidate);
void write_to_buffer(eObj* e, cObj cCandidate);

/*
 *  Parse fasta database.
 */

extern void Tempest::search_fasta_database()
{
    int i;
    int iProtein;
    int iLengthProtein;
    int iLengthBuffer;
    int iProteinBuffer = INITIAL_PROTEIN_BUFFER;
    char sLineBuffer[LENGTH_INPUT_BUFFER+1];
    char sProteinReference[MAX_LENGTH_REFERENCE+1];
    char *sProteinSequence;
    FILE *fp;
    eObj *e;
    float fProgress;

    // Setup
    // allocate working memory
    if (0 == (sProteinSequence = (char*) malloc(iProteinBuffer*sizeof(char)))) {
        fprintf(stderr, "\nERROR\tUnable to allocate memory for the active protein sequence.\n");
        Tempest::tempest_exit(EXIT_FAILURE);
    }

    // calculate peptide mass range
    float fMaxAAModMass = 0.0f;
    for (i=0; i<Tempest::params.iNumMods; i++) {
        if (Tempest::params.tMods[i].cSymbol && Tempest::params.tMods[i].dMassDiff > fMaxAAModMass)
            fMaxAAModMass = Tempest::params.tMods[i].dMassDiff;
    }
    float fMaxNtermModMass = 0.0f;
    for (i=0; i<Tempest::params.numNtermModSites; i++)
        if (Tempest::params.ntermModMasses[i] > fMaxNtermModMass)
            fMaxNtermModMass = Tempest::params.ntermModMasses[i];
    float fMaxCtermModMass = 0.0f;
    for (i=0; i<Tempest::params.numCtermModSites; i++)
        if (Tempest::params.ctermModMasses[i] > fMaxCtermModMass)
            fMaxCtermModMass = Tempest::params.ctermModMasses[i];

    float maxMaxTermModMass = fMaxNtermModMass > fMaxCtermModMass ? fMaxNtermModMass : fMaxCtermModMass;
    float minMaxTermModMass = fMaxNtermModMass < fMaxCtermModMass ? fMaxNtermModMass : fMaxCtermModMass;

    float fMaxMassDiff = 0.0;
    int remainingMods = Tempest::params.iModificationsMax;
    if (remainingMods && maxMaxTermModMass > fMaxAAModMass) {
        fMaxMassDiff += maxMaxTermModMass;
        remainingMods -= 1;
    }
    if (remainingMods && minMaxTermModMass > fMaxAAModMass) {
        fMaxMassDiff += minMaxTermModMass;
        remainingMods -= 1;
    }
    fMaxMassDiff += fMaxAAModMass * remainingMods;

    if (Tempest::params.bPrecursorTolerancePPM) {
        Tempest::data.fMinPrecursorMass -= Tempest::data.fMinPrecursorMass*Tempest::params.fPrecursorTolerance/1000000.0f;
        Tempest::data.fMaxPrecursorMass += Tempest::data.fMaxPrecursorMass*Tempest::params.fPrecursorTolerance/1000000.0f;
    }
    else {
        Tempest::data.fMinPrecursorMass -= Tempest::params.fPrecursorTolerance;
        Tempest::data.fMaxPrecursorMass += Tempest::params.fPrecursorTolerance;
    }
    Tempest::data.dMinPeptideMass = Tempest::data.fMinPrecursorMass - fMaxMassDiff;
    Tempest::data.dMaxPeptideMass = Tempest::data.fMaxPrecursorMass;
    
    // open fasta file
    // command line overrides param file
    char* fastaFile = Tempest::args.sFasta == NULL ? Tempest::params.sFasta : Tempest::args.sFasta; 
    if (0 == (fp = (FILE *) fopen(fastaFile, "r"))) {
        fprintf(stderr, "\nERROR\tUnable to open fasta database: %s\n", fastaFile);
        Tempest::tempest_exit(EXIT_FAILURE);
    }

    // count proteins
    Tempest::data.iNumProteins = count_references(fp);
    
    if (Tempest::data.iNumProteins == 0) {
        fprintf(stderr, "\nERROR\tNo proteins found in fasta database: %s\n", Tempest::args.sFasta);
        Tempest::tempest_exit(EXIT_FAILURE);
    }

    // progress
    if (Tempest::args.iPrintLevel && !Tempest::args.bPrintCandidates) {
        printf(" » Digesting %d proteins...     ", Tempest::data.iNumProteins);
    }
    
    // Get memory for the protein references
    if (0 == (Tempest::data.sProteinReferences = (char *) malloc(Tempest::data.iNumProteins * MAX_LENGTH_REFERENCE * sizeof(char)))) {
        fprintf(stderr, "\nERROR\tUnable to allocate memory for protein references.\n");
        Tempest::tempest_exit(EXIT_FAILURE);
    }
    
    // parse fasta line by line
    iProtein = 0;
    iLengthProtein = 0;
    sProteinSequence[0] = '\0';
    while (fgets(sLineBuffer, LENGTH_INPUT_BUFFER, fp)) {
        if (strlen(sLineBuffer) <= 0) continue;
        
        // reference line
        if ('>' == sLineBuffer[0]) {
            // check for previous protein
            if (iLengthProtein > 0) {
                // process
                strcpy(&Tempest::data.sProteinReferences[iProtein * MAX_LENGTH_REFERENCE * sizeof(char)], sProteinReference);
                digest_protein(iProtein++, iLengthProtein, sProteinSequence);

                if (Tempest::args.iPrintLevel && !Tempest::args.bPrintCandidates) {
                    fProgress = (100.0f * iProtein) / Tempest::data.iNumProteins;
                    if (fProgress - floor(fProgress) < 0.01f) printf("\b\b\b\b%*d%%", 3, (int) fProgress);
                    fflush(0);
                }
            
                // reset
                iLengthProtein = 0;
                sProteinSequence[0] = '\0';
            }

            // get next reference
            parse_reference(sLineBuffer, sProteinReference);
        }

        // protein sequence line
        else {
            //get length
            iLengthBuffer = (int) strlen(sLineBuffer);
            
            //strip newlines
            while(iLengthBuffer>0 && !isalnum(sLineBuffer[iLengthBuffer-1])) {
                iLengthBuffer--;
            }

            if (iLengthBuffer <= 0)
                continue;

            // check new sequence length against buffer size
            if (iLengthProtein + iLengthBuffer > iProteinBuffer) {
                //grow protein sequence buffer
                iProteinBuffer *= 2;
                if (0 == (sProteinSequence = (char*) realloc(sProteinSequence, iProteinBuffer*sizeof(char)))) {
                    fprintf(stderr, "\nERROR\tCould not allocate memory for protein %s (length %d)\n", sProteinReference, iProteinBuffer);
                    Tempest::tempest_exit(EXIT_FAILURE);
                }
            }

            // concatenate to sequence
            strncat(sProteinSequence,sLineBuffer,iLengthBuffer);
            iLengthProtein += iLengthBuffer;
        }
    }

    // check for errors
    if (ferror(fp)) {
        fprintf(stderr, "\nERROR\tUnable to read from %s, while reading database.\n", Tempest::args.sFasta);
        Tempest::tempest_exit(EXIT_FAILURE);
    }
    
    // Final protein
    strncpy(&Tempest::data.sProteinReferences[iProtein * MAX_LENGTH_REFERENCE],sProteinReference,MAX_LENGTH_REFERENCE);
    digest_protein(iProtein++, iLengthProtein, sProteinSequence);
    
    if (Tempest::args.iPrintLevel && !Tempest::args.bPrintCandidates) {
        fProgress = (100.0f * iProtein) / Tempest::data.iNumProteins;
        printf("\b\b\b\b%*d%%", 3, (int) fProgress);
    }
        
    // Score Remaining Candidates
    for(i=0;i<Tempest::data.iNumPrecursorBins;i++) {
        for (e = Tempest::data.eScanIndex[i]; e; e = e->next) {
            //printf("%d\t%d\n", gpu_info.iNumScoringKernels, i);
            if (e->iNumBufferedCandidates > 0)
                e->device->scoreCandidates(e);
        }
    }
    
    // close fasta
    fclose(fp);

    // summary
    if (Tempest::args.iPrintLevel && !Tempest::args.bPrintCandidates) {
        printf("\n");
    }

    //cleanup
    free(sProteinSequence);
}

/*
 *  Digest protein sequence.
 */

void digest_protein(int iProtein, int iLengthProtein, char *sProteinSequence)
{
    int iStart, iEnd, iPeptideLength, iEndLength, iNextNterm, iBin, iPattern;
    int iNumAAModTargets, iModMass, iNumMods;
    double dPeptideMass;
    double dBaseMass;
    unsigned int uiMod;
    unsigned int uiMaxModPattern;
    bool bMatched;
    eObj* e;
    bool bNterm;
    int  iDigestSites;
    int modInds[MAX_MOD_TARGETS_PER_PEPTIDE];

    // Initialize some variables.
    dBaseMass = Tempest::params.dPeptideNtermMass + Tempest::params.dPeptideCtermMass + H_MASS + OH_MASS + PROTON_MASS;
    iEndLength = iLengthProtein - 1;
    iStart = 0;
    bNterm = 1;

    // loop through protein
    while (iStart < iLengthProtein) {
        // reset
        iEnd = iStart;
        dPeptideMass = dBaseMass;
        if (iStart==0) dPeptideMass += Tempest::params.dProteinNtermMass;
        iDigestSites = 0;
        iNextNterm = 0;
        iNumAAModTargets = 0;

        for (iEnd=iStart; iEnd < iLengthProtein; iEnd++) {
            unsigned char aa = sProteinSequence[iEnd];
            if (aa < 'A' || aa > 'Z') {
                fprintf(stderr, "\nERROR: Invalid amino acid symbol '%c' in protein sequence (ref=%s)\n", aa, &Tempest::data.sProteinReferences[iProtein * MAX_LENGTH_REFERENCE * sizeof(char)]);
                Tempest::tempest_exit(EXIT_FAILURE);
            }

            iPeptideLength = iEnd - iStart + 1;

            // too long?
            if (iPeptideLength >= Tempest::params.iMaxPeptideLength) {
                break;
            }
            
            // update mass
            dPeptideMass += Tempest::params.dMassAA[aa];
            if (iEnd == iEndLength) dPeptideMass += Tempest::params.dProteinCtermMass;
            
            // too large (no spectra to match it)?
            if (dPeptideMass > Tempest::data.dMaxPeptideMass) {
                break;
            }
            
            // mod target?
            if (Tempest::params.numAAModSites[aa]) {
                // too many mod targets?
                if (iNumAAModTargets >= MAX_MOD_TARGETS_PER_PEPTIDE) {
                    break;
                }
                modInds[iNumAAModTargets] = iEnd - iStart;
                iNumAAModTargets++;
            }
            
            // cleavage site?
            if (Tempest::params.iDigestSpecificity > 0) {
                bool isCleavageSite = false;
                if (Tempest::params.iDigestOffset == 0)
                    isCleavageSite = (iEnd == iEndLength) || (Tempest::params.bDigestSites[sProteinSequence[iEnd+1]] && !Tempest::params.bDigestNoSites[aa]);
                else
                    isCleavageSite = (iEnd == iEndLength) || (Tempest::params.bDigestSites[aa] && !Tempest::params.bDigestNoSites[sProteinSequence[iEnd+1]]);
               
                if (isCleavageSite) {    
                    // set next iStart value
                    if (!iNextNterm) {
                        iNextNterm = iEnd + 1;
                    }

                    // too many missed cleavages?
                    if (iDigestSites > Tempest::params.iDigestMaxMissedCleavages) {
                        break;
                    }

                    iDigestSites++;
                }
                else if (Tempest::params.iDigestSpecificity == 2 || (Tempest::params.iDigestSpecificity == 1 && bNterm == 0)) {
                    continue;
                }
            }
            
            // too short or too small?
            if (iPeptideLength < Tempest::params.iMinPeptideLength) {
                continue;
            }

            //too small?
            if (dPeptideMass < Tempest::data.dMinPeptideMass) {
                continue;
            }
            
            // Potential Peptide found
            Tempest::data.lNumPeptides += 1;

            cObj init_cCandidate;
            strncpy((char*)init_cCandidate.sPeptide, &sProteinSequence[iStart], (size_t)iPeptideLength);
            init_cCandidate.decoy = 0;
            init_cCandidate.sPeptide[iPeptideLength] = '\0';
            init_cCandidate.ntermMod = 0;
            init_cCandidate.ctermMod = 0;
            init_cCandidate.iProtein       = iProtein;
            init_cCandidate.iPeptideLength = iPeptideLength;
            init_cCandidate.fPeptideMass   = float(dPeptideMass);
            init_cCandidate.cBefore        = (iStart == 0)        ? '-' : sProteinSequence[iStart-1];
            init_cCandidate.cAfter         = (iEnd == iEndLength) ? '-' : sProteinSequence[iEnd+1];

            gen_candidates(init_cCandidate, Tempest::params.numNtermModSites>0, iStart==0, Tempest::params.numCtermModSites>0, iEnd==iEndLength, modInds, iNumAAModTargets, 0);
        }

        //find next n-terminal residue

        // Full: next digest site
        if (Tempest::params.iDigestSpecificity == 2) {
            if (iNextNterm > 0) {
                iStart = iNextNterm;
            }
            else if (Tempest::params.iDigestOffset == 0) {
                for (iStart=iEnd+1;
                     iStart < iLengthProtein && (!Tempest::params.bDigestSites[sProteinSequence[iStart]] || Tempest::params.bDigestNoSites[sProteinSequence[iStart-1]]);
                     iStart++);
            }
            else {
                for (iStart=iEnd+1;
                     iStart < iLengthProtein && (!Tempest::params.bDigestSites[sProteinSequence[(iStart-1)]] || Tempest::params.bDigestNoSites[sProteinSequence[iStart]]);
                     iStart++);
            }
        }

        // Partial or Nonspecific: next residue
        else {
            iStart++;
            bNterm = Tempest::params.bDigestSites[sProteinSequence[(iStart-1)]] && !Tempest::params.bDigestNoSites[sProteinSequence[iStart]];
        }
    }
}
    
/*
 * Count reference lines in a fasta file.
 */
int count_references(FILE *fp) {
    int n = 0;
    char c;
    char cPrevious = '\n';
    
    while (EOF != (c = getc (fp))) {
        if (c == '>' && cPrevious == '\n') n++;
        cPrevious = c;
    }
    rewind(fp);
    return n;
}

/*
 * Parse a (truncated) protein reference from a fasta description line.
 *
 *  e.g.
 *  description line:  >P31946|1433B_HUMAN 14-3-3 protein beta/alpha - Homo sapiens
 *  protein reference: P31946|1433B_HUMAN
 *
 *
 *      sLine       null-terminated protein description from fasta file (with a leading '>').
 *      sReference  pointer to pre-allocated memory for the parsed protein reference.
 */

void parse_reference(char* sLine, char* sReference) {
    static int iLengthRef;
    
    //look for first space
    for (iLengthRef=1; iLengthRef<MAX_LENGTH_REFERENCE && sLine[iLengthRef+1]!=' '; iLengthRef++);

    //copy
    strncpy(sReference, sLine+1, iLengthRef);
    sReference[iLengthRef] = '\0';
}

void gen_candidates(cObj cCandidate, bool ntermMod, bool proteinNterm, bool ctermMod, bool proteinCterm, int* modInds, int modIndsLeft, int modCount) {
    if (modCount >= Tempest::params.iModificationsMax)
        store_candidate(cCandidate);
    else if (modIndsLeft > 0) {
        int modInd = *modInds;
        unsigned char aa = cCandidate.sPeptide[modInd];
        for (int i=0; i<Tempest::params.numAAModSites[aa]; i++) {
            cObj new_cCandidate = cCandidate;
            //memcpy(&new_cCandidate, &cCandidate, sizeof(cObj));
            new_cCandidate.sPeptide[modInd] = Tempest::toMod(aa, i+1);
            new_cCandidate.fPeptideMass += Tempest::params.fModValues[aa][i];
            //printf("%c[%d] = %f\n", aa, i, sPeptide[modInd], Tempest::fModValues[aa][i]);
            gen_candidates(new_cCandidate, ntermMod, proteinNterm, ctermMod, proteinCterm, modInds+1, modIndsLeft-1, modCount+1);
        }
        gen_candidates(cCandidate, ntermMod, proteinNterm, ctermMod, proteinCterm, modInds+1, modIndsLeft-1, modCount);
    }
    else if (ntermMod) {
        //cObj new_cCandidate;
        for (int i=0; i<Tempest::params.numNtermModSites; i++) {
            if (!proteinNterm && Tempest::params.ntermModProtein[i])
                continue;
            //memcpy(&new_cCandidate, &cCandidate, sizeof(cObj));
            cObj new_cCandidate = cCandidate;
            new_cCandidate.ntermMod = i+1;
            new_cCandidate.fPeptideMass += Tempest::params.ntermModMasses[i];
            gen_candidates(new_cCandidate, 0, proteinNterm, ctermMod, proteinCterm, modInds, modIndsLeft, modCount+1);
        }
        gen_candidates(cCandidate, 0, proteinNterm, ctermMod, proteinCterm, modInds, modIndsLeft, modCount);
    }
    else if (ctermMod) {
        //cObj new_cCandidate;
        for (int i=0; i<Tempest::params.numCtermModSites; i++) {
            if (!proteinCterm && Tempest::params.ctermModProtein[i])
                continue;
            //memcpy(&new_cCandidate, &cCandidate, sizeof(cObj));
            cObj new_cCandidate = cCandidate;
            new_cCandidate.ctermMod = i+1;
            new_cCandidate.fPeptideMass += Tempest::params.ctermModMasses[i];  
            gen_candidates(new_cCandidate, ntermMod, proteinNterm, 0, proteinCterm, modInds, modIndsLeft, modCount+1);
        }
        gen_candidates(cCandidate, ntermMod, proteinNterm, 0, proteinCterm, modInds, modIndsLeft, modCount);
    }
    else 
        store_candidate(cCandidate);
}

void store_candidate(cObj cCandidate) {
    float fModMass = cCandidate.fPeptideMass;
    if (fModMass < Tempest::data.fMinPrecursorMass || fModMass > Tempest::data.fMaxPrecursorMass)
        return;
    int iModMass = (int) roundf(fModMass / Tempest::data.fPrecursorBinWidth);
    for (int iBin=iModMass-1; iBin<=iModMass+1 && iBin<=Tempest::data.iNumPrecursorBins; iBin++) {
        for (eObj* e = Tempest::data.eScanIndex[iBin]; e != 0; e = e->next) {
            // check mass error
            if (Tempest::params.bPrecursorTolerancePPM) {
                if (1000000.0f*fabs(fModMass - (float) e->dPrecursorMass) / fModMass > Tempest::params.fPrecursorTolerance)
                    continue;
            }
            else {
                if (fabs(fModMass - (float) e->dPrecursorMass) > Tempest::params.fPrecursorTolerance)
                    continue;
            }

            write_to_buffer(e, cCandidate);

            if (Tempest::params.decoySearch) {
                cObj decoyCandidate = cCandidate;
                decoyCandidate.decoy = 1;
                //reverse sequence while keeping enzymatic cleavage sites fixed
                int start = 0;
                int end = cCandidate.iPeptideLength-1;
                if (Tempest::params.iDigestSpecificity > 0)
                    if (Tempest::params.iDigestOffset == 0 && cCandidate.cBefore != '-' && Tempest::params.bDigestSites[cCandidate.sPeptide[start]])
                        start++;
                    else if (Tempest::params.iDigestOffset == 1 && cCandidate.cAfter != '-' && Tempest::params.bDigestSites[cCandidate.sPeptide[end]])
                        end--;
                for (; end>=0; start++, end--)
                    decoyCandidate.sPeptide[start] = cCandidate.sPeptide[end];

                write_to_buffer(e, decoyCandidate);
            }
        }
    }
}

void write_to_buffer(eObj* e, cObj cCandidate) {
    Tempest::data.lNumPSMs += 1;
    if (e->iNumBufferedCandidates == 0) {
        clWaitForEvents(1, &(e->clEventSent));
        if (Tempest::config.profile) {
          cl_ulong start;
          cl_ulong end;
          int err;
          err = clGetEventProfilingInfo(e->clEventSent, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, NULL);
          err |= clGetEventProfilingInfo(e->clEventSent, CL_PROFILING_COMMAND_END,   sizeof(cl_ulong), &end,   NULL);
          if (err == 0)
          e->device->totalSendTime += (end-start);
          }
        clReleaseEvent(e->clEventSent);
    }
    e->candidateBuffer[e->iNumBufferedCandidates] = cCandidate;
    //memcpy(e->candidateBuffer+e->iNumBufferedCandidates, &cCandidate, sizeof(cObj));
    e->iNumCandidates++;
    e->iNumBufferedCandidates++;
    if (e->iNumBufferedCandidates == e->candidateBufferSize) {
        //printf("%d\t%d\n", gpu_info.iNumScoringKernels, iBin);
        e->device->scoreCandidates(e);
    }
}
