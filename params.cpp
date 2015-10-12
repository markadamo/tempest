#include "tempest.hpp"

#ifdef MAC
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

void append_mod(char, double, char);

extern void parse_params() {
    FILE* fp;
    char sLine[STRING_SIZE];
    char sParam[STRING_SIZE];
    char sTemp[STRING_SIZE];
    int iTemp;

    char sLocation[STRING_SIZE];
    char cSymbol;
    double dMassDiff;

    float *fMzRange;

    char *c;

    if((fp = (FILE*) fopen(args.sParams, "r")) == NULL) {
        fprintf( stderr, "ERROR\tUnable to open params file %s: %s\n", args.sParams, strerror(errno) );
        tempest_exit(EXIT_FAILURE);
    }
    
    // read params line by line
    while (fgets(sLine, sizeof(sLine), fp)) {
        if (sscanf(sLine, "%s", sParam) == -1) continue;
        if (sParam[0] == '#') continue;

        if (strcmp(sParam, "digest_definition") == 0) {
            params.sDigestSites = (char*) malloc(26*sizeof(char));
            params.sDigestNoSites = (char*) malloc(26*sizeof(char));
            if (sscanf(sLine, "digest_definition %s %s %d \n", params.sDigestSites, params.sDigestNoSites, &params.iDigestOffset) != 3) {
                fprintf(stderr, "\nERROR: Could not parse digest definition\n\t>%s", sLine);
                tempest_exit(EXIT_FAILURE);
            }
        }

        else if (strcmp(sParam, "digest_specificity") == 0) {
            if (sscanf(sLine, "digest_specificity %d \n", &params.iDigestSpecificity) != 1) {
                fprintf(stderr, "\nERROR: Could not parse digest specificity\n\t>%s", sLine);
                tempest_exit(EXIT_FAILURE);
            }

            if (params.iDigestSpecificity < 0 || params.iDigestSpecificity > 2) {
                fprintf(stderr, "\nERROR: Invalid digest specificity: %d\n\t>%s", params.iDigestSpecificity, sLine);
                tempest_exit(EXIT_FAILURE); 
            }
        }

        else if (strcmp(sParam, "digest_missed_cleavages") == 0) {
            if (sscanf(sLine, "digest_missed_cleavages %d \n", &params.iDigestMaxMissedCleavages) != 1) {
                fprintf(stderr, "\nERROR: Could not parse max missed cleavages\n\t>%s", sLine);
                tempest_exit(EXIT_FAILURE);
            }

            if (params.iDigestMaxMissedCleavages < 0) {
                fprintf(stderr, "\nERROR: Invalid max missed cleavages: %d\n\t>%s", params.iDigestMaxMissedCleavages, sLine);
                tempest_exit(EXIT_FAILURE); 
            }
        }

        else if (strcmp(sParam, "digest_length") == 0) {
            if (sscanf(sLine, "digest_length %d - %d \n", &params.iMinPeptideLength, &params.iMaxPeptideLength) != 2) {
                fprintf(stderr, "\nERROR: Could not parse digest length range\n\t>%s", sLine);
                tempest_exit(EXIT_FAILURE);
            }

            if (params.iMinPeptideLength < MIN_PEPTIDE_LENGTH) {
                fprintf(stderr, "\nERROR: Minimum digest length must be greater than %d\n\t>%s", MIN_PEPTIDE_LENGTH, sLine);
                tempest_exit(EXIT_FAILURE);
            }

            if (params.iMaxPeptideLength > MAX_PEPTIDE_LENGTH) {
                fprintf(stderr, "\nERROR: Maximum digest length must be less than %d\n\t>%s", MAX_PEPTIDE_LENGTH, sLine);
                tempest_exit(EXIT_FAILURE); 
            }

            if (params.iMinPeptideLength > params.iMaxPeptideLength) {
                fprintf(stderr, "\nERROR: Invalid digest length range (min must be less than max)\n\t>%s", sLine);
                tempest_exit(EXIT_FAILURE); 
            }
        }

        else if (strcmp(sParam, "fixed_modification") == 0) {
            if (sscanf(sLine, "fixed_modification %s %lf \n", sLocation, &dMassDiff) != 2) {
                fprintf(stderr, "\nERROR: Could not parse fixed modification\n\t>%s", sLine);
                tempest_exit(EXIT_FAILURE);
            }

            //validate location and update masses
            if (strcmp(sLocation, "nterm") == 0 || strcmp(sLocation, "peptide-nterm") == 0) {
                params.dPeptideNtermMass = dMassDiff;
            }
            else if (strcmp(sLocation, "cterm") == 0 || strcmp(sLocation, "peptide-cterm") == 0) {
                params.dPeptideCtermMass = dMassDiff;
            }
            else if (strcmp(sLocation, "protein-nterm") == 0) {
                params.dProteinNtermMass = dMassDiff;
            }
            else if (strcmp(sLocation, "protein-cterm") == 0) {
                params.dProteinCtermMass = dMassDiff;
            }
            else {
                for (c=sLocation; *c; c++) {
                    if (!isupper(*c)) {
                        fprintf(stderr, "\nERROR: Invalid modification location: %s\n\t>%s", sLocation, sLine);
                        tempest_exit(EXIT_FAILURE);
                    }

                    dMassAA[aa2i(*c)] += dMassDiff;
                    dMassAA[aa2i(tolower(*c))] += dMassDiff;
                    append_mod(*c, dMassDiff, '\0');
                }
            }
        }

        else if (strcmp(sParam, "variable_modification") == 0) {
            // parse with or without NL (default 0.0)
            if (sscanf(sLine, "variable_modification %s %c %lf \n", sLocation, &cSymbol, &dMassDiff) != 3) {
                fprintf(stderr, "\nERROR: Could not parse variable modification\n\t>%s", sLine);
                tempest_exit(EXIT_FAILURE);
            }

            //validate symbol
            if (isalpha(cSymbol)) {
                fprintf(stderr, "\nERROR: Invalid modification symbol: %c\n\t>%s", cSymbol, sLine);
                tempest_exit(EXIT_FAILURE);
            }
            else if (!isprint(cSymbol)) {
                fprintf(stderr, "\nERROR: Invalid modification symbol (nonprintable): ascii %d\n\t>%s", (short) cSymbol, sLine);
                tempest_exit(EXIT_FAILURE);
            }

            //validate location and update masses/symbols
            if (strcmp(sLocation, "nterm") == 0 || strcmp(sLocation, "peptide-nterm") == 0) {
                params.dVariableNtermMassDiff = dMassDiff;
                params.cVariableNtermSymbol = cSymbol;
            }
            else if (strcmp(sLocation, "cterm") == 0 || strcmp(sLocation, "peptide-cterm") == 0) {
                params.dVariableCtermMassDiff = dMassDiff;
                params.cVariableCtermSymbol = cSymbol;
            }
            else {
                for (c=sLocation; *c; c++) {
                    if (!isupper(*c)) {
                        fprintf(stderr, "\nERROR: Invalid modification location: %s\n\t>%s", sLocation, sLine);
                        tempest_exit(EXIT_FAILURE);
                    }

                    if (cModSites[aa2i(*c)]) {
                        fprintf(stderr, "\nERROR: Multiple mods specificied for %c\n", *c);
                        tempest_exit(EXIT_FAILURE);
                    }

                    cModSites[aa2i(*c)] = cSymbol;
                    fModValues[aa2i(*c)] = dMassDiff;
                    dMassAA[aa2i(tolower(*c))] += dMassDiff;
                    append_mod(*c, dMassDiff, cSymbol);
                }
            }
        }

        else if (strcmp(sParam, "variable_modifications_max") == 0) {
            if (sscanf(sLine, "variable_modifications_max %d \n", &params.iModificationsMax) != 1) {
                fprintf(stderr, "\nERROR: Could not parse max variable modifications\n\t>%s", sLine);
                tempest_exit(EXIT_FAILURE);
            }

            if (params.iModificationsMax < 0) {
                fprintf(stderr, "\nERROR: Invalid max variable modifications: %d\n\t>%s", params.iModificationsMax, sLine);
                tempest_exit(EXIT_FAILURE); 
            }
        }

        else if (strcmp(sParam, "precursor_tolerance") == 0) {
            if (sscanf(sLine, "precursor_tolerance %f %s \n", &params.fPrecursorTolerance, sTemp) != 2) {
                fprintf(stderr, "\nERROR: Could not parse precursor tolerance\n\t>%s", sLine);
                tempest_exit(EXIT_FAILURE);
            }

            if (params.fPrecursorTolerance <= 0) {
                fprintf(stderr, "\nERROR: Invalid precursor tolerance: %f\n\t>%s", params.fPrecursorTolerance, sLine);
                tempest_exit(EXIT_FAILURE); 
            }

            if (strcmp(sTemp, "DA") == 0 || strcmp(sTemp, "Da") == 0 || strcmp(sTemp, "da") == 0) {
                params.bPrecursorTolerancePPM = 0;
            }
            else if (strcmp(sTemp, "PPM") == 0 || strcmp(sTemp, "ppm") == 0) {
                params.bPrecursorTolerancePPM = 1;
            }
            else {
                fprintf(stderr, "\nERROR: Invalid precursor tolerance units: %s\n\t>%s", sTemp, sLine);
                tempest_exit(EXIT_FAILURE);
            }
        }

        else if (strcmp(sParam, "neutral_loss") == 0) {
            // parse with or without NL (default 0.0)
            double dWeighting;
            if (sscanf(sLine, "neutral_loss %s %lf %lf \n", sLocation, &dMassDiff, &dWeighting) != 3) {
                fprintf(stderr, "\nERROR: Could not parse neutral loss\n\t>%s", sLine);
                tempest_exit(EXIT_FAILURE);
            }

            //validate location and update masses/symbols
            if (strcmp(sLocation, "nterm") == 0 || strcmp(sLocation, "peptide-nterm") == 0) {
                params.numNtermNL += 1;
                nlValue nlv;
                nlv.massDelta = (float)dMassDiff;
                nlv.weighting = (float)dWeighting;
                nlValuesNterm.push_back(nlv);
            }
            else if (strcmp(sLocation, "cterm") == 0 || strcmp(sLocation, "peptide-cterm") == 0) {
                params.numCtermNL += 1;
                nlValue nlv;
                nlv.massDelta = (float)dMassDiff;
                nlv.weighting = (float)dWeighting;
                nlValuesCterm.push_back(nlv);
            }
            else {
                bool hasAA['z'-'A'+1] = {0};
                int numSites = 0;
                for (c=sLocation; *c; c++) {
                    if (!isalnum(*c)) {
                        bool matchedMod;
                        for (int i=0; i<params.iNumMods; i++) {
                            if (params.tMods[i].cSymbol == *c) {
                                hasAA[aa2i(tolower(params.tMods[i].cAminoAcid))] = 1;
                                numSites += 1;
                                matchedMod = true;
                            }
                        }
                        if (!matchedMod)
                            printf("\nWARNING: Neutral loss modification symbol (%c) is not used in any specified variable mods.\n", *c);
                    }
                    else if (isalpha(*c) && isupper(*c)) {
                        hasAA[aa2i(*c)] = 1;
                        numSites += 1;
                    }                  
                    else {
                        fprintf(stderr, "\nERROR: Neutral loss location (%c) is not a valid amino acid or modification symbol\n", *c);
                        tempest_exit(EXIT_FAILURE);
                    } 
                }
                if (numSites) {
                    params.numAANL += 1;
                    nlValue nlv;
                    nlv.numSites = numSites;
                    nlv.massDelta = dMassDiff;
                    nlv.weighting = dWeighting;
                    memcpy(nlv.hasAA, hasAA, ('z'-'A'+1)*sizeof(bool));
                    nlValuesAA.push_back(nlv);
                }
            }
        }

        else if (strcmp(sParam, "msms_fragment_tolerance") == 0) {
            if (sscanf(sLine, "msms_fragment_tolerance %f %s \n", &params.fFragmentTolerance, sTemp) != 2) {
                fprintf(stderr, "\nERROR: Could not parse MS/MS fragment tolerance\n\t>%s", sLine);
                tempest_exit(EXIT_FAILURE);
            }

            if (params.fFragmentTolerance <= 0) {
                fprintf(stderr, "\nERROR: Invalid MS/MS fragment tolerance: %f\n\t>%s", params.fFragmentTolerance, sLine);
                tempest_exit(EXIT_FAILURE); 
            }

            if (strcmp(sTemp, "MZ") == 0 || strcmp(sTemp, "mz") == 0 || strcmp(sTemp, "M/Z") == 0 || strcmp(sTemp, "m/z") == 0) {
                params.bFragmentTolerancePPM = 0;
            }
            else if (strcmp(sTemp, "PPM") == 0 || strcmp(sTemp, "ppm") == 0) {
                fprintf(stderr, "\nPPM Fragment Tolerance is unsupported\n");
                tempest_exit(EXIT_FAILURE);
                params.bFragmentTolerancePPM = 1;
            }
            else {
                fprintf(stderr, "\nERROR: Invalid MS/MS fragment tolerance units: %s\n\t>%s", sTemp, sLine);
                tempest_exit(EXIT_FAILURE);
            }
        }

        else if (strcmp(sParam, "msms_remove_precursor") == 0) {
            if (sscanf(sLine, "msms_remove_precursor %d \n", &iTemp) != 1) {
                fprintf(stderr, "\nERROR: Could not parse MS/MS remove precursor flag\n\t>%s", sLine);
                tempest_exit(EXIT_FAILURE);
            }

            if (iTemp == 1) {
                params.bRemovePrecursors = 1;
            }
            else if (iTemp != 0) {
                fprintf(stderr, "\nERROR: Invalid flag for MS/MS remove precursor (expected 0 or 1): %d\n\t>%s", iTemp, sLine);
                tempest_exit(EXIT_FAILURE); 
            }
        }

        else if (strcmp(sParam, "msms_remove_mz_range") == 0) {
            params.iNumRemoveMzRanges += 1;
            if (0 == (params.fRemoveMzRanges = (float*) realloc(params.fRemoveMzRanges, 2*params.iNumRemoveMzRanges*sizeof(float)))) {
                fprintf(stderr, "\nERROR: Could not allocate memory for %d remove M/Z ranges.\n", params.iNumRemoveMzRanges);
                tempest_exit(EXIT_FAILURE);
            }

            fMzRange = params.fRemoveMzRanges+2*(params.iNumRemoveMzRanges-1);
            if (sscanf(sLine, "msms_remove_mz_range %f - %f \n", fMzRange, fMzRange+1) != 2) {
                fprintf(stderr, "\nERROR: Could not parse remove m/z range\n\t>%s", sLine);
                tempest_exit(EXIT_FAILURE);
            }

            if (fMzRange[0] < 0.0f) {
                fprintf(stderr, "\nERROR: Invalid remove m/z range (min cannot be less than 0)\n\t>%s", sLine);
                tempest_exit(EXIT_FAILURE); 
            }

            if (fMzRange[0] > fMzRange[1]) {
                fprintf(stderr, "\nERROR: Invalid remove m/z range (min cannot be greater than max)\n\t>%s", sLine);
                tempest_exit(EXIT_FAILURE); 
            }
        }

        else if (strcmp(sParam, "xcorr_transform") == 0) {
            if (sscanf(sLine, "xcorr_transform %d \n", &iTemp) != 1) {
                fprintf(stderr, "\nERROR: Could not parse xcorr transform flag\n\t>%s", sLine);
                tempest_exit(EXIT_FAILURE);
            }

            if (iTemp == 1) {
                params.bCrossCorrelation = 1;
            }

            else if (iTemp == 0) {
                params.bCrossCorrelation = 0;
            }

            else {
                fprintf(stderr, "\nERROR: Invalid flag for xcorr_transform (expected 0 or 1): %d\n\t>%s", iTemp, sLine);
                tempest_exit(EXIT_FAILURE); 
            }
        }

        else if (strcmp(sParam, "theoretical_flanking") == 0) {
            if (sscanf(sLine, "theoretical_flanking %d \n", &iTemp) != 1) {
                fprintf(stderr, "\nERROR: Could not parse theoretical_flanking flag\n\t>%s", sLine);
                tempest_exit(EXIT_FAILURE);
            }

            if (iTemp == 0) {
                params.bTheoreticalFlanking = 0;
            }
            else if (iTemp != 1) {
                fprintf(stderr, "\nERROR: Invalid flag for theoretical_flanking (expected 0 or 1): %d\n\t>%s", iTemp, sLine);
                tempest_exit(EXIT_FAILURE); 
            }
        }

        else if (strcmp(sParam, "num_output_psms") == 0) {
            if (sscanf(sLine, "num_output_psms %d \n", &params.iNumOutputPSMs) != 1) {
                fprintf(stderr, "\nERROR: Could not parse num output PSMs\n\t>%s", sLine);
                tempest_exit(EXIT_FAILURE);
            }

            if (params.iNumOutputPSMs < 1) {
                fprintf(stderr, "\nERROR: Invalid num output PSMs: %d\n\t>%s", params.iNumOutputPSMs, sLine);
                tempest_exit(EXIT_FAILURE); 
            }
        }

        else if (strcmp(sParam, "fix_delta_score") == 0) {
            if (sscanf(sLine, "fix_delta_score %d \n", &iTemp) != 1) {
                fprintf(stderr, "\nERROR: Could not parse fix delta score flag\n\t>%s", sLine);
                tempest_exit(EXIT_FAILURE);
            }

            if (iTemp == 1) {
                params.bFixDeltaScore = 1;
            }
            else if (iTemp != 0) {
                fprintf(stderr, "\nERROR: Invalid flag for fix delta score (expected 0 or 1): %d\n\t>%s", iTemp, sLine);
                tempest_exit(EXIT_FAILURE); 
            }
        }

        else {
            fprintf(stderr, "\nWARNING: Unknown search paramater: '%s'\n", sParam);
            tempest_exit(EXIT_FAILURE);
        }
    }
}

void append_mod(char cAminoAcid, double dMassDiff, char cSymbol) {
    params.iNumMods += 1;
    if (0 == (params.tMods = (mod_t*) realloc(params.tMods, params.iNumMods*sizeof(mod_t)))) {
        fprintf(stderr, "\nERROR: Could not allocate memory for %d mods.\n", params.iNumMods);
        tempest_exit(EXIT_FAILURE);
    }

    params.tMods[params.iNumMods-1].cAminoAcid = cAminoAcid;
    params.tMods[params.iNumMods-1].cSymbol = cSymbol;
    params.tMods[params.iNumMods-1].dMassDiff = dMassDiff;
}
