#include "tempest.hpp"

void append_mod(char, int, double, char);

extern void Tempest::parse_params() {
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

    if((fp = (FILE*) fopen(Tempest::args.sParams, "r")) == NULL) {
        fprintf( stderr, "ERROR\tUnable to open params file %s: %s\n", Tempest::args.sParams, strerror(errno) );
        Tempest::tempest_exit(EXIT_FAILURE);
    }
    
    // read params line by line
    while (fgets(sLine, sizeof(sLine), fp)) {
        if (sscanf(sLine, "%s", sParam) == -1) continue;
        if (sParam[0] == '#') continue;

        if (strcmp(sParam, "digest_definition") == 0) {
            Tempest::params.sDigestSites = (char*) malloc(26*sizeof(char));
            Tempest::params.sDigestNoSites = (char*) malloc(26*sizeof(char));
            if (sscanf(sLine, "digest_definition %s %s %d \n", Tempest::params.sDigestSites, Tempest::params.sDigestNoSites, &Tempest::params.iDigestOffset) != 3) {
                fprintf(stderr, "\nERROR: Could not parse digest definition\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }
        }

        else if (strcmp(sParam, "digest_specificity") == 0) {
            if (sscanf(sLine, "digest_specificity %d \n", &Tempest::params.iDigestSpecificity) != 1) {
                fprintf(stderr, "\nERROR: Could not parse digest specificity\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }

            if (Tempest::params.iDigestSpecificity < 0 || Tempest::params.iDigestSpecificity > 2) {
                fprintf(stderr, "\nERROR: Invalid digest specificity: %d\n\t>%s", Tempest::params.iDigestSpecificity, sLine);
                Tempest::tempest_exit(EXIT_FAILURE); 
            }
        }

        else if (strcmp(sParam, "digest_missed_cleavages") == 0) {
            if (sscanf(sLine, "digest_missed_cleavages %d \n", &Tempest::params.iDigestMaxMissedCleavages) != 1) {
                fprintf(stderr, "\nERROR: Could not parse max missed cleavages\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }

            if (Tempest::params.iDigestMaxMissedCleavages < 0) {
                fprintf(stderr, "\nERROR: Invalid max missed cleavages: %d\n\t>%s", Tempest::params.iDigestMaxMissedCleavages, sLine);
                Tempest::tempest_exit(EXIT_FAILURE); 
            }
        }

        else if (strcmp(sParam, "digest_length") == 0) {
            if (sscanf(sLine, "digest_length %d %d \n", &Tempest::params.iMinPeptideLength, &Tempest::params.iMaxPeptideLength) != 2) {
                fprintf(stderr, "\nERROR: Could not parse digest length range\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }

            if (Tempest::params.iMinPeptideLength < MIN_PEPTIDE_LENGTH) {
                fprintf(stderr, "\nERROR: Minimum digest length must be greater than %d\n\t>%s", MIN_PEPTIDE_LENGTH, sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }

            if (Tempest::params.iMaxPeptideLength > MAX_PEPTIDE_LENGTH) {
                fprintf(stderr, "\nERROR: Maximum digest length must be less than %d\n\t>%s", MAX_PEPTIDE_LENGTH, sLine);
                Tempest::tempest_exit(EXIT_FAILURE); 
            }

            if (Tempest::params.iMinPeptideLength > Tempest::params.iMaxPeptideLength) {
                fprintf(stderr, "\nERROR: Invalid digest length range (min must be less than max)\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE); 
            }
        }

        else if (strcmp(sParam, "digest_mass_range") == 0) {
            if (sscanf(sLine, "digest_mass_range %f %f \n", &Tempest::params.minPeptideMass, &Tempest::params.maxPeptideMass) != 2) {
                fprintf(stderr, "\nERROR: Could not parse digest_mass_range\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }

            if (Tempest::params.minPeptideMass < 0 || Tempest::params.maxPeptideMass < 0) {
                fprintf(stderr, "\nERROR: Negative value(s) in digest_mass_range:\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE); 
            }

            if (Tempest::params.minPeptideMass > Tempest::params.maxPeptideMass) {
                fprintf(stderr, "\nERROR: Min peptide mass greater than max peptide mass in digest_mass_range:\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE); 
            }
        }

        else if (strcmp(sParam, "fixed_modification") == 0) {
            if (sscanf(sLine, "fixed_modification %s %lf \n", sLocation, &dMassDiff) != 2) {
                fprintf(stderr, "\nERROR: Could not parse fixed modification\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }

            //validate location and update masses
            if (strcmp(sLocation, "nterm") == 0 || strcmp(sLocation, "peptide-nterm") == 0) {
                Tempest::params.dPeptideNtermMass = dMassDiff;
            }
            else if (strcmp(sLocation, "cterm") == 0 || strcmp(sLocation, "peptide-cterm") == 0) {
                Tempest::params.dPeptideCtermMass = dMassDiff;
            }
            else if (strcmp(sLocation, "protein-nterm") == 0) {
                Tempest::params.dProteinNtermMass = dMassDiff;
            }
            else if (strcmp(sLocation, "protein-cterm") == 0) {
                Tempest::params.dProteinCtermMass = dMassDiff;
            }
            else {
                for (c=sLocation; *c; c++) {
                    if (!isupper(*c)) {
                        fprintf(stderr, "\nERROR: Invalid modification location: %s\n\t>%s", sLocation, sLine);
                        Tempest::tempest_exit(EXIT_FAILURE);
                    }
                    for (int offset=0; ('Z'+offset)<256; offset+=32) {
                        Tempest::dMassAA[*c+offset] += dMassDiff;
                    }
                    append_mod(*c, -1, dMassDiff, '\0');
                }
            }
        }

        else if (strcmp(sParam, "variable_modification") == 0) {
            // parse with or without NL (default 0.0)
            if (sscanf(sLine, "variable_modification %s %c %lf \n", sLocation, &cSymbol, &dMassDiff) != 3) {
                fprintf(stderr, "\nERROR: Could not parse variable modification\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }

            //validate symbol
            if (isalpha(cSymbol)) {
                fprintf(stderr, "\nERROR: Invalid modification symbol: %c\n\t>%s", cSymbol, sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }
            else if (!isprint(cSymbol)) {
                fprintf(stderr, "\nERROR: Invalid modification symbol (nonprintable): ascii %d\n\t>%s", (short) cSymbol, sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }

            //validate location and update masses/symbols
            if (strcmp(sLocation, "nterm") == 0 || strcmp(sLocation, "peptide-nterm") == 0) {
                if (Tempest::numNtermModSites >= 5) {
                    fprintf(stderr, "\nERROR: Too many mods specified for N-terminus (max 5)\n");
                    Tempest::tempest_exit(EXIT_FAILURE);
                }
                Tempest::ntermModMasses[Tempest::numNtermModSites] = dMassDiff;
                Tempest::ntermModSymbols[Tempest::numNtermModSites] = cSymbol;
                Tempest::numNtermModSites += 1;
                
            }
            else if (strcmp(sLocation, "cterm") == 0 || strcmp(sLocation, "peptide-cterm") == 0) {
                if (Tempest::numCtermModSites >= 5) {
                    fprintf(stderr, "\nERROR: Too many mods specified for C-terminus (max 5)\n");
                    Tempest::tempest_exit(EXIT_FAILURE);
                }
                Tempest::ctermModMasses[Tempest::numCtermModSites] = dMassDiff;
                Tempest::ctermModSymbols[Tempest::numCtermModSites] = cSymbol;
                Tempest::numCtermModSites += 1;
            }
            else {
                for (c=sLocation; *c; c++) {
                    if (!isupper(*c)) {
                        fprintf(stderr, "\nERROR: Invalid modification location: %s\n\t>%s", sLocation, sLine);
                        Tempest::tempest_exit(EXIT_FAILURE);
                    }

                    if (Tempest::numAAModSites[*c] >= 5) {
                        fprintf(stderr, "\nERROR: Too many mods specified for %c (max 5)\n", *c);
                        Tempest::tempest_exit(EXIT_FAILURE);
                    }

                    int modInd = Tempest::numAAModSites[*c];
                    Tempest::cModSites[*c][modInd] = cSymbol;
                    Tempest::fModValues[*c][modInd] = dMassDiff;
                    Tempest::dMassAA[Tempest::toMod(*c, modInd+1)] += dMassDiff;
                    append_mod(*c, modInd, dMassDiff, cSymbol);
                    Tempest::numAAModSites[*c] += 1;
                }
            }
        }

        else if (strcmp(sParam, "variable_modifications_max") == 0) {
            if (sscanf(sLine, "variable_modifications_max %d \n", &Tempest::params.iModificationsMax) != 1) {
                fprintf(stderr, "\nERROR: Could not parse max variable modifications\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }

            if (Tempest::params.iModificationsMax < 0) {
                fprintf(stderr, "\nERROR: Invalid max variable modifications: %d\n\t>%s", Tempest::params.iModificationsMax, sLine);
                Tempest::tempest_exit(EXIT_FAILURE); 
            }
        }

        else if (strcmp(sParam, "precursor_tolerance") == 0) {
            if (sscanf(sLine, "precursor_tolerance %f %s \n", &Tempest::params.fPrecursorTolerance, sTemp) != 2) {
                fprintf(stderr, "\nERROR: Could not parse precursor tolerance\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }

            if (Tempest::params.fPrecursorTolerance <= 0) {
                fprintf(stderr, "\nERROR: Invalid precursor tolerance: %f\n\t>%s", Tempest::params.fPrecursorTolerance, sLine);
                Tempest::tempest_exit(EXIT_FAILURE); 
            }

            if (strcmp(sTemp, "DA") == 0 || strcmp(sTemp, "Da") == 0 || strcmp(sTemp, "da") == 0) {
                Tempest::params.bPrecursorTolerancePPM = 0;
            }
            else if (strcmp(sTemp, "PPM") == 0 || strcmp(sTemp, "ppm") == 0) {
                Tempest::params.bPrecursorTolerancePPM = 1;
            }
            else {
                fprintf(stderr, "\nERROR: Invalid precursor tolerance units: %s\n\t>%s", sTemp, sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }
        }

        else if (strcmp(sParam, "neutral_loss") == 0) {
            // parse with or without NL (default 0.0)
            double dWeighting;
            if (sscanf(sLine, "neutral_loss %s %lf %lf \n", sLocation, &dMassDiff, &dWeighting) != 3) {
                fprintf(stderr, "\nERROR: Could not parse neutral loss\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }

            //validate location and update masses/symbols
            if (strcmp(sLocation, "nterm") == 0 || strcmp(sLocation, "peptide-nterm") == 0) {
                Tempest::params.numNtermNL += 1;
                nlValue nlv;
                nlv.massDelta = (float)dMassDiff;
                nlv.weighting = (float)dWeighting;
                nlv.modNum = 0;
                Tempest::nlValuesNterm.push_back(nlv);
            }
            else if (strcmp(sLocation, "cterm") == 0 || strcmp(sLocation, "peptide-cterm") == 0) {
                Tempest::params.numCtermNL += 1;
                nlValue nlv;
                nlv.massDelta = (float)dMassDiff;
                nlv.weighting = (float)dWeighting;
                nlv.modNum = 0;
                Tempest::nlValuesCterm.push_back(nlv);
            }
            else {
                bool hasAA[256] = {0};
                int numSites = 0;
                for (c=sLocation; *c; c++) {
                    if (!isalnum(*c)) {
                        bool matchedMod;
                        //MEA: cases exclusive?
                        for (int i=0; i<Tempest::numNtermModSites; i++) {
                            if (*c == Tempest::ntermModSymbols[i]) {
                                Tempest::params.numNtermNL += 1;
                                nlValue nlv;
                                nlv.massDelta = (float)dMassDiff;
                                nlv.weighting = (float)dWeighting;
                                nlv.modNum = i+1;
                                Tempest::nlValuesNterm.push_back(nlv);
                                matchedMod = true;
                            }
                        }
                        for (int i=0; i<Tempest::numCtermModSites; i++) {
                            if (*c == Tempest::ctermModSymbols[i]) {
                                Tempest::params.numCtermNL += 1;
                                nlValue nlv;
                                nlv.massDelta = (float)dMassDiff;
                                nlv.weighting = (float)dWeighting;
                                nlv.modNum = i+1;
                                Tempest::nlValuesCterm.push_back(nlv);
                                matchedMod = true;
                            }
                        }
                        for (int i=0; i<Tempest::params.iNumMods; i++) {
                            if (Tempest::params.tMods[i].cSymbol == *c) {
                                hasAA[Tempest::toMod(Tempest::params.tMods[i].cAminoAcid, Tempest::params.tMods[i].modInd+1)] = 1;
                                numSites += 1;
                                matchedMod = true;
                            }
                        }
                        if (!matchedMod) {
                            printf("\nWARNING: Neutral loss modification symbol (%c) is not used in any specified variable mods.\n", *c);
                            printf("         Note: specify variable mod neutral losses after variable_mod definition.\n");
                        }
                    }
                    else if (isalpha(*c) && isupper(*c)) {
                        hasAA[*c] = 1;
                        numSites += 1;
                    }
                    else {
                        fprintf(stderr, "\nERROR: Neutral loss location (%c) is not a valid amino acid or modification symbol\n", *c);
                        Tempest::tempest_exit(EXIT_FAILURE);
                    } 
                }
                if (numSites) {
                    Tempest::params.numAANL += 1;
                    nlValue nlv;
                    nlv.numSites = numSites;
                    nlv.massDelta = dMassDiff;
                    nlv.weighting = dWeighting;
                    memcpy(nlv.hasAA, hasAA, 256*sizeof(bool));
                    Tempest::nlValuesAA.push_back(nlv);
                }
            }
        }

        else if (strcmp(sParam, "use_a_ions") == 0) {
            if (sscanf(sLine, "use_a_ions %d \n", &iTemp) != 1) {
                fprintf(stderr, "\nERROR: Could not parse use_a_ions flag:\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }
            if (iTemp != 0 && iTemp != 1) {
                fprintf(stderr, "\nERROR: Invalid flag for use_a_ions (expected 0 or 1):\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }
            Tempest::params.useAIons = iTemp;
        }
        else if (strcmp(sParam, "use_b_ions") == 0) {
            if (sscanf(sLine, "use_b_ions %d \n", &iTemp) != 1) {
                fprintf(stderr, "\nERROR: Could not parse use_b_ions flag:\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }
            if (iTemp != 0 && iTemp != 1) {
                fprintf(stderr, "\nERROR: Invalid flag for use_b_ions (expected 0 or 1):\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }
            Tempest::params.useBIons = iTemp;
        }
        else if (strcmp(sParam, "use_c_ions") == 0) {
            if (sscanf(sLine, "use_c_ions %d \n", &iTemp) != 1) {
                fprintf(stderr, "\nERROR: Could not parse use_c_ions flag:\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }
            if (iTemp != 0 && iTemp != 1) {
                fprintf(stderr, "\nERROR: Invalid flag for use_c_ions (expected 0 or 1):\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }
            Tempest::params.useCIons = iTemp;
        }
        else if (strcmp(sParam, "use_x_ions") == 0) {
            if (sscanf(sLine, "use_x_ions %d \n", &iTemp) != 1) {
                fprintf(stderr, "\nERROR: Could not parse use_x_ions flag:\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }
            if (iTemp != 0 && iTemp != 1) {
                fprintf(stderr, "\nERROR: Invalid flag for use_x_ions (expected 0 or 1):\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }
            Tempest::params.useXIons = iTemp;
        }
        else if (strcmp(sParam, "use_y_ions") == 0) {
            if (sscanf(sLine, "use_y_ions %d \n", &iTemp) != 1) {
                fprintf(stderr, "\nERROR: Could not parse use_y_ions flag:\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }
            if (iTemp != 0 && iTemp != 1) {
                fprintf(stderr, "\nERROR: Invalid flag for use_y_ions (expected 0 or 1):\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }
            Tempest::params.useYIons = iTemp;
        }
        else if (strcmp(sParam, "use_z_ions") == 0) {
            if (sscanf(sLine, "use_z_ions %d \n", &iTemp) != 1) {
                fprintf(stderr, "\nERROR: Could not parse use_z_ions flag:\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }
            if (iTemp != 0 && iTemp != 1) {
                fprintf(stderr, "\nERROR: Invalid flag for use_z_ions (expected 0 or 1):\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }
            Tempest::params.useZIons = iTemp;
        }

        else if (strcmp(sParam, "fragment_tolerance") == 0) {
            if (sscanf(sLine, "fragment_tolerance %f %s \n", &Tempest::params.fFragmentTolerance, sTemp) != 2) {
                fprintf(stderr, "\nERROR: Could not parse MS/MS fragment tolerance\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }

            if (Tempest::params.fFragmentTolerance <= 0) {
                fprintf(stderr, "\nERROR: Invalid MS/MS fragment tolerance: %f\n\t>%s", Tempest::params.fFragmentTolerance, sLine);
                Tempest::tempest_exit(EXIT_FAILURE); 
            }

            if (strcmp(sTemp, "MZ") == 0 || strcmp(sTemp, "mz") == 0 || strcmp(sTemp, "M/Z") == 0 || strcmp(sTemp, "m/z") == 0) {
                Tempest::params.bFragmentTolerancePPM = 0;
            }
            else if (strcmp(sTemp, "PPM") == 0 || strcmp(sTemp, "ppm") == 0) {
                fprintf(stderr, "\nPPM Fragment Tolerance is unsupported\n");
                Tempest::tempest_exit(EXIT_FAILURE);
                Tempest::params.bFragmentTolerancePPM = 1;
            }
            else {
                fprintf(stderr, "\nERROR: Invalid MS/MS fragment tolerance units: %s\n\t>%s", sTemp, sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }
        }

        else if (strcmp(sParam, "fragment_bin_offset") == 0) {
            if (sscanf(sLine, "fragment_bin_offset %f \n", &Tempest::params.fragmentBinOffset) != 1) {
                fprintf(stderr, "\nERROR: Could not parse fragment_bin_offset:\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }
            if (Tempest::params.fragmentBinOffset < 0 || Tempest::params.fragmentBinOffset >= 1) {
                fprintf(stderr, "\nERROR: Invalid value for fragment_bin_offset (expected value between 0 and 1):\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE); 
            }
        }

        else if (strcmp(sParam, "msms_remove_precursor") == 0) {
            if (sscanf(sLine, "msms_remove_precursor %d \n", &iTemp) != 1) {
                fprintf(stderr, "\nERROR: Could not parse MS/MS remove precursor flag\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }

            if (iTemp == 1) {
                Tempest::params.bRemovePrecursors = 1;
            }
            else if (iTemp != 0) {
                fprintf(stderr, "\nERROR: Invalid flag for MS/MS remove precursor (expected 0 or 1): %d\n\t>%s", iTemp, sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }
        }

        else if (strcmp(sParam, "msms_remove_mz_range") == 0) {
            Tempest::params.iNumRemoveMzRanges += 1;
            if (0 == (Tempest::params.fRemoveMzRanges = (float*) realloc(Tempest::params.fRemoveMzRanges, 2*Tempest::params.iNumRemoveMzRanges*sizeof(float)))) {
                fprintf(stderr, "\nERROR: Could not allocate memory for %d remove M/Z ranges.\n", Tempest::params.iNumRemoveMzRanges);
                Tempest::tempest_exit(EXIT_FAILURE);
            }

            fMzRange = Tempest::params.fRemoveMzRanges+2*(Tempest::params.iNumRemoveMzRanges-1);
            if (sscanf(sLine, "msms_remove_mz_range %f %f \n", fMzRange, fMzRange+1) != 2) {
                fprintf(stderr, "\nERROR: Could not parse remove m/z range\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }

            if (fMzRange[0] < 0.0f) {
                fprintf(stderr, "\nERROR: Invalid remove m/z range (min cannot be less than 0):\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE); 
            }

            if (fMzRange[0] > fMzRange[1]) {
                fprintf(stderr, "\nERROR: Invalid remove m/z range (min cannot be greater than max):\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE); 
            }
        }

        else if (strcmp(sParam, "normalization_regions") == 0) {
            if (sscanf(sLine, "normalization_regions %d \n", &iTemp) != 1) {
                fprintf(stderr, "\nERROR: Could not parse normalization_regions:\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }
            if (iTemp > 0)
                Tempest::params.numNormRegions = iTemp;
            else {
                fprintf(stderr, "\nERROR: Invalid value for normalization_regions (expected value >0):\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE); 
            }
        }

        else if (strcmp(sParam, "intensity_threshold") == 0) {
            if (sscanf(sLine, "intensity_threshold %f \n", &Tempest::params.intensityThreshold) != 1) {
                fprintf(stderr, "\nERROR: Could not parse intensity_threshold:\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }
            if (Tempest::params.intensityThreshold < 0 || Tempest::params.intensityThreshold > 1) {
                fprintf(stderr, "\nERROR: Invalid value for intensity_threshold (expected value between 0 and 1):\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE); 
            }
        }

        else if (strcmp(sParam, "xcorr_transform_width") == 0) {
            if (sscanf(sLine, "xcorr_transform_width %d \n", &Tempest::params.xcorrTransformWidth) != 1) {
                fprintf(stderr, "\nERROR: Could not parse xcorr_transform_width:\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }
            if (Tempest::params.xcorrTransformWidth < 0) {
                fprintf(stderr, "\nERROR: Expected non-negative value for xcorr_transform_width :\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE); 
            }
        }

        else if (strcmp(sParam, "flanking_intensity") == 0) {
            if (sscanf(sLine, "flanking_intensity %f \n", &Tempest::params.flankingIntensity) != 1) {
                fprintf(stderr, "\nERROR: Could not parse flanking_intensity value:\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }

            else if (Tempest::params.flankingIntensity < 0 || Tempest::params.flankingIntensity >= 1) {
                fprintf(stderr, "\nERROR: Invalid value for flanking_intensity (expected value between 0 and 1):\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE); 
            }
        }

        else if (strcmp(sParam, "num_internal_psms") == 0) {
            if (sscanf(sLine, "num_internal_psms %d \n", &Tempest::params.numInternalPSMs) != 1) {
                fprintf(stderr, "\nERROR: Could not parse num_internal_psms:\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }

            if (Tempest::params.numInternalPSMs < 2) {
                fprintf(stderr, "\nERROR: Invalid num_internal_psms: %d\n\t>%s", Tempest::params.numInternalPSMs, sLine);
                Tempest::tempest_exit(EXIT_FAILURE); 
            }
        }

        else if (strcmp(sParam, "num_output_psms") == 0) {
            if (sscanf(sLine, "num_output_psms %d \n", &Tempest::params.numOutputPSMs) != 1) {
                fprintf(stderr, "\nERROR: Could not parse num_output_psms:\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }

            if (Tempest::params.numOutputPSMs < 1) {
                fprintf(stderr, "\nERROR: Invalid num_output_psms: %d\n\t>%s", Tempest::params.numOutputPSMs, sLine);
                Tempest::tempest_exit(EXIT_FAILURE); 
            }
        }

        else if (strcmp(sParam, "backbone_delta_score") == 0) {
            if (sscanf(sLine, "backbone_delta_score %d \n", &iTemp) != 1) {
                fprintf(stderr, "\nERROR: Could not parse backbone_delta_score flag\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }

            if (iTemp == 1) {
                Tempest::params.bFixDeltaScore = 1;
            }
            else if (iTemp != 0) {
                fprintf(stderr, "\nERROR: Invalid flag for backbone_delta_score (expected 0 or 1): %d\n\t>%s", iTemp, sLine);
                Tempest::tempest_exit(EXIT_FAILURE); 
            }
        }

        else {
            fprintf(stderr, "\nWARNING: Unknown search parameter: '%s'\n", sParam);
            //Tempest::tempest_exit(EXIT_FAILURE);
        }
    }
    
    if (Tempest::params.numInternalPSMs < Tempest::params.numOutputPSMs)
        Tempest::params.numOutputPSMs = Tempest::params.numInternalPSMs;
        
}

void append_mod(char cAminoAcid, int modInd, double dMassDiff, char cSymbol) {
    Tempest::params.iNumMods += 1;
    if (0 == (Tempest::params.tMods = (mod_t*) realloc(Tempest::params.tMods, Tempest::params.iNumMods*sizeof(mod_t)))) {
        fprintf(stderr, "\nERROR: Could not allocate memory for %d mods.\n", Tempest::params.iNumMods);
        Tempest::tempest_exit(EXIT_FAILURE);
    }

    Tempest::params.tMods[Tempest::params.iNumMods-1].cAminoAcid = cAminoAcid;
    Tempest::params.tMods[Tempest::params.iNumMods-1].modInd = modInd;
    Tempest::params.tMods[Tempest::params.iNumMods-1].cSymbol = cSymbol;
    Tempest::params.tMods[Tempest::params.iNumMods-1].dMassDiff = dMassDiff;
}
