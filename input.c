
#include "tempest.h"

void usage();
void version();
void set_short_argument(char, char*);
void set_long_argument(char*, char*);

/*
 * Parse command line input.
 */

extern void parse_input(int argc, char **argv)
{
    int i;
    char short_argument;
    char long_argument[STRING_SIZE];
    char value[STRING_SIZE];
    FILE* fp;
    
    // Loop through command line arguments
    for (i=1; i<argc; i++) {

        // long format [--argument] or [--argument value]
        if (sscanf(argv[i], "--%s", long_argument) == 1) {
            if (i+1 == argc || argv[i+1][0] == '-') set_long_argument(long_argument, 0);
            else set_long_argument(long_argument, argv[++i]);
        }

        // short format [-X] or [-X value]
        else if (strlen(argv[i]) == 2 && sscanf(argv[i], "-%c", &short_argument) == 1) {
            if (i+1 == argc || argv[i+1][0] == '-') set_short_argument(short_argument, 0);
            else set_short_argument(short_argument, argv[++i]);
        }

        // short format [-Xvalue]
        else if (strlen(argv[i]) > 2 && sscanf(argv[i], "-%c%s", &short_argument, value) == 2) {
            set_short_argument(short_argument, value);
        }

        else {
            fprintf(stderr, "\nERROR\tCould not parse argument [%s]\n", argv[i]);
            tempest_exit(EXIT_FAILURE);
        }
    }
    
    // MzXML
    if (0 == args.sSpectra) {
        fprintf(stderr, "\nERROR\tNo spectra file given.\n");
        usage();
        tempest_exit(EXIT_FAILURE);
    }
    
    if ((fp = (FILE *) fopen(args.sSpectra, "r")) == NULL) {
        fprintf(stderr, "\nERROR\tUnable to open spectra file %s: %s\n", args.sSpectra, strerror(errno));
        tempest_exit(EXIT_FAILURE);
    }
    fclose(fp);
    
    // Fasta database
    if (0 == args.sFasta) {
        fprintf(stderr, "\nERROR\tNo fasta database given.\n");
        usage();
        tempest_exit(EXIT_FAILURE);
    }
    
    if (0 == (fp = (FILE *) fopen(args.sFasta, "r"))) {
        fprintf(stderr, "\nERROR\tUnable to open fasta database %s: %s\n", args.sFasta, strerror(errno));
        tempest_exit(EXIT_FAILURE);
    }
    fclose(fp);
    
    // Params
    if (0 == args.sParams) {
        args.sParams = strdup_s(DEFAULT_PARAMS_FILE);
    }
    
    if (0 == (fp = (FILE *) fopen(args.sParams, "r"))) {
        fprintf(stderr, "\nERROR\tUnable to open params file %s: %s\n", args.sParams, strerror(errno));
        tempest_exit(EXIT_FAILURE);
    }
    fclose(fp);

    // Out file
    if (!(args.sOut || sscanf(args.sSpectra, "%s.mzXML", args.sOut) == 1 || sscanf(args.sSpectra, "%s.mzxml", args.sOut) == 1)) {
        strcpy(args.sOut, args.sSpectra);
    }
}

/*
 * Usage message - prints message and quits.
 */

void usage()
{
    printf("\nUSAGE tempest <-m mzXML> <-f FASTA database>\n");
    printf("\nARGUMENTS\n");
    printf("    -i <PATH>              spectra file (mzXML)\n");
    printf("    -f <PATH>              fasta database\n");
    printf("    -p <PATH>              params file (default: tempest.params)\n");
    printf("    -o <PATH>              output file (default: <mzxml path>.csv)\n");
    printf("    -d <N>                 OpenCL platform:device to use (default: 0:0)\n");
    printf("    -l <N>                 print level (verbosity)\n");
    printf("    -h                     help\n");
    printf("    -v                     version\n");
    printf("\nMISC\n");
    printf("    --split-output         Write results for each observed spectrum to a separate output file.\n");
    printf("    --print-candidates     Print digested candidate peptides to stdout.\n");
    printf("    --no-digest            Process MS/MS, but do not digest fasta database.\n");
    printf("    --no-match             Process MS/MS and digest fasta, but do not assign candidates.\n");
    printf("    --no-score             Process MS/MS and assign candidates, but do not produce similarity scores.\n");
    printf("\nGPU CONFIGURATION\n");
    printf("    --force-no-gpu         Do not use CUDA-enabled GPUs, even if available.\n");
    printf("    --force-shared         Use shared GPU memory, even if ...\n");
    printf("    --force-no-prebuild    Do not prebuild MS/MS spectra, even if space allows.\n");
    printf("    --buffer-size <N>      Num candidates to collect before scoring (default %d).\n", DEFAULT_CANDIDATE_BUFFER_SIZE);
    printf("    --block-dim <N>        Num kernel blocks when scoring candidates (default %d).\n", DEFAULT_BLOCKDIM_SCORE);
    printf("\n");
}

/*
 * Version Info - prints version and quits.
 */

void version()
{
    printf("Tempest v%s %s\n", VERSION_STRING, LICENSE_STRING);
}

/*
 * The Major Workhorse of common - set variables from command line input
 */
void set_short_argument(char arg, char *value)
{
    char* token;
    switch(toupper(arg)) {
    case 'I': args.sSpectra = strdup_s(value);  break;
    case 'F': args.sFasta = strdup_s(value);    break;
    case 'O': args.sOut = strdup_s(value);      break;
    case 'P': args.sParams = strdup_s(value);   break;
    case 'D':
        token = strtok(value, ",.;:|-/_");
        config.iPlatform = atoi(token);
        token = strtok(NULL, ",.;:|-/_");
        config.iDevice = atoi(token);
        break;
    case 'S':
        token = strtok(value, ",.;:|-/_");
        config.minScan = atoi(token);
        token = strtok(NULL, ",.;:|-/_");
        config.maxScan = atoi(token);
        break;
    case 'L': args.iPrintLevel = atoi(value);   break;
    case 'V': version();
        tempest_exit(EXIT_SUCCESS);     break;
    case 'H': usage();
        tempest_exit(EXIT_SUCCESS);     break;
    default:
        fprintf(stderr, "\nERROR\tUnknown argument (%c)\n", arg);
        usage();
        tempest_exit(EXIT_FAILURE);
    }
}

void set_long_argument(char* arg, char* value)
{
    if (strcasecmp(arg, "device-info") == 0) {
        print_device_info();
        tempest_exit(EXIT_SUCCESS);
    }
    
    if (strcasecmp(arg, "split-output") == 0) {
        args.bSplitOutput = 1;
    }

    else if (strcasecmp(arg, "print-candidates") == 0) {
        args.bPrintCandidates = 1;
    }

    else if (strcasecmp(arg, "no-digest") == 0) {
        args.bNoDigest = 1;
    }

    else if (strcasecmp(arg, "no-match") == 0) {
        args.bNoMatch = 1;
    }

    else if (strcasecmp(arg, "no-score") == 0) {
        args.bNoScore = 1;
    }

    else if (strcasecmp(arg, "force-no-gpu") == 0) {
        config.bForceNoGPU = 1;
    }

    else if (strcasecmp(arg, "force-shared") == 0) {
        config.bForceShared = 1;
    }

    else if (strcasecmp(arg, "force-no-prebuild") == 0) {
        config.bForceNoPrebuild = 1;
    }

    else if (strcasecmp(arg, "buffer-size") == 0) {
        if (value) {
            config.iCandidateBufferSize = atoi(value);
        }
        else {
            fprintf(stderr, "\nERROR\tArgument --buffer-size requires a value\n");
            usage();
            tempest_exit(EXIT_FAILURE);
        }
    }

    else if (strcasecmp(arg, "block-dim") == 0) {
        if (value) {
            config.iScoreBlockDim = atoi(value);
        }
        else {
            fprintf(stderr, "\nERROR\tArgument --block-dim requires a value\n");
            usage();
            tempest_exit(EXIT_FAILURE);
        }

        if (config.iScoreBlockDim > MAX_BLOCKDIM) {
            fprintf(stderr, "\nERROR\tInvalid block dimension: %d (max block dimension is %d)\n", config.iScoreBlockDim, MAX_BLOCKDIM);
            usage();
            tempest_exit(EXIT_FAILURE);
        }
    }

    else {
        fprintf(stderr, "\nERROR\tUnknown argument (%s)\n", arg);
        usage();
        tempest_exit(EXIT_FAILURE);
    }
}
