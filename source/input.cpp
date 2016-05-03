#include "tempest.hpp"

void usage();
void version();
void set_short_argument(char, char*);
void set_long_argument(char*, char*);

/*
 * Parse command line input.
 */

extern void Tempest::parse_input(int argc, char **argv) {
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
            Tempest::tempest_exit(EXIT_FAILURE);
        }
    }
    
    // Method
    if (0 == Tempest::args.sParams) {
        Tempest::args.sParams = DEFAULT_PARAMS_FILE;
    }
    if (0 == (fp = (FILE *) fopen(Tempest::args.sParams, "r"))) {
        fprintf(stderr, "\nERROR\tUnable to open method file %s: %s\n", Tempest::args.sParams, strerror(errno));
        Tempest::tempest_exit(EXIT_FAILURE);
    }
    fclose(fp);

    // Config
    if (0 == Tempest::args.sConfig) {
        Tempest::args.sConfig = DEFAULT_CONFIG_FILE;
    }
    if (0 == (fp = (FILE *) fopen(Tempest::args.sConfig, "r"))) {
        fprintf(stderr, "\nERROR\tUnable to open config file %s: %s\n", Tempest::args.sConfig, strerror(errno));
        Tempest::tempest_exit(EXIT_FAILURE);
    }
    fclose(fp);

    // Out file
    if (!Tempest::args.sOut)
        Tempest::args.sOut = strdup_s(Tempest::args.sSpectra);
}

/*
 * Usage message - prints message and quits.
 */

extern void Tempest::usage() {
    printf("\nUSAGE tempest [options] -i <spectra_file> -f <FASTA_file>\n");
    printf("\nOPTIONS\n");
    printf("    -i <PATH>              spectra file (override method\n");
    printf("    -f <PATH>              fasta database (override method)\n");
    printf("    -m <PATH>              method file (default: tempest.method)\n");
    printf("    -c <PATH>              config file (default: tempest.config)\n");
    printf("    -o <PATH>              output file (default: <spectra path>.tsv)\n");
    printf("    -d <N:N>               OpenCL platform:device to use (override config)\n");
    printf("    -s <N:N>               Scan range (min:max MS2 scan IDs)\n");
    printf("    -l <N>                 print level (verbosity)\n");
    printf("    -h                     help\n");
    printf("    -v                     version\n");
    printf("\nMISC\n");
    printf("    --device-info          Print OpenCL device map.\n");
    //printf("    --split-output         Write results for each observed spectrum to a separate output file.\n");
    //printf("    --print-candidates     Print digested candidate peptides to stdout.\n");
    printf("\n");
}


// Version Info - prints version and quits.
void version() {
    printf("Tempest v%s %s\n", VERSION_STRING, LICENSE_STRING);
}

// Parse command line arguments
void set_short_argument(char arg, char *value) {
    char* token;
    switch(toupper(arg)) {
    case 'I': Tempest::args.sSpectra = Tempest::strdup_s(value); break;
    case 'F': Tempest::args.sFasta   = Tempest::strdup_s(value); break;
    case 'O': Tempest::args.sOut     = Tempest::strdup_s(value); break;
    case 'M': Tempest::args.sParams  = Tempest::strdup_s(value); break;
    case 'C': Tempest::args.sConfig  = Tempest::strdup_s(value); break;
    case 'D':
        token = strtok(value, ",.:|-/_");
        Tempest::config.iPlatform = atoi(token);
        while (token = strtok(NULL, ",.:|-/_"))
            Tempest::config.iDevices.push_back(atoi(token));
        Tempest::args.deviceOverride = 1;
        break;
    case 'S':
        token = strtok(value, ",.:|-/_");
        Tempest::config.minScan = atoi(token);
        token = strtok(NULL, ",.:|-/_");
        Tempest::config.maxScan = atoi(token);
        break;
    case 'L': Tempest::args.iPrintLevel = atoi(value); break;
    case 'V': version();
        Tempest::tempest_exit(EXIT_SUCCESS); break;
    case 'H': Tempest::usage();
        Tempest::tempest_exit(EXIT_SUCCESS); break;
    default:
        fprintf(stderr, "\nERROR\tUnknown argument (%c)\n", arg);
        Tempest::usage();
        Tempest::tempest_exit(EXIT_FAILURE);
    }
}

void set_long_argument(char* arg, char* value) {
    if (strcasecmp(arg, "device-info") == 0) {
        Tempest::print_device_info();
        Tempest::tempest_exit(EXIT_SUCCESS);
    }

    else if (strcasecmp(arg, "print-candidates") == 0) {
        Tempest::args.bPrintCandidates = 1;
    }

    else {
        fprintf(stderr, "\nERROR\tUnknown argument (%s)\n", arg);
        Tempest::usage();
        Tempest::tempest_exit(EXIT_FAILURE);
    }
}
