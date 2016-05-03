#include "tempest.hpp"

extern void Tempest::parse_config() {
    FILE* fp;
    char sLine[STRING_SIZE];
    char sParam[STRING_SIZE];
    char sTemp[STRING_SIZE];

    if((fp = (FILE*) fopen(Tempest::args.sConfig, "r")) == NULL) {
        fprintf( stderr, "ERROR\tUnable to open config file %s: %s\n", Tempest::args.sParams, strerror(errno) );
        Tempest::tempest_exit(EXIT_FAILURE);
    }

    while (fgets(sLine, sizeof(sLine), fp)) {
        //change > to null terminator to stop parsing line at comment symbol
        for (int i=0; i<STRING_SIZE; i++) {
            if (sLine[i] == '>')
                sLine[i] = '\0';
            if (sLine[i] == '\0')
                break;
        }
        if (sscanf(sLine, "%s", sParam) == -1) continue;
        

        if (strcmp(sParam, "platform") == 0) {
            if (Tempest::args.deviceOverride)
                continue;
            if (Tempest::config.iPlatform != -1) {
                fprintf(stderr, "\nERROR: Multiple platform definitions:\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }
            if (sscanf(sLine, "platform %d \n", &Tempest::config.iPlatform) != 1) {
                fprintf(stderr, "\nERROR: Could not parse platform definition:\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }
        }

        else if (strcmp(sParam, "device") == 0) {
            if (Tempest::args.deviceOverride)
                continue;
            int deviceID;
            if (sscanf(sLine, "device %d \n", &deviceID) != 1) {
                fprintf(stderr, "\nERROR: Could not parse device definition:\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }
            Tempest::config.iDevices.push_back(deviceID);
        }

        else if (strcmp(sParam, "max_host_mem") == 0) {
            if (Tempest::args.hostMemOverride)
                continue;
            float maxMem;
            if (sscanf(sLine, "max_host_mem %f %s\n", &maxMem, sTemp) != 2) {
                fprintf(stderr, "\nERROR: Could not parse max_host_mem definition:\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }
            if (maxMem < 0) {
                fprintf(stderr, "\nERROR: Negative value in parse max_host_mem definition:\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }
            char unit = tolower(sTemp[0]);
            if (unit!='b' && unit!='k' && unit!='m' && unit!='g') {
                fprintf(stderr, "\nERROR: Invalid units in parse max_host_mem definition:\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }
            int factor;
            if      (unit == 'k') factor = KB;
            else if (unit == 'm') factor = MB;
            else if (unit == 'g') factor = GB;
            else                  factor = 1;
            Tempest::config.maxHostMem = (size_t)(maxMem*factor);
        }
        
        else if (strcmp(sParam, "max_device_mem") == 0) {
            if (Tempest::args.deviceMemOverride)
                continue;
            float maxMem;
            if (sscanf(sLine, "max_device_mem %f %s\n", &maxMem, sTemp) != 2) {
                fprintf(stderr, "\nERROR: Could not parse max_device_mem definition\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }
            if (maxMem < 0) {
                fprintf(stderr, "\nERROR: Invalid value in parse max_device_mem definition\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }
            char unit = tolower(sTemp[0]);
            if (unit!='b' && unit!='k' && unit!='m' && unit!='g') {
                fprintf(stderr, "\nERROR: Invalid units in parse max_device_mem definition\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }
            int factor;
            if      (unit == 'k') factor = KB;
            else if (unit == 'm') factor = MB;
            else if (unit == 'g') factor = GB;
            else                  factor = 1;
            Tempest::config.maxDeviceMem = (size_t)(maxMem*factor);
        }

        else if (strcmp(sParam, "min_work_size") == 0) {
            if (sscanf(sLine, "min_work_size %d \n", &Tempest::config.minWorkSize) != 1) {
                fprintf(stderr, "\nERROR: Could not parse min_work_size:\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }
        }

        else if (strcmp(sParam, "parallel_reduce") == 0) {
            int temp;
            if (sscanf(sLine, "parallel_reduce %d\n", &temp) != 1) {
                fprintf(stderr, "\nERROR: Could not parse parallel_reduce flag:\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }
            if (temp != 0 && temp != 1) {
                fprintf(stderr, "\nERROR: Invalid value in parallel_reduce flag:\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }
            else
                Tempest::config.parallelReduce = temp;
        }

        else if (strcmp(sParam, "profile") == 0) {
            int temp;
            if (sscanf(sLine, "profile %d\n", &temp) != 1) {
                fprintf(stderr, "\nERROR: Could not parse profile flag:\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }
            if (temp != 0 && temp != 1) {
                fprintf(stderr, "\nERROR: Invalid value in profile flag:\n\t>%s", sLine);
                Tempest::tempest_exit(EXIT_FAILURE);
            }
            else
                Tempest::config.profile = temp;
        }
            

	else
	  fprintf(stderr, "\nWARNING: unknown config parameter: '%s'\n", sParam);
    }
}
