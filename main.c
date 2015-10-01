/*
  Tempest: CUDA-enabled Database Search for Tandem Mass Spectrometry
*/

#include "tempest.h"

#ifdef MAC
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

int main(int argc, char **argv) 
{
    initialize_globals();
    parse_input(argc,argv);
    parse_params();
    setup_globals();
    initialize_device();
    collect_msms_spectra();
    setup_device();
    setup_constant_memory();
    search_fasta_database();
    write_psms();
    write_log();

    if (PROFILE) {
        printf("Total memset time:       %fs\n", (float)totalMemsetTime/1000000000);
        printf("Total build time:        %fs\n", (float)totalBuildTime/1000000000);
        printf("Total transform time:    %fs\n", (float)totalTransformTime/1000000000);
        printf("Total send time:         %fs\n", (float)totalSendTime/1000000000);
        printf("Total score time:        %fs\n", (float)totalScoreTime/1000000000);
        printf("Total reduce time:       %fs\n", (float)totalReduceTime/1000000000);
        printf("Scoring kernel launches: %ld\n", scoreKernelLaunches);
    }
    
    tempest_exit(EXIT_SUCCESS);
    return EXIT_SUCCESS;
}

/* 
 * Cleanup and Exit.
 */

extern void tempest_exit(int EXIT_FLAG)
{
    fflush(0);
    //cleanup_device();
    //cleanup_globals();
    exit(EXIT_FLAG);
}


/* End of File */
