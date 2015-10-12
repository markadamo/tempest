#include "tempest.hpp"

#ifdef MAC
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

std::vector<Device*> Tempest::devices;
cl_context clContextl;

int main(int argc, char **argv) 
{
    initialize_globals();
    parse_input(argc,argv);
    parse_params();
    setup_globals();
    collect_msms_spectra();
    for (int i=0; i<config.iDevices.size(); i++)
        Tempest::devices.push_back(new Device(config.iPlatform, i, KERNEL_FILE, 0, eScans.size()));
    search_fasta_database();
    for (int i=0; i<config.iDevices.size(); i++)
        //wait for all kernels to finish before reading results
        Tempest::devices[i]->finish();
    write_psms();
    write_log();

    if (PROFILE) {
        for (int i=0; i<Tempest::devices.size(); i++)
            Tempest::devices[i]->printProfilingData();
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
