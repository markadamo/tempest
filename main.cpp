#include "tempest.hpp"

#ifdef MAC
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

std::vector<Device*> Tempest::devices;
cl_context clContextl;

int main(int argc, char** argv) {
    Tempest::initialize_globals();
    
    Tempest::parse_input(argc,argv);

    Tempest::parse_params();

    Tempest::parse_config();

    Tempest::setup_globals();
    
    for (int i=0; i<Tempest::config.iDevices.size(); i++)
        Tempest::devices.push_back(new Device(Tempest::config.iPlatform, i));
    
    Tempest::collect_msms_spectra();
    
    for (int i=0; i<Tempest::config.iDevices.size(); i++)
        Tempest::devices[i]->setup(0, Tempest::data.eScans.size());
    
    Tempest::search_fasta_database();
    
    for (int i=0; i<Tempest::config.iDevices.size(); i++)
        //wait for all kernels to finish before reading results
        Tempest::devices[i]->finish();
    
    Tempest::write_psms();
    
    Tempest::write_log(argc, argv);

    if (Tempest::config.profile) {
        for (int i=0; i<Tempest::devices.size(); i++)
            Tempest::devices[i]->printProfilingData();
    }
    
    Tempest::tempest_exit(EXIT_SUCCESS);
    return EXIT_SUCCESS;
}

/* 
 * Cleanup and Exit.
 */

extern void Tempest::tempest_exit(int EXIT_FLAG) {
    fflush(0);
    //cleanup_device();
    //cleanup_globals();
    exit(EXIT_FLAG);
}


/* End of File */
