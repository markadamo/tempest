
#include "tempest.h"

#ifdef MAC
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif


#ifdef _WIN32
#define PADDING 3
#else
#define PADDING 2
#endif

/*
 * Checks for OpenCL errors. Prints message and error string and exits on error.
 */

extern void check_cl_error(char* file, int line, int err, char* sMessage) {
    if (err < 0) {
        fprintf(stderr, "\n");
        fprintf(stderr, "OpenCL Error: %s@%d\t%s\n", file, line, sMessage);
        fprintf(stderr, "          \t%s (%d)\n\n", get_error_string(err), err);
        tempest_exit(EXIT_FAILURE);
    }
}

const char * get_error_string(cl_int err){
    switch(err){
    case   0: return "CL_SUCCESS";
    case  -1: return "CL_DEVICE_NOT_FOUND";
    case  -2: return "CL_DEVICE_NOT_AVAILABLE";
    case  -3: return "CL_COMPILER_NOT_AVAILABLE";
    case  -4: return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
    case  -5: return "CL_OUT_OF_RESOURCES";
    case  -6: return "CL_OUT_OF_HOST_MEMORY";
    case  -7: return "CL_PROFILING_INFO_NOT_AVAILABLE";
    case  -8: return "CL_MEM_COPY_OVERLAP";
    case  -9: return "CL_IMAGE_FORMAT_MISMATCH";
    case -10: return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
    case -11: return "CL_BUILD_PROGRAM_FAILURE";
    case -12: return "CL_MAP_FAILURE";

    case -30: return "CL_INVALID_VALUE";
    case -31: return "CL_INVALID_DEVICE_TYPE";
    case -32: return "CL_INVALID_PLATFORM";
    case -33: return "CL_INVALID_DEVICE";
    case -34: return "CL_INVALID_CONTEXT";
    case -35: return "CL_INVALID_QUEUE_PROPERTIES";
    case -36: return "CL_INVALID_COMMAND_QUEUE";
    case -37: return "CL_INVALID_HOST_PTR";
    case -38: return "CL_INVALID_MEM_OBJECT";
    case -39: return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
    case -40: return "CL_INVALID_IMAGE_SIZE";
    case -41: return "CL_INVALID_SAMPLER";
    case -42: return "CL_INVALID_BINARY";
    case -43: return "CL_INVALID_BUILD_OPTIONS";
    case -44: return "CL_INVALID_PROGRAM";
    case -45: return "CL_INVALID_PROGRAM_EXECUTABLE";
    case -46: return "CL_INVALID_KERNEL_NAME";
    case -47: return "CL_INVALID_KERNEL_DEFINITION";
    case -48: return "CL_INVALID_KERNEL";
    case -49: return "CL_INVALID_ARG_INDEX";
    case -50: return "CL_INVALID_ARG_VALUE";
    case -51: return "CL_INVALID_ARG_SIZE";
    case -52: return "CL_INVALID_KERNEL_ARGS";
    case -53: return "CL_INVALID_WORK_DIMENSION";
    case -54: return "CL_INVALID_WORK_GROUP_SIZE";
    case -55: return "CL_INVALID_WORK_ITEM_SIZE";
    case -56: return "CL_INVALID_GLOBAL_OFFSET";
    case -57: return "CL_INVALID_EVENT_WAIT_LIST";
    case -58: return "CL_INVALID_EVENT";
    case -59: return "CL_INVALID_OPERATION";
    case -60: return "CL_INVALID_GL_OBJECT";
    case -61: return "CL_INVALID_BUFFER_SIZE";
    case -62: return "CL_INVALID_MIP_LEVEL";
    case -63: return "CL_INVALID_GLOBAL_WORK_SIZE";
    default: return "Unknown OpenCL error";
    }
}

/*
 * Print information about OpenCL platforms/devices available on the system.
 */

extern void print_device_info() {
    cl_platform_id platforms[100];
    cl_uint platforms_n = 0;
    int err;

    err = clGetPlatformIDs(100, platforms, &platforms_n);
    check_cl_error(__FILE__, __LINE__, err, "Unable to enumerate OpenCL platforms");

    //printf("=== %d OpenCL platform(s) found: ===\n", platforms_n);   
    for (int i=0; i<platforms_n; i++)
        {
            char buffer[10240];
            printf("  ╔═\e[1mPlatform %d\e[22m\n", i);
            err |= clGetPlatformInfo(platforms[i], CL_PLATFORM_VERSION, 10240, buffer, NULL);
            printf("  ╟ VERSION = %s\n", buffer);
            err |= clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, 10240, buffer, NULL);
            printf("  ╟ NAME    = %s\n", buffer);
            err |= clGetPlatformInfo(platforms[i], CL_PLATFORM_VENDOR, 10240, buffer, NULL);
            printf("  ╟ VENDOR  = %s\n", buffer);
            check_cl_error(__FILE__, __LINE__, err, "Unable to get OpenCL platform attributes");
            //err = clGetPlatformInfo(platforms[i], CL_PLATFORM_EXTENSIONS, 10240, buffer, NULL);
            //printf("  EXTENSIONS = %s\n", buffer);

            cl_device_id devices[100];
            cl_uint devices_n = 0;
            err = clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, 100, devices, &devices_n);
            check_cl_error(__FILE__, __LINE__, err, "Unable to enumerate OpenCL devices");

            //printf("      === %d OpenCL device(s) found on platform %d: ===\n", i);
            for (int j=0; j<devices_n; j++)
                {
                    char buffer[10240];
                    cl_uint buf_uint;
                    cl_ulong buf_ulong;
                    printf("  ║ ┌─\e[1mDevice %d (%d:%d)\e[22m\n", j, i, j);
                    err = clGetDeviceInfo(devices[j], CL_DEVICE_NAME, sizeof(buffer), buffer, NULL);
                    printf("  ║ ├ NAME                = %s\n", buffer);
                    err |= clGetDeviceInfo(devices[j], CL_DEVICE_VENDOR, sizeof(buffer), buffer, NULL);
                    printf("  ║ ├ VENDOR              = %s\n", buffer);
                    err |= clGetDeviceInfo(devices[j], CL_DEVICE_VERSION, sizeof(buffer), buffer, NULL);
                    printf("  ║ ├ DEVICE_VERSION      = %s\n", buffer);
                    err |= clGetDeviceInfo(devices[j], CL_DRIVER_VERSION, sizeof(buffer), buffer, NULL);
                    printf("  ║ ├ DRIVER_VERSION      = %s\n", buffer);
                    err |= clGetDeviceInfo(devices[j], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(buf_uint), &(buf_uint), NULL);
                    printf("  ║ ├ MAX_COMPUTE_UNITS   = %u\n", (unsigned int)buf_uint);
                    err |= clGetDeviceInfo(devices[j], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(buf_ulong), &(buf_ulong), NULL);
                    printf("  ║ ├ MAX_WORK_GROUP_SIZE = %llu\n", (unsigned long long)buf_ulong);
                    err |= clGetDeviceInfo(devices[j], CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(buf_uint), &buf_uint, NULL);
                    printf("  ║ ├ MAX_CLOCK_FREQUENCY = %u\n", (unsigned int)buf_uint);
                    err |= clGetDeviceInfo(devices[j], CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(buf_ulong), &buf_ulong, NULL);
                    printf("  ║ ├ GLOBAL_MEM_SIZE     = %llu\n", (unsigned long long)buf_ulong);
                    err |= clGetDeviceInfo(devices[j], CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(buf_ulong), &buf_ulong, NULL);
                    printf("  ║ ├ MAX_MEM_ALLOC_SIZE  = %llu\n", (unsigned long long)buf_ulong);
                    printf("  ║ └─■\n");
                    check_cl_error(__FILE__, __LINE__, err, "Unable to get OpenCL device attributes");
                }
            printf("  ╚═■\n");
        }
    printf("\n");

    if (platforms_n == 0) {
        printf("  No OpenCL platforms found!\n");
        printf("  Make sure OpenCL drivers are installed for your hardware.\n");
    }
}

/*
 *  Duplicate a string - wraps malloc and strcpy with error checking.
 */

extern char* strdup_s(char* s) {
    int n = strlen(s);
    char* d = (char*) malloc(n + PADDING);
    if (0 == s) {
        fprintf(stderr, "\nERROR! Unable to allocate memory.\n");
        tempest_exit(EXIT_FAILURE);
    }
    strcpy(d,s);
    d[n+1] = '\0';
    return d;
}

/*
 * Calculates combinations (binomial coefficient)
 */

extern int n_choose_k(int n, int k)
{
    int i, c;
    
    if (k > n-k) k = n-k;
    for (i=0, c=1; i<k; i++) {
        c = c * (n-i);
        c = c / (i+1);
    }

    return c;
}

#ifdef _WIN32

/*
 * roundf - based on http://forums.belution.com/en/cpp/000/050/13.shtml 
 */

extern float roundf(float value)
{
    if (value < 0) return (float) -(floor(-value + 0.5));
    else return  (float) floor(value + 0.5);
}

#endif

/*
 * Round up to nearest multiple (useful for global/local work dimensions)
 */

extern long mround(long num, int multiple) {
    return num < multiple ? multiple : ((num + multiple - 1) / multiple) * multiple;
}

/*
 * convert to string binary
 */
    
extern const char *byte_to_binary(unsigned int x)
{
    static char b[17];
    unsigned int z;
    b[0] = '\0';
    
    for (z = 32768; z > 0; z >>= 1) {
        strcat(b, ((x & z) == z) ? "1" : "0");
    }

    return b;
}

/*
 * Count set bits ("Hamming Weight")
 */

/*extern int count_set_bits(unsigned int v) {
  v = v - ((v >> 1) & 0x55555555);                         // reuse input as temporary
  v = (v & 0x33333333) + ((v >> 2) & 0x33333333);          // temp
  return ((v + ((v >> 4) & 0xF0F0F0F)) * 0x1010101) >> 24; // count
  }*/

extern int count_set_bits(unsigned int v) {
    unsigned int c;
    for (c=0; v; c++) v &= v-1;
    return c;
}

/*
 *
 */

extern char* modpeptide(const char* sPeptide, unsigned int uiModPattern) {
    static int i;
    static char *c;
    static unsigned int uiMod;
    static char *sModPeptide=0;

    if (sModPeptide == 0) sModPeptide = (char*) malloc((MAX_PEPTIDE_LENGTH + params.iModificationsMax) * sizeof(char));

    uiMod = 1;
    
    // advance uiMod bit for variable nterm
    if (params.cVariableNtermSymbol) uiMod <<= 1;

    for (c=sModPeptide, i=0; *c; c++) {
        // aa
        sModPeptide[i++] = *c;
        
        // aa symbol
        if (cModSites[aa2i(*c)]) {
            if (uiMod & uiModPattern) sModPeptide[i++] = cModSites[aa2i(*c)];
            uiMod <<= 1;
        }

        // nterm symbol
        if (c==sModPeptide && params.cVariableNtermSymbol && (uiModPattern & (unsigned int) 1)) {
            sModPeptide[i++] = params.cVariableNtermSymbol;
        }
    }

    // cterm symbol
    if (params.cVariableCtermSymbol && (uiMod & uiModPattern)) {
        sModPeptide[i++] = params.cVariableCtermSymbol;
    }

    sModPeptide[i] = '\0';

    return sModPeptide;
}

extern char* modnpeptide(const char* sPeptide, int iLen, unsigned int uiModPattern) {
    static int i, j;
    static unsigned int uiMod;
    static char *sModPeptide=0;

    if (sModPeptide == 0) sModPeptide = (char*) malloc((MAX_PEPTIDE_LENGTH + params.iModificationsMax + 1) * sizeof(char));

    uiMod=1;

    // advance uiMod bit for variable nterm
    if (params.cVariableNtermSymbol) uiMod <<= 1;

    for (i=0, j=0; i<iLen; i++) {
        // aa
        sModPeptide[j++] = sPeptide[i];
        
        // aa symbol
        if (cModSites[aa2i(sPeptide[i])]) {
            if (uiMod & uiModPattern) sModPeptide[j++] = cModSites[aa2i(sPeptide[i])];
            uiMod <<= 1;
        }

        // nterm symbol
        if (i==0 && params.cVariableNtermSymbol && (uiModPattern & (unsigned int) 1)) {
            sModPeptide[j++] = params.cVariableNtermSymbol;
        }
    }
    
    // cterm symbol
    if (params.cVariableCtermSymbol && (uiMod & uiModPattern)) {
        sModPeptide[j++] = params.cVariableCtermSymbol;
    }

    sModPeptide[j] = '\0';

    return sModPeptide;
}

extern const char* get_setup_mode(enum SETUP_MODE e) {
    switch (e) {
    case NO_PREBUILD:  return "No Prebuild (copy only)";
    case SEQ_PREBUILD: return "Sequential Prebuild";
    case PAR_PREBUILD: return "Parallel Prebuild";
    }
    return "Unknown";
}

extern const char* get_score_mode(enum SCORE_MODE e) {
    switch (e) {
    case SCORE:                    return "Prebuilt Score";
    case BUILD_SCORE:              return "Build, Score";
    case BUILD_TRANS_SCORE:        return "Build, Transform, Score";
    case SCORE_SHARED:             return "Prebuilt Score - Shared";
    case BUILD_SCORE_SHARED:       return "Build+Score - Shared";
    case BUILD_TRANS_SCORE_SHARED: return "Build+Transform+Score - Shared";
    }
    return "Unknown";
}

extern unsigned long hash(char *s) {
    unsigned long hash = 5381;
    int c;

    while ((c = *s)) {
        hash = ((hash << 5) + hash) + c; // hash * 33 + c
    }

    return hash;
}

extern unsigned long hashn(char *s, int n) {
    unsigned long hash = 5381;
    int i;

    for (i=0; i<n; i++) {
        hash = ((hash << 5) + hash) + s[i]; // hash * 33 + c
    }

    return hash;
}
