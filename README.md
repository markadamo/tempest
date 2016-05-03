# tempest
OpenCL MS/MS database search

Kernel source code resides in kernels.cl, but is compiled into the executable as a byte array in kernels.h. Following modification of kernels.cl, navigate to the 'source' directory and run the following command before rebuilding with 'make':
       
       xxd -i kernels.cl > kernels.h

