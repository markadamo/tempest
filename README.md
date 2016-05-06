# tempest
OpenCL MS/MS database search

The included executable has been compiled for 64-bit linux, compatible with glibc version 2.12 or higher. 

+----------+
| Building |
+----------+
In order to build, Tempest requires the MSToolkit source in a directory called 'mstoolkit'. If you have git, run the following command from the tempest root directory:

    git clone https://github.com/mhoopmann/mstoolkit

Alternatively, MSToolkit can be downloaded as a .zip from https://github.com/mhoopmann/mstoolkit and extracted to the tempest directory.
Once MSToolkit has been cloned/extracted to 'mstoolkit', run 'make' from the tempest root directory (where the Makefile is).

 ----------

Kernel source code resides in kernels.cl, but is compiled into the executable as a byte array in kernels.h. Following modification of kernels.cl, navigate to the 'source' directory and run the following command before rebuilding with 'make':

    xxd -i kernels.cl > kernels.h

