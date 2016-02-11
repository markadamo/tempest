CC=g++
override CFLAGS=-g -O3 #-fsanitize=address
SOURCES=$(wildcard *.cpp)
#OBJECTS=$(SOURCES:.c=.o)
MSTOOLKIT=MSToolkit
INC_DIRS=MSToolkit/include
LIB_DIRS=$(MSTOOLKIT)
LIBS=-lmstoolkitlite #-lasan
EXECUTABLE=tempest

# Check for 32-bit vs 64-bit
PROC_TYPE = $(strip $(shell uname -m | grep 64))

# Check for Mac OS
OS = $(shell uname -s 2>/dev/null | tr [:lower:] [:upper:])
DARWIN = $(strip $(findstring DARWIN, $(OS)))

# MacOS System
ifneq ($(DARWIN),)
	CFLAGS += -DMAC
	LIBS += -framework OpenCL

	ifeq ($(PROC_TYPE),)
		CFLAGS+=-arch i386
	else
		CFLAGS+=-arch x86_64
	endif
else

# Linux OS
LIBS += -lm -lOpenCL
ifeq ($(PROC_TYPE),)
	CFLAGS+=-m32
else
	CFLAGS+=-m64
endif
# Check for Linux-AMD
ifdef AMDAPPSDKROOT
   INC_DIRS += $(AMDAPPSDKROOT)/include
	ifeq ($(PROC_TYPE),)
		LIB_DIRS += $(AMDAPPSDKROOT)/lib/x86
	else
		LIB_DIRS += $(AMDAPPSDKROOT)/lib/x86_64
	endif
else

# Check for Linux-Nvidia
ifdef CUDA
	INC_DIRS += $(CUDA)/include
else
	INC_DIRS += /usr/local/cuda/include
endif

endif
endif

#all: $(SOURCES) $(EXECUTABLE)

all:
	cd MSToolkit; make; cd ..
	$(CC) $(CFLAGS) $(SOURCES) $(INC_DIRS:%=-I%) $(LIB_DIRS:%=-L%) $(LIBS) -o $(EXECUTABLE)

print-%:
	@echo $* = $($*)

clean:
	rm *.o

