# AMReX
COMP = gnu
PRECISION = DOUBLE

# Profiling
PROFILE       = FALSE
TINY_PROFILE  = FALSE
COMM_PROFILE  = FALSE
TRACE_PROFILE = FALSE
MEM_PROFILE   = FALSE
USE_GPROF     = FALSE

# Performance
USE_MPI   = TRUE
USE_OMP   = FALSE
USE_CUDA  = FALSE
USE_HIP   = FALSE
USE_SYCL  = FALSE

# Debugging
DEBUG = FALSE

TEST = TRUE
USE_ASSERTION = TRUE

# Pytorch (set environment variable PYTORCH_ROOT = /path/to/torch)
ifeq ($(USE_CUDA),TRUE)
  PYTORCH_ROOT ?= ../libtorch_cuda
else
  PYTORCH_ROOT ?= ../libtorch_cpu
endif
TORCH_LIBPATH = $(PYTORCH_ROOT)/lib

ifeq ($(USE_CUDA),TRUE)
  TORCH_LIBS = -ltorch -ltorch_cpu -lc10 -lc10_cuda -lcuda 
else
  TORCH_LIBS = -ltorch -ltorch_cpu -lc10
endif

INCLUDE_LOCATIONS += $(PYTORCH_ROOT)/include \
                     $(PYTORCH_ROOT)/include/torch/csrc/api/include
LIBRARY_LOCATIONS += $(TORCH_LIBPATH)

DEFINES += -D_GLIBCXX_USE_CXX11_ABI=0
ifeq ($(USE_CUDA),TRUE)
  LDFLAGS += -Xlinker "--no-as-needed,-rpath $(TORCH_LIBPATH) $(TORCH_LIBS)"
else
  LDFLAGS += -Wl,--no-as-needed,-rpath=$(TORCH_LIBPATH) $(TORCH_LIBS)
endif

# GNU Make
Bpack := ./Make.package
Blocs := .

ERF_HOME := ../../..
ERF_PROBLEM_DIR = $(ERF_HOME)/Exec/MoistRegTests/Bubble
include $(ERF_HOME)/Exec/Make.ERF
