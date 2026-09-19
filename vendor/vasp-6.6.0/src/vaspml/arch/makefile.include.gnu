# --------
# Compiler
# --------
CXX             = mpic++

# --------------
# Compiler flags
# --------------
CXXFLAGS       := -O3
CXXFLAGS       += -g
CXXFLAGS       += -std=c++17 -pedantic-errors
CXXFLAGS       += -Wall -Wextra
# Extra flags required for testing and coverage.
#CXXFLAGS       += --coverage -fno-default-inline -fno-inline -fno-inline-small-functions -fno-elide-constructors

# -------------------
# Include directories
# -------------------
INCLUDE        := -I$(OPENBLAS_ROOT)/include

# -------
# Linking
# -------
LD              = $(CXX)
LDFLAGS        := -L$(OPENBLAS_ROOT)/lib -lopenblas

# ------------
# Preprocessor
# ------------
CPP_OPTIONS    :=
#CPP_OPTIONS    += -DVASPML_DEBUG_LEVEL=3

# --------------------
# Dependency generator
# --------------------
CPP_DEP         = $(CXX)

# --------------
# Archiving tool
# --------------
AR              = ar
ARFLAGS         = rcs

# -----------------------
# Additional make options
# -----------------------
#MAKE_OPTIONS   += --no-color
#MAKE_OPTIONS   += --no-logo

# ------------------------------
# Support for GRACE force fields
# (optional)
# ------------------------------
#CPP_OPTIONS     += -DVASPML_ENABLE_GRACE
#CPPFLOW_ROOT    ?= /path/to/your/cppflow/installation
## Auto-detect tensorflow root directory from python if available from current environment.
## Alternatively, enter /path/to/your/tensorflow/installation manually.
#TENSORFLOW_ROOT ?= $(shell dirname `python -c 'import tensorflow as tf; print(tf.__file__)' 2>/dev/null`)
#INCLUDE         += -I$(CPPFLOW_ROOT)/include -I$(TENSORFLOW_ROOT)/include
#LDFLAGS         += $(TENSORFLOW_ROOT)/libtensorflow_cc.so.2 $(TENSORFLOW_ROOT)/libtensorflow_framework.so.2
#LDFLAGS         += -Wl,-rpath,$(TENSORFLOW_ROOT)

# -----------------------------------
# Boost library
# (optional, only needed for testing)
# -----------------------------------
BOOST_ROOT     ?=
BOOST_INCLUDE   = -isystem${BOOST_ROOT}/include
BOOST_LIB       = -L${BOOST_ROOT}/lib -lboost_unit_test_framework
