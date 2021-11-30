# Makefile for Song.
# Created by Guido Walter Pettinari on 2010.

# =========================================================================
# =                         Compilers and flags                           =
# =========================================================================

# Your C compiler and library tool
# To compile on a vanilla Mac, you might need to disable openmp
# support, see "Parallelisation flags" below.
CC = gcc
AR = ar rv

# Optimization flags
OPTFLAG = -O3

# Debug flags
CFLAGS += -g
CFLAGS += -DDEBUG

# Allow C99 syntax
CFLAGS += -std=c99

# Uncomment to inhibit all warning messages
CFLAGS   += -w

# Uncomment to enable most warning messages
# CFLAGS += -Wall -Wno-unused -Wno-logical-not-parentheses

# Parallelisation flags
# Uncomment both rows if you not need parallel support.
# If you are compiling on a Mac, in order to enable parallel support
# you'll need to either manually install openmp for Clang (see
# https://stackoverflow.com/a/39843038/2972183) or to use a standard
# implementation of GCC (ex. GNU GCC using Homebrew).

# Header files and libraries
INCLUDES = -I../include -I../$(CLASS_DIR)/include
LDFLAGS = -lm

CFLAGS += -fopenmp
LDFLAGS += -fopenmp


# =========================================================================
# =                        Directories locations                          =
# =========================================================================

# Build directory
MDIR := .
WRKDIR = $(MDIR)/build

# CLASS subfolder
CLASS_DIR = $(MDIR)/class.git

# Source files, including those in CLASS subfolder
vpath %.c source:tools:main:test:$(CLASS_DIR)/source:$(CLASS_DIR)/tools:$(CLASS_DIR)/main:$(CLASS_DIR)/test
vpath %.o build:$(CLASS_DIR)/build
vpath .base build

.base:
	if ! [ -a $(WRKDIR) ]; then mkdir $(WRKDIR) ; mkdir $(WRKDIR)/lib; fi;
	touch build/.base

# Make sure to compile CLASS with support for SONG
CFLAGS += -DWITH_SONG1 # support for SONG features requiring a first-order computation
CFLAGS += -DWITH_SONG2 # support for SONG features requiring a second-order computation

# Pass to the code CLASS location
CFLAGS += -D__CLASSDIR__='"$(CLASS_DIR)"'


# =========================================================================
# =                          Compilation rules                            =
# =========================================================================

%.o:  %.c .base
	cd $(WRKDIR); $(CC) $(OPTFLAG) $(CFLAGS) $(INCLUDES) -c ../$< -o $*.o


# ==========================================================================
# =                             Object files                               =
# ==========================================================================

# Source files also present in CLASS
TOOLS = growTable.o dei_rkck.o sparse.o evolver_rkck.o  evolver_ndf15.o \
arrays.o parser.o quadrature.o hyperspherical.o common.o
SOURCE_CLASS = input.o background.o thermodynamics.o perturbations.o primordial.o \
nonlinear.o transfer.o spectra.o lensing.o bessel.o bispectra.o fisher.o output.o
INPUT = input.o
BACKGROUND = background.o
THERMO = thermodynamics.o
PERTURBATIONS = perturbations.o
TRANSFER = transfer.o
PRIMORDIAL = primordial.o
SPECTRA = spectra.o
NONLINEAR = trg.o nonlinear.o
LENSING = lensing.o
OUTPUT = output.o

# Source files exclusive of SONG
SONG_TOOLS = $(TOOLS) utility.o song_tools.o slatec_3j_C.o mesh_interpolation.o binary.o
INPUT2 = input2.o
PERTURBATIONS2 = perturbations2.o
BESSEL = bessel.o
BESSEL2 = bessel2.o
TRANSFER2 = transfer2.o
BISPECTRA = bispectra.o
BISPECTRA2 = bispectra2.o
SPECTRA2 = spectra2.o
FISHER = fisher.o


# ==========================================================================
# =                                Targets                                 =
# ==========================================================================

# Rule to be executed when make is called with no parameters
default: class song print_params print_sources1 print_sources2 print_transfers1 print_transfers2 print_bispectra

# CLASS executables
.PHONY: libclass.a class test_background test_thermodynamics test_perturbations test_transfer classy tar
libclass.a class test_background test_thermodynamics test_perturbations test_transfer classy tar: 
	cd $(CLASS_DIR);\
	export WITH_SONG1=1;\
	export WITH_SONG2=1;\
	make $@; mv $@ ..

# SONG executables
SONG = song.o
PRINT_PARAMS = print_params.o
PRINT_SOURCES1 = print_sources1.o
PRINT_SOURCES2 = print_sources2.o
PRINT_TRANSFERS1 = print_transfers1.o
PRINT_TRANSFERS2 = print_transfers2.o
PRINT_BISPECTRA = print_bispectra.o
PRINT_BACKGROUND = print_background.o
PRINT_THERMO = print_thermo.o
PRINT_COUPLINGS = print_couplings.o
PRINT_TAU_SONG = print_tau_song.o
PRINT_K_SONG = print_k_song.o

song: $(SONG_TOOLS) $(SOURCE_CLASS) $(INPUT2) $(PERTURBATIONS2)\
	$(BESSEL2) $(TRANSFER2) $(BISPECTRA2) $(SPECTRA2) $(OUTPUT) $(SONG)
	$(CC) $(LDFLAGS) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_params: $(SONG_TOOLS) $(SOURCE_CLASS) $(INPUT2) $(PRINT_PARAMS)
	$(CC) $(LDFLAGS) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_sources1: $(SONG_TOOLS) $(SOURCE_CLASS) $(INPUT2) $(PERTURBATIONS2) \
	$(PRINT_SOURCES1)
	$(CC) $(LDFLAGS) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_sources2: $(SONG_TOOLS) $(SOURCE_CLASS) $(INPUT2) $(PERTURBATIONS2) \
	$(PRINT_SOURCES2)
	$(CC) $(LDFLAGS) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_transfers1: $(SONG_TOOLS) $(SOURCE_CLASS) $(INPUT2) $(PERTURBATIONS2) $(PRINT_TRANSFERS1)
	$(CC) $(LDFLAGS) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_transfers2: $(SONG_TOOLS) $(SOURCE_CLASS) $(INPUT2) $(PERTURBATIONS2)\
	$(BESSEL2) $(TRANSFER2) $(PRINT_TRANSFERS2)
	$(CC) $(LDFLAGS) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_bispectra: $(SONG_TOOLS) $(SOURCE_CLASS) $(INPUT2) $(PERTURBATIONS2)\
	$(BESSEL2) $(TRANSFER2) $(BISPECTRA2) $(SPECTRA2) $(OUTPUT) $(PRINT_BISPECTRA)
	$(CC) $(LDFLAGS) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_background: $(SONG_TOOLS) $(SOURCE_CLASS) $(INPUT2) $(PRINT_BACKGROUND)
	$(CC) $(LDFLAGS) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_thermo: $(SONG_TOOLS) $(SOURCE_CLASS) $(INPUT2) $(PRINT_THERMO)
	$(CC) $(LDFLAGS) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_couplings: $(SONG_TOOLS) $(SOURCE_CLASS) $(PRINT_COUPLINGS)
	$(CC) $(LDFLAGS) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_tau_song: $(SONG_TOOLS) $(SOURCE_CLASS) $(INPUT2) $(PERTURBATIONS2) $(PRINT_TAU_SONG)
	$(CC) $(LDFLAGS) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_k_song: $(SONG_TOOLS) $(SOURCE_CLASS) $(INPUT2) $(PERTURBATIONS2) $(PRINT_K_SONG)
	$(CC) $(LDFLAGS) -o  $@ $(addprefix build/,$(notdir $^)) -lm

clean: .base
	rm -rf $(WRKDIR);
	cd $(CLASS_DIR); $(MAKE) clean;


# THE EXECUTABLES BELOW THIS POINT ARE BROKEN

PRINT_K = print_k.o
PRINT_INITIAL_CONDITIONS_SONG = print_initial_conditions_song.o
PRINT_KERNEL = print_kernel.o

print_k: $(SONG_TOOLS) $(SOURCE_CLASS) $(INPUT2) $(PRINT_K)
	$(CC) $(LDFLAGS) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_k_song: $(SONG_TOOLS) $(SOURCE_CLASS) $(PRINT_K_SONG)
	$(CC) $(LDFLAGS) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_initial_conditions_song: $(SONG_TOOLS) $(SOURCE_CLASS) $(PERTURBATIONS2) \
	$(PRINT_INITIAL_CONDITIONS_SONG)
	$(CC) $(LDFLAGS) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_kernel: $(SONG_TOOLS) $(SOURCE_CLASS) $(PERTURBATIONS2) $(PRINT_KERNEL)
	$(CC) $(LDFLAGS) -o  $@ $(addprefix build/,$(notdir $^)) -lm


