# Makefile for Song.
# Created by Guido Walter Pettinari on 2010.

# =========================================================================
# =                         Compilers and flags                           =
# =========================================================================

# Your C compiler and library tool
CC			 = gcc
AR       = ar rv

# Optimization flags
OPTFLAG  = -O4

# Compilation flags
CFLAGS   += -g
# CFLAGS   += -w
# CFLAGS   += -Wall -Wno-unused
CFLAGS   += -std=c99
CFLAGS   += -fopenmp
CFLAGS   += -DDEBUG

# Header files and libraries
INCLUDES 	 = -I../include -I../$(CLASS_DIR)/include
LIBRARIES  = -fopenmp -lm


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
CFLAGS += -DWITH_SONG_SUPPORT
CFLAGS += -DWITH_BISPECTRA

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
	cd $(CLASS_DIR);
	export WITH_BISPECTRA=1;
	export WITH_SONG_SUPPORT=1;
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

song: $(SONG_TOOLS) $(SOURCE_CLASS) $(INPUT2) $(PERTURBATIONS2)\
	$(BESSEL2) $(TRANSFER2) $(BISPECTRA2) $(SPECTRA2) $(OUTPUT) $(SONG)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_params: $(SONG_TOOLS) $(SOURCE_CLASS) $(INPUT2) $(PRINT_PARAMS)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_sources1: $(SONG_TOOLS) $(SOURCE_CLASS) $(INPUT2) $(PERTURBATIONS2) \
	$(PRINT_SOURCES1)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_sources2: $(SONG_TOOLS) $(SOURCE_CLASS) $(INPUT2) $(PERTURBATIONS2) \
	$(PRINT_SOURCES2)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_transfers1: $(SONG_TOOLS) $(SOURCE_CLASS) $(INPUT2) $(PERTURBATIONS2) $(PRINT_TRANSFERS1)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_transfers2: $(SONG_TOOLS) $(SOURCE_CLASS) $(INPUT2) $(PERTURBATIONS2)\
	$(BESSEL2) $(TRANSFER2) $(PRINT_TRANSFERS2)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_bispectra: $(SONG_TOOLS) $(SOURCE_CLASS) $(INPUT2) $(PERTURBATIONS2)\
	$(BESSEL2) $(TRANSFER2) $(BISPECTRA2) $(SPECTRA2) $(OUTPUT) $(PRINT_BISPECTRA)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_background: $(SONG_TOOLS) $(SOURCE_CLASS) $(INPUT2) $(PRINT_BACKGROUND)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_thermo: $(SONG_TOOLS) $(SOURCE_CLASS) $(INPUT2) $(PRINT_THERMO)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_couplings: $(SONG_TOOLS) $(SOURCE_CLASS) $(PRINT_COUPLINGS)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

clean: .base
	rm -rf $(WRKDIR);
	cd $(CLASS_DIR); $(MAKE) clean;

# THE EXECUTABLES BELOW THIS POINT ARE BROKEN AND NEED TO BE UPDATED TO THE
# LATEST VERSION OF SONG!

PRINT_K = print_k.o
PRINT_K_SONG = print_k_song.o
PRINT_TAU_SONG = print_tau_song.o
PRINT_INITIAL_CONDITIONS_SONG = print_initial_conditions_song.o
PRINT_KERNEL = print_kernel.o

print_k: $(SONG_TOOLS) $(SOURCE_CLASS) $(INPUT2) $(PRINT_K)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_k_song: $(SONG_TOOLS) $(SOURCE_CLASS) $(PRINT_K_SONG)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_tau_song: $(SONG_TOOLS) $(SOURCE_CLASS) $(PERTURBATIONS2) $(PRINT_TAU_SONG)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_initial_conditions_song: $(SONG_TOOLS) $(SOURCE_CLASS) $(PERTURBATIONS2) \
	$(PRINT_INITIAL_CONDITIONS_SONG)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm	

print_kernel: $(SONG_TOOLS) $(SOURCE_CLASS) $(PERTURBATIONS2) $(PRINT_KERNEL)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm


