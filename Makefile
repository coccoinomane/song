# Makefile for Song.
# Created by Guido Walter Pettinari on 2010.

# =========================================================================
# =                         Compilers and flags                           =
# =========================================================================

# Your C compiler and library tool
CC			 = gcc 
AR       = ar rv

# Optimization flags
OPTFLAG    = -O

# Compilation flags
CFLAGS     += -g
# CFLAGS     += -w
CFLAGS     += -std=c99
CFLAGS     += -fopenmp # comment for compiling without parallel support
# CFLAGS     += -DVALGRIND

# Header files and libraries
INCLUDES 					  = -I../include -I../$(CLASS_DIR)/include
LIBRARIES           = -fopenmp -lm


# =========================================================================
# =                        Directories locations                          =
# =========================================================================

# Build directory
MDIR := .
WRKDIR = $(MDIR)/build

# CLASS subfolder
CLASS_DIR = $(MDIR)/class.git

# Source files
vpath %.c source:tools:main:test:$(CLASS_DIR)/source:$(CLASS_DIR)/tools:$(CLASS_DIR)/main:$(CLASS_DIR)/test
# vpath %.c source:tools:main:test:source/no_longer_needed:tools/no_longer_needed:main/no_longer_needed :$(CLASS_DIR)/source:$(CLASS_DIR)/tools:$(CLASS_DIR)/main:$(CLASS_DIR)/test
vpath %.o build:$(CLASS_DIR)/build
vpath .base build

.base:
	if ! [ -a $(WRKDIR) ]; then mkdir $(WRKDIR) ; mkdir $(WRKDIR)/lib; fi;
	touch build/.base

# Make sure to compile CLASS with support for SONG
CFLAGS += -DWITH_SONG_SUPPORT
CFLAGS += -DWITH_BISPECTRA


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
nonlinear.o transfer.o spectra.o lensing.o bessel.o bispectra.o # fisher.o
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
SONG_TOOLS = $(TOOLS) song_tools.o slatec_3j_C.o mesh_interpolation.o
INPUT2 = input2.o
PERTURBATIONS2 = perturbations2.o
BESSEL = bessel.o
BESSEL2 = bessel2.o
TRANSFER2 = transfer2.o
BISPECTRA = bispectra.o
BISPECTRA2 = bispectra2.o
FISHER = fisher.o
UTILITY = utility.o


# ==========================================================================
# =                                Targets                                 =
# ==========================================================================

clean: .base
	rm -rf $(WRKDIR);
	cd $(CLASS_DIR); $(MAKE) clean;

# Rule to be executed when make is called with no parameters
default: class song print_params print_sources1 print_sources2 print_transfers1 print_transfers2 print_bispectra

# CLASS executables
libclass.a class test_background test_thermodynamics test_perturbations test_transfer classy tar: 
	cd $(CLASS_DIR); export WITH_BISPECTRA=1; export WITH_SONG_SUPPORT=1; make $@

# SONG executables
TEST_INPUT = test_input.o
SONG = song.o
PRINT_K = print_k.o
PRINT_K_2ND_ORDER = print_k_2nd_order.o
PRINT_TAU_2ND_ORDER = print_tau_2nd_order.o
PRINT_PARAMS = print_params.o
PRINT_SOURCES1 = print_sources1.o
PRINT_SOURCES2 = print_sources2.o
PRINT_TRANSFERS1 = print_transfers1.o
PRINT_TRANSFERS2 = print_transfers2.o
PRINT_CLS = print_cls.o
PRINT_BISPECTRA = print_bispectra.o
PRINT_BACKGROUND = print_background.o
PRINT_THERMO = print_thermo.o
PRINT_INITIAL_CONDITIONS_2ND_ORDER = print_initial_conditions_2nd_order.o
PRINT_KERNEL = print_kernel.o
PRINT_MATRICES = print_matrices.o
TEST_COUPLINGS = test_couplings.o

song: $(SONG_TOOLS) $(SOURCE_CLASS) $(INPUT2) $(PERTURBATIONS2)\
	$(BESSEL) $(BESSEL2) $(TRANSFER2) $(BISPECTRA) $(BISPECTRA2)\
	$(FISHER) $(UTILITY) $(SONG)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_bispectra: $(SONG_TOOLS) $(SOURCE_CLASS) $(INPUT2) $(PERTURBATIONS2)\
	$(BESSEL) $(BESSEL2) $(TRANSFER2) $(BISPECTRA) test_bispectra.o
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
	$(BESSEL) $(BESSEL2) $(TRANSFER2) $(PRINT_TRANSFERS2)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_cls: $(SONG_TOOLS) $(SOURCE_CLASS) $(BESSEL)\
	$(BISPECTRA) $(FISHER) $(PRINT_CLS)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

# print_bispectra: $(SONG_TOOLS) $(SOURCE_CLASS) $(INPUT2) $(PERTURBATIONS2)\
# 	$(BESSEL) $(BESSEL2) $(TRANSFER2) $(BISPECTRA) $(BISPECTRA2) $(FISHER)\
# 	$(UTILITY) $(PRINT_BISPECTRA)
# 	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_bispectra: $(SONG_TOOLS) $(SOURCE_CLASS) $(INPUT2) $(PERTURBATIONS2)\
	$(BESSEL) $(BESSEL2) $(TRANSFER2) $(BISPECTRA)\
	$(UTILITY) $(PRINT_BISPECTRA)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_background: $(SONG_TOOLS) $(SOURCE_CLASS) $(INPUT2) $(PRINT_BACKGROUND)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_thermo: $(SONG_TOOLS) $(SOURCE_CLASS) $(INPUT2) $(PRINT_THERMO)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_k: $(SONG_TOOLS) $(SOURCE_CLASS) $(INPUT2) $(PRINT_K)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_k_2nd_order: $(SONG_TOOLS) $(SOURCE_CLASS) $(PRINT_K_2ND_ORDER)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_tau_2nd_order: $(SONG_TOOLS) $(SOURCE_CLASS) $(PERTURBATIONS2) $(PRINT_TAU_2ND_ORDER)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_initial_conditions_2nd_order: $(SONG_TOOLS) $(SOURCE_CLASS) $(PERTURBATIONS2) \
	$(PRINT_INITIAL_CONDITIONS_2ND_ORDER)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm	

print_kernel: $(SONG_TOOLS) $(SOURCE_CLASS) $(PERTURBATIONS2) $(PRINT_KERNEL)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_couplings: $(SONG_TOOLS) $(SOURCE_CLASS) $(TEST_COUPLINGS)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

