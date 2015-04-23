# Makefile for Song.
# Created by Guido Walter Pettinari on 2010.

# =========================================================================
# =                         Compilers and flags                           =
# =========================================================================

# Your C compiler and library tool
CC			 = gcc -g
AR       = ar rv

# Optimization flags
OPTFLAG    = -O

# Openmp flag (comment for compiling without openmp)
CFLAGS     = -std=c99
CFLAGS     = -fopenmp -std=c99

# Header files and libraries
INCLUDES 					  = -I../include
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
vpath %.o build:$(CLASS_DIR)/build
vpath .base build

.base:
	if ! [ -a $(WRKDIR) ]; then mkdir $(WRKDIR) ; mkdir $(WRKDIR)/lib; fi;
	touch build/.base

# Tell CLASS to include all SONG related stuff
export WITH_BISPECTRA = 1
export WITH_SONG_SUPPORT = 1

# =========================================================================
# =                          Compilation rules                            =
# =========================================================================

%.o:  %.c .base
	cd $(WRKDIR); $(CC) $(OPTFLAG) $(CFLAGS) $(INCLUDES) -c ../$< -o $*.o


# ==========================================================================
# =                             Object files                               =
# ==========================================================================

# Source files also present in CLASS
TOOLS = growTable.o dei_rkck.o sparse.o evolver_rkck.o evolver_ndf15.o arrays.o parser.o quadrature.o
INPUT = input.o
BACKGROUND = background.o
THERMO = thermodynamics.o
PERTURBATIONS = perturbations.o
TRANSFER = transfer.o
PRIMORDIAL = primordial.o
SPECTRA = spectra.o
BISPECTRA = bispectra.o
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
BISPECTRA2 = bispectra2.o
FISHER = fisher.o
MAINS = mains.o



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
	cd $(CLASS_DIR); make $@

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

song: $(SONG_TOOLS) $(INPUT) $(INPUT2) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(PERTURBATIONS2)\
	$(BESSEL) $(BESSEL2) $(TRANSFER) $(TRANSFER2) $(PRIMORDIAL) $(SPECTRA) $(BISPECTRA) $(BISPECTRA2)\
	$(FISHER) $(NONLINEAR) $(LENSING) $(OUTPUT) $(MAINS) $(SONG)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_bispectra: $(SONG_TOOLS) $(INPUT) $(INPUT2) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(PERTURBATIONS2)\
	$(BESSEL) $(BESSEL2) $(TRANSFER) $(TRANSFER2) $(PRIMORDIAL) $(SPECTRA) $(BISPECTRA) test_bispectra.o
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_input: $(INPUT) $(INPUT2) $(TEST_INPUT) $(BACKGROUND) $(SONG_TOOLS)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_params: $(SONG_TOOLS) $(INPUT) $(INPUT2) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(PRINT_PARAMS)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_sources1: $(SONG_TOOLS) $(INPUT) $(INPUT2) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(PERTURBATIONS2) \
	$(TRANSFER) $(TRANSFER2) $(PRINT_SOURCES1)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_sources2: $(SONG_TOOLS) $(INPUT) $(INPUT2) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(PERTURBATIONS2) \
	$(PRINT_SOURCES2)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_transfers1: $(SONG_TOOLS) $(INPUT) $(INPUT2) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(PERTURBATIONS2)\
	$(BESSEL) $(BESSEL2) $(TRANSFER) $(TRANSFER2) $(PRINT_TRANSFERS1)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_transfers2: $(SONG_TOOLS) $(INPUT) $(INPUT2) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(PERTURBATIONS2)\
	$(BESSEL) $(BESSEL2) $(TRANSFER) $(TRANSFER2) $(PRINT_TRANSFERS2)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_cls: $(SONG_TOOLS) $(INPUT) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(BESSEL) $(TRANSFER) $(PRIMORDIAL)\
	$(SPECTRA) $(BISPECTRA) $(FISHER) $(LENSING) $(NONLINEAR) $(OUTPUT) $(PRINT_CLS)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_bispectra: $(SONG_TOOLS) $(INPUT) $(INPUT2) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(PERTURBATIONS2)\
	$(BESSEL) $(BESSEL2) $(TRANSFER) $(TRANSFER2) $(PRIMORDIAL) $(SPECTRA) $(BISPECTRA) $(BISPECTRA2) $(FISHER)\
	$(FISHER) $(NONLINEAR) $(LENSING) $(OUTPUT)	$(MAINS) $(PRINT_BISPECTRA)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_background: $(SONG_TOOLS) $(INPUT) $(INPUT2) $(BACKGROUND) $(THERMO) $(PRINT_BACKGROUND)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_thermo: $(SONG_TOOLS) $(INPUT) $(INPUT2) $(BACKGROUND) $(THERMO) $(PRINT_THERMO)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_k: $(SONG_TOOLS) $(INPUT) $(INPUT2) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(PRINT_K)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_k_2nd_order: $(SONG_TOOLS) $(INPUT) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(PRINT_K_2ND_ORDER)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_tau_2nd_order: $(SONG_TOOLS) $(INPUT) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(PRINT_TAU_2ND_ORDER)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_initial_conditions_2nd_order: $(SONG_TOOLS) $(INPUT) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) \
	$(PRINT_INITIAL_CONDITIONS_2ND_ORDER)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm	

print_kernel: $(SONG_TOOLS) $(INPUT) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(PRINT_KERNEL)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

print_matrices: $(SONG_TOOLS) $(INPUT) $(BACKGROUND) $(THREEJ) $(PRINT_MATRICES)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm	

test_couplings: $(SONG_TOOLS) $(INPUT) $(BACKGROUND) $(TEST_COUPLINGS)
	$(CC) $(LIBRARIES) -o  $@ $(addprefix build/,$(notdir $^)) -lm

