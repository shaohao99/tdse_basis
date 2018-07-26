# $Id: Makefile 447 2006-08-31 16:38:09Z fhe $
include ../common.mk

#include ${petsc_dir}/conf/base

.PHONY: all

# make single simple rule for creating each target 
# to have all source files listed in one rule
# to make tar file out of them
#
# here specify the default target:

INCDIR=../$(INCDIRNAME)
FFTWDIR=$(FFTWDIRNAME)

headers=$(addprefix $(INCDIR)/, $(headerfilenames))

all: libbasis.a

libbasis.a: $(basisfilenames:.cpp=.o)
	ar rs $@ $^

# how to build a .o file from a .cpp file
%.o: %.cpp
	$(MPICompiler) $(CFLAGS) -I$(INCDIR) -I$(gslInclude) $(petsc_compile_flag) -c $<
#	$(MPICompiler) $(CFLAGS) -I$(INCDIR) -I$(gslInclude) -I$(omp_include) $(petsc_compile_flag) -c $<
#       $(MPICompiler) -I$(INCDIR) -I$(FFTWDIR) -I$(dirCoulombwave) $(PETSC_KSP_LIB) -c $<

# auto-generate the dependency tree
include $(subst .cpp,.d,$(parallelfilenames))

%.d: %.cpp
	$(MPICompiler) -M -I$(INCDIR) -I$(gslInclude) -I$(omp_include) $(petsc_compile_flag) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

# auto-depend done

doc:
	$(DOXYGEN)

cleaning:
	@rm -f *.o *.a *.d

.PHONY: all clean
