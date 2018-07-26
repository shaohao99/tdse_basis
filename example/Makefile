# $Id: Makefile 558 2007-06-27 18:03:31Z arvid $

npsflib_path=/home/shaohao/program-basis/npsflib-basis

include $(npsflib_path)/common.mk
#include $(petsc_dir)/conf/base

#petsc_link_flag=$(PETSC_SYS_LIB_COMPLEX) $(PETSC_VEC_LIB_COMPLEX) $(PETSC_MAT_LIB_COMPLEX) $(PETSC_DM_LIB_COMPLEX) $(PETSC_KSP_LIB_COMPLEX) $(PETSC_SNES_LIB_COMPLEX) $(PETSC_LIB_COMPLEX)

INCDIR=$(npsflib_path)/$(INCDIRNAME)
LIBDIR=$(npsflib_path)/$(LIBDIRNAME)
LIBSRCDIR=$(npsflib_path)/$(LIBSRCDIRNAME)

headers=$(addprefix $(INCDIR)/, $(headerfilenames))
libraryfiles=$(addprefix $(LIBSRCDIR)/, $(libraryfilenames))

libbasis.a: $(LIBDIR)/libbasis.a

#$(LIBDIR)/libnpsf.a:
#	$(MAKE) -C $(LIBDIR)

main-tdse: main-tdse.cpp libbasis.a
	$(MPICompiler) $(CFLAGS) -I$(INCDIR) main-tdse.cpp -o $@ -L$(LIBDIR) -lbasis -g $(petsc_link_flag) -I$(gslInclude) $(gsllibdir)/libgsl.a
	rm *.o
	$(call gettar)
	$(call lib_msg, $@)

main-dme: main-dme.cpp libbasis.a
	$(MPICompiler) $(CFLAGS) -I$(INCDIR) main-dme.cpp -o $@ -L$(LIBDIR) -lbasis -g $(petsc_link_flag) -I$(gslInclude) $(gsllibdir)/libgsl.a
	rm *.o
	$(call gettar)
	$(call lib_msg, $@)

main-cwf: main-cwf.cpp libbasis.a
	$(MPICompiler) $(CFLAGS) -I$(INCDIR) main-cwf.cpp -o $@ -L$(LIBDIR) -lbasis -g $(petsc_link_flag) -I$(gslInclude) $(gsllibdir)/libgsl.a
	rm *.o
	$(call gettar)
	$(call lib_msg, $@)

doc:
	$(DOXYGEN)

clean: 	
	rm -f *.o
	rm -f $(npsflib_path)/lib/*.a
	rm -f $(npsflib_path)/libsrc/*.a
	rm -f $(npsflib_path)/libsrc/*.o
	rm -f $(npsflib_path)/libsrc/*.d


.PHONY: all clean

