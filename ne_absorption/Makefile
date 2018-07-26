# $Id: Makefile 558 2007-06-27 18:03:31Z arvid $

#npsflib_path=/users/becker/shaohao/npsf/npsflib-sh-photon
npsflib_path=/home/shaohao/program-basis/npsflib-basis

include $(npsflib_path)/common.mk

INCDIR=$(npsflib_path)/$(INCDIRNAME)
LIBDIR=$(npsflib_path)/$(LIBDIRNAME)
LIBSRCDIR=$(npsflib_path)/$(LIBSRCDIRNAME)

headers=$(addprefix $(INCDIR)/, $(headerfilenames))
libraryfiles=$(addprefix $(LIBSRCDIR)/, $(libraryfilenames))

libbasis.a: $(LIBDIR)/libbasis.a

#$(LIBDIR)/libnpsf.a:
#	$(MAKE) -C $(LIBDIR)


tdse_grasp: tdse_grasp.cpp libbasis.a
	$(MPICompiler) $(CFLAGS) -I$(INCDIR) tdse_grasp.cpp -o $@ -L$(LIBDIR) -lbasis -g $(petsc_link_flag) -I$(gslInclude) $(gsllibdir)/libgsl.a
	rm *.o
	$(call gettar)
	$(call lib_msg, $@)

main-dme: main-dme.cpp libbasis.a
	$(MPICompiler) $(CFLAGS) -I$(INCDIR) main-dme.cpp -o $@ -L$(LIBDIR) -lbasis -g $(petsc_link_flag) -I$(gslInclude) $(gsllibdir)/libgsl.a
	$(call gettar)
	$(call lib_msg, $@)

main-cwf: main-cwf.cpp libbasis.a
	$(MPICompiler) $(CFLAGS) -I$(INCDIR) main-cwf.cpp -o $@ -L$(LIBDIR) -lbasis -g $(petsc_link_flag) -I$(gslInclude) $(gsllibdir)/libgsl.a
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

