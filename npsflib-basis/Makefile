# $Id: Makefile 512 2006-11-21 15:22:16Z goerke $
include common.mk

.PHONY: all doc test dist distclean message example

relpathheaders=$(addprefix $(INCDIR)/, $(headers))
relpathlibraryfiles=$(addprefix $(LIBSRCDIR)/, $(libraryfiles))

# how to build the stuff in the lib directory:
SUBDIRS = lib
.PHONY: $(SUBDIRS)
     
#$(SUBDIRS):
#	$(MAKE) -C $@

new:
#	$(MAKE) -C $(LIBSRCDIRNAME)
	$(MAKE) -C $(LIBBASIS)
	$(MAKE) clean -C $(SUBDIRS)
	$(MAKE) -C $(SUBDIRS)
 
clean:
	rm $(LIBBASIS)/*.o
	rm $(LIBBASIS)/*.a
	rm $(LIBDIRNAME)/*.a
#	$(MAKE) $@ -C $(SUBDIRS)
#	$(MAKE) $@ -C $(LIBSRCDIRNAME)
    
#-----------------End of general section---------------------------------------

# make single simple rule for creating each target 
# to have all source files listed in one rule
# to make tar file out of them
#
# here specify the default target:
all: example

# what target shall be build as an example ?
example: sim/car3D/car
# and how to create it:
sim/car3D/car:
	$(MAKE) -C sim/car3D car

# example simulations
example_car3D_filenames = car.cfg Main_cylindrical_3D.cpp
car3Ddir=sim/car3D
example_cyl3D_filenames = example.cfg Main_cylindrical_laser.cpp
cyl3Ddir=sim/cyl3D
# and all files together
simulations=$(addprefix $(car3Ddir), $(example_car3D_filenames)) $(addprefix $(cyl3Ddir), $(example_cyl3D_filenames))
distfiles=$(addprefix $(INCDIR),$(headerfilenames)) $(addprefix libsrc/,$(libraryfilenames)) ${simulations} ${docfilenames} Makefile 

SRCDIRS = libsrc include sim
  
tarball:
#this is command to make tarfile
#- in front says not to stop proceeding if it fails
#shopt -s nullglob says to shell to remove wildcards if don't match
# first remove file in case there is another file hardlinked to it.
	@echo Creating tarball...
	-shopt -s nullglob; rm -f all_sources.tar.gz; tar -czf all_sources.tar.gz --exclude="*.[o,a]" Makefile *.{cpp,h,c,f,cfg} $(SRCDIRS)

test:	
	test/check

doc:
	$(DOXYGEN)

dist:	
	svn export http://planet.pks.mpg.de/svn/npsflib/tags/npsflib-0.2 dist
	tar -cf npsflib-0.2.tar dist && gzip -f npsflib-0.2.tar

distclean: 
	rm -f -r dist

