SHELL=/bin/bash

npsflibdir=/home/shaohao/program-basis/npsflib-basis

MATLABDIR=/usr/local/matlab/2008b
MATCPPOPTS=-I$(MATLABDIR)/extern/include -L$(MATLABDIR)/bin/glnxa64
MATLIBDIR=$(npsflibdir)/Matlab2008blib
MATLIBS=-lmat -lmx -lut -licuio -licuuc -licudata -licui18n -lhdf5 -lmwfl -lexpat -lboost_thread-gcc41-mt-1_34_1 -lboost_signals-gcc41-mt-1_34_1  # for matlab2008b

#OPTFLAGS=-O4 -march=k8 -m64 -msse -msse2 -fomit-frame-pointer -mfpmath=sse -funroll-loops -fopenmp
#OPTFLAGS=-O2 -march=k8 -m64 -fomit-frame-pointer -mfpmath=sse -funroll-loops
OPTFLAGS=-O3 -m64 -fomit-frame-pointer -funroll-loops
CFLAGS= -fno-common $(OPTFLAGS)

LIBSRCDIRNAME=libsrc
LIBDIRNAME=lib
INCDIRNAME=include

FFTWDIRNAME=/usr/local/packages/fftw/3.2/gcc-4.3.2-mvapich-1.1/include   # for fftw-3.2
FFTWLIBDIR=/usr/local/packages/fftw/3.2/gcc-4.3.2-mvapich-1.1/lib

gslInclude=/usr/local/packages/gsl/1.9/intel-11.1-mvapich-1.1/include
gsllibdir=/usr/local/packages/gsl/1.9/intel-11.1-mvapich-1.1/lib

#omp_include=/usr/local/packages/openmpi/1.3.4/intel-11.1/include

dirCoulombwave=$(npsflibdir)/Coulomb_wave  # for Coulomb wave

MPICompiler=mpicxx # mpicc # mpicxx
LIBBASIS=libbasis
basisfilenames=BasisPropagate.cpp BasisLaser.cpp #PETSc.cpp

petsc_dir=/usr/local/packages/petsc/3.4.0/intel-11.1-mvapich-1.1
petsc_include=$(petsc_dir)/include
petsc_lib=$(petsc_dir)/lib

#mpich_dir=/usr/local/packages/mvapich/1.1/intel-11.1
#mpich_include=$(mpich_dir)/include
#mpich_lib=$(mpich_dir)/lib

#cxx_dir=/usr/lib/gcc/x86_64-redhat-linux/3.4.6

petsc_compile_flag = -DPETSC_USE_COMPLEX -DPETSC_CLANGUAGE_CXX -fPIC -wd1572 -g   -I/usr/local/packages/petsc/3.4.0/intel-11.1-mvapich-1.1/include -I/usr/local/packages/petsc/3.4.0/intel-11.1-mvapich-1.1/include -I/usr/X11R6/include -I/usr/include    -D__INSDIR__=src/vec/vec/examples/tutorials/

petsc_link_flag = -DPETSC_USE_COMPLEX -DPETSC_CLANGUAGE_CXX -fPIC -wd1572 -g -Wl,-rpath,/usr/local/packages/petsc/3.4.0/intel-11.1-mvapich-1.1/lib -L/usr/local/packages/petsc/3.4.0/intel-11.1-mvapich-1.1/lib  -lpetsc -Wl,-rpath,/usr/local/packages/lapack/3.4.2/intel-11.1/lib -L/usr/local/packages/lapack/3.4.2/intel-11.1/lib -llapack -lblas -Wl,-rpath,/usr/X11R6/lib64 -L/usr/X11R6/lib64 -lX11 -Wl,-rpath,/usr/lib -L/usr/lib -lexpat -lpthread -Wl,-rpath,/usr/local/ofed/lib64 -L/usr/local/ofed/lib64 -Wl,-rpath,/usr/local/packages/mvapich/1.1/intel-11.1/lib/shared -L/usr/local/packages/mvapich/1.1/intel-11.1/lib/shared -Wl,-rpath,/usr/local/packages/mvapich/1.1/intel-11.1/lib -L/usr/local/packages/mvapich/1.1/intel-11.1/lib -Wl,-rpath,/usr/local/compilers/Intel/intel_cc_11.1/lib/intel64 -L/usr/local/compilers/Intel/intel_cc_11.1/lib/intel64 -Wl,-rpath,/home/compilers/GNU/gcc-4.3.2/lib/gcc/x86_64-unknown-linux-gnu/4.3.2 -L/home/compilers/GNU/gcc-4.3.2/lib/gcc/x86_64-unknown-linux-gnu/4.3.2 -Wl,-rpath,/home/compilers/GNU/gcc-4.3.2/lib/gcc -L/home/compilers/GNU/gcc-4.3.2/lib/gcc -Wl,-rpath,/home/compilers/GNU/gcc-4.3.2/lib64 -L/home/compilers/GNU/gcc-4.3.2/lib64 -Wl,-rpath,/home/compilers/GNU/gcc-4.3.2/lib -L/home/compilers/GNU/gcc-4.3.2/lib -lmpichf90nc -lmpichfarg -Wl,-rpath,/usr/local/compilers/Intel/intel_fc_11.1/lib/intel64 -L/usr/local/compilers/Intel/intel_fc_11.1/lib/intel64 -lifport -lifcore -lm -lm -lpmpich++ -lstdc++ -lpmpich++ -lstdc++ -ldl -lmpich -libverbs -libumad -lpthread -lrt -limf -lsvml -lipgo -ldecimal -lgcc_s -lirc -lirc_s -ldl -I/usr/local/packages/petsc/3.4.0/intel-11.1-mvapich-1.1/include



define lib_msg
	@echo "--------------------------------------------------------------"
	@echo "    To run $(1), please set (just one time)"
	@echo "    For bash:"
	@echo 
	@echo '  source /data/npsf/share/setenv.sh'
	@echo ""
	@echo "    For csh:"
	@echo 
	@echo '  source /data/npsf/share/setenv.csh'
	@echo
	@echo "(you may also consider putting that line into either .bash_profile or .login)"
endef

define gettar
	@echo "(should be) Creating tarball..."
	@echo "(but it's actually not done at this time.. fix it in common.mk)"
endef

