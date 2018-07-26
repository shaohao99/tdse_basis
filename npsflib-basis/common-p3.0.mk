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

omp_include=/usr/local/packages/openmpi/1.3.4/intel-11.1/include

dirCoulombwave=$(npsflibdir)/Coulomb_wave  # for Coulomb wave

MPICompiler=mpicxx # mpicc # mpicxx
LIBBASIS=libbasis
basisfilenames=BasisPropagate.cpp BasisLaser.cpp #PETSc.cpp

petsc_dir=/usr/local/packages/petsc/3.0.0.p10/intel-11.1-mvapich-1.1
#petsc_dir= /usr/local/packages/petsc/3.3-p6/intel-11.1-mvapich2-1.4
petsc_include=$(petsc_dir)/include
petsc_lib=$(petsc_dir)/lib

mpich_dir=/usr/local/packages/mvapich/1.1/intel-11.1
mpich_include=$(mpich_dir)/include
mpich_lib=$(mpich_dir)/lib

cxx_dir=/usr/lib/gcc/x86_64-redhat-linux/3.4.6

petsc_compile_flag = -DPETSC_USE_COMPLEX -DPETSC_CLANGUAGE_CXX -fPIC -xW -I/usr/local/packages/petsc/3.0.0.p10/intel-11.1-mvapich-1.1/include -I/usr/local/packages/SuperLU_DIST/2.3/intel-11.1-mvapich-1.1/include -I/usr/local/packages/ParMetis/3.1.1/intel-11.1-mvapich-1.1/include -I/usr/local/packages/hypre/2.4.0b/intel-11.1-mvapich-1.1/include -I/usr/local/packages/mvapich/1.1/intel-11.1/include
#petsc_compile_flag = -DPETSC_USE_COMPLEX -DPETSC_CLANGUAGE_CXX -fPIC -xW -I/usr/local/packages/petsc/3.0.0.p10/intel-11.1-mvapich-1.1/include -I/usr/local/packages/mvapich/1.1/intel-11.1/include

petsc_link_flag = -DPETSC_USE_COMPLEX -DPETSC_CLANGUAGE_CXX -fPIC -xW -I/usr/local/packages/petsc/3.0.0.p10/intel-11.1-mvapich-1.1/include -Wl,-rpath,/usr/local/packages/petsc/3.0.0.p10/intel-11.1-mvapich-1.1/lib -L/usr/local/packages/petsc/3.0.0.p10/intel-11.1-mvapich-1.1/lib -lpetscsnes -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc -lpetsccontrib -lpetscts  -Wl,-rpath,/usr/local/packages/SuperLU_DIST/2.3/intel-11.1-mvapich-1.1/lib -L/usr/local/packages/SuperLU_DIST/2.3/intel-11.1-mvapich-1.1/lib -lsuperlu_dist_2.3 -Wl,-rpath,/usr/local/packages/ParMetis/3.1.1/intel-11.1-mvapich-1.1/lib -L/usr/local/packages/ParMetis/3.1.1/intel-11.1-mvapich-1.1/lib -lparmetis -lmetis -Wl,-rpath,/usr/local/packages/hypre/2.4.0b/intel-11.1-mvapich-1.1/lib -L/usr/local/packages/hypre/2.4.0b/intel-11.1-mvapich-1.1/lib -lHYPRE -lpmpich++ -lstdc++ -Wl,-rpath,/usr/local/ofed/lib64 -Wl,-rpath,/usr/local/packages/lapack/3.2/intel-11.1-mvapich-1.1/lib -L/usr/local/packages/lapack/3.2/intel-11.1-mvapich-1.1/lib -llapack -lblas -Wl,-rpath,/usr/local/ofed/lib64 -L/usr/local/ofed/lib64 -Wl,-rpath,/usr/local/packages/mvapich/1.1/intel-11.1/lib/shared -L/usr/local/packages/mvapich/1.1/intel-11.1/lib/shared -Wl,-rpath,/usr/local/packages/mvapich/1.1/intel-11.1/lib -L/usr/local/packages/mvapich/1.1/intel-11.1/lib -ldl -lmpichf90nc -lmpichfarg -lmpich -libverbs -libumad -lpthread -lrt -Wl,-rpath,/usr/local/compilers/Intel/intel_fc_11.1/lib/intel64 -L/usr/local/compilers/Intel/intel_fc_11.1/lib/intel64 -Wl,-rpath,/usr/lib/gcc/x86_64-redhat-linux/3.4.6 -L/usr/lib/gcc/x86_64-redhat-linux/3.4.6 -lifport -lifcore -limf -lsvml -lm -lipgo -lirc -lgcc_s -lirc_s -lm -lpmpich++ -lstdc++ -lpmpich++ -lstdc++ -ldl
#petsc_link_flag = -DPETSC_USE_COMPLEX -DPETSC_CLANGUAGE_CXX -fPIC -xW -I/usr/local/packages/petsc/3.0.0.p10/intel-11.1-mvapich-1.1/include -Wl,-rpath,/usr/local/packages/petsc/3.0.0.p10/intel-11.1-mvapich-1.1/lib -L/usr/local/packages/petsc/3.0.0.p10/intel-11.1-mvapich-1.1/lib -lpetscsnes -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc -lpetsccontrib -lpetscts  -Wl,-rpath,/usr/local/packages/SuperLU_DIST/2.3/intel-11.1-mvapich-1.1/lib -L/usr/local/packages/SuperLU_DIST/2.3/intel-11.1-mvapich-1.1/lib -lsuperlu_dist_2.3 -Wl,-rpath,/usr/local/packages/ParMetis/3.1.1/intel-11.1-mvapich-1.1/lib -L/usr/local/packages/ParMetis/3.1.1/intel-11.1-mvapich-1.1/lib -lparmetis -lmetis -Wl,-rpath,/usr/local/packages/hypre/2.4.0b/intel-11.1-mvapich-1.1/lib -L/usr/local/packages/hypre/2.4.0b/intel-11.1-mvapich-1.1/lib -lHYPRE -lpmpich++ -lstdc++ -Wl,-rpath,/usr/local/ofed/lib64 -Wl,-rpath,/usr/local/packages/lapack/3.2/intel-11.1-mvapich-1.1/lib -L/usr/local/packages/lapack/3.2/intel-11.1-mvapich-1.1/lib -llapack -lblas -Wl,-rpath,/usr/local/ofed/lib64 -L/usr/local/ofed/lib64 -Wl,-rpath,/usr/local/packages/mvapich/1.1/intel-11.1/lib/shared -L/usr/local/packages/mvapich/1.1/intel-11.1/lib/shared -Wl,-rpath,/usr/local/packages/mvapich/1.1/intel-11.1/lib -L/usr/local/packages/mvapich/1.1/intel-11.1/lib -ldl -lmpichf90nc -lmpichfarg -lmpich -libverbs -libumad -lpthread -lrt -Wl,-rpath,/usr/local/compilers/Intel/intel_fc_11.1/lib/intel64 -L/usr/local/compilers/Intel/intel_fc_11.1/lib/intel64 -Wl,-rpath,/usr/lib/gcc/x86_64-redhat-linux/3.4.6 -L/usr/lib/gcc/x86_64-redhat-linux/3.4.6 -lifport -lifcore -limf -lsvml -lm -lipgo -lirc -lgcc_s -lirc_s -lm -lpmpich++ -lstdc++ -lpmpich++ -lstdc++ -ldl

#petsc_compile_flag = -DPETSC_USE_COMPLEX -DPETSC_CLANGUAGE_CXX -ansi -Wwrite-strings -Wno-strict-aliasing -g3 -I$(petsc_include) -I$(mpich_include)
#petsc_link_flag = -Wwrite-strings -Wno-strict-aliasing -g3  -Wl,-rpath,/usr/local/petsc-3.0.0-p9/lib -L/usr/local/petsc-3.0.0-p9/lib -lpetscsnes -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc   -lX11 -llapack -lblas -L/usr/local/mpich2-1.0.5p4/lib -L/usr/local/intel/mkl_10.1.0.015/lib/em64t -L/usr/local/intel/fce_11.1.046/lib/intel64 -L/central/intel/mkl_10.1.0.015/lib/em64t -L/central/intel/fc_11.1.046/lib/intel64 -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2 -ldl -lmpich -lrt -luuid -lgcc_s -lmpichf90 -lifport -lifcore -limf -lsvml -lm -lipgo -lirc -lpthread -lirc_s -lm -lmpichcxx -lstdc++ -lmpichcxx -lstdc++ -ldl -lmpich -lrt -luuid -lgcc_s -ldl  $(petsc_include) $(petsc_lib) -DPETSC_USE_COMPLEX -DPETSC_CLANGUAGE_CXX



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

