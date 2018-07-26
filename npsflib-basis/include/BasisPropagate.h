#ifndef BasisPropagate_H
#define BasisPropagate_H

#include <math.h>
#include <complex>
#define complexd complex<double>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

#include "petscts.h"
#include "petscksp.h"
#include "petscvec.h"
#include "petscmat.h"
#include "petscis.h"
#include "petscsys.h"

using namespace std; // to save the std:: in front of vector and cout

namespace Basis_Propagate{

 //const double pi = 3.141592653589793238463383;
 //const complexd I=complexd(0.,1.); 
 //PetscScalar  P_real1=complexd(1.,0.);
 //PetscScalar  P_realm1=complexd(-1.,0.);
 //PetscScalar  P_real2=complexd(2.,0.);

 void Open_outFile( ofstream &filename, string name );
 void Open_inFile( ifstream &filename, string name );

 int index3( int k , int j , int i, int nj, int ni );
 int index2( int a , int b, int nb );
 int index_wf( int a, vector<int> na, int b );
 void int_to_string(int i, string &a);

 void Spherical_grids ( double R0, double dr, double dtheta, double dphi, vector<double> &r, vector<double> &theta, vector<double> &phi);
 void Basis_struct ( int n_max, double k_min, double k_max, double dk, vector<double> &k, vector<int> &n_basis);
 void DME_analytical_linear (vector<double> r, int n_max, vector<int> n_basis, vector<double> k, Vec &diag, Mat &H_dipole, int normalize);
 // void DME_linear (int n_basis, vector<int> l, vector<int> nwf, vector<double> r_max, vector<double> &r, vector<double> &wf, double r_max_new, double dr, Mat &H_dipole);
 void DME_linear (int n_basis, vector<int> l, vector<int> nwf, vector<double> r_max, vector<double> &r, vector<double> &wf, double r_max_new, double dr, string folder);
 void Set_DME (int n_basis, vector<int> l, Mat &H_dipole, string folder);
 void Set_DME_GRASP (int n_basis, int n_dip, vector<int> idx_high, vector<int> idx_low, vector<double> dipole, Mat &H_dipole);
 //void PropagateCN_linear_PETSc ( PetscScalar cn_factor, PetscScalar Et, Vec diag, Mat H_dipole, KSP &ksp, PC &pc, Vec &psi );
  void PropagateCN_linear_PETSc ( PetscScalar cn_factor, PetscScalar Et, Vec diag, Mat H_dipole, KSP &ksp, PC &pc, Vec &psi, PetscViewer &viewer_test, PetscMPIInt rank );
 //void PropagateCN_linear_PETSc ( PetscScalar cn_factor, PetscScalar Et, Vec diag, Mat H_dipole, KSP &ksp, PC &pc, Vec &psi, Mat &H_temp);

 void Create_Vec (Vec &vec, int nsize);
 void Get_Local_Info ( Vec vec, vector<int> &nlocal);
 void Create_Mat ( Mat &matrix, vector<int> nlocal, PetscInt n, PetscInt m );
 void Assembly_vector ( Vec &vec );
 void Assembly_matrix_final ( Mat &matrix );
 void Assembly_matrix_flush ( Mat &matrix );
 void Create_Viewer (PetscViewer &viewer, string filename );

 double Obs_Norm ( Vec psi  );
 PetscReal Obs_Norm_sqrt ( Vec psi  );
 void Normalize( Vec &psi );
 double Obs_Energy ( Vec psi, Vec E_eigen );

 inline complexd Gamma_Lanczos (complexd z);
 complexd Y_lm_plgndr ( double theta, double phi, int l, int m);
 double plgndr(int l, int m, double x);
 inline double Factorial(int n);

  void Mom2D_cwf ( double k_max, double dk, int n_theta, vector<double> wf_new, double r_max, double dr, int n0, int num_basis, int l_max, vector<int> l_basis, vector<int> n_basis, vector< vector<double> > c_nl, vector<double> &k, vector<double> &theta, vector<double> &P_mom, vector<double> &dP_dk_dkrho);


  void DME_index_AIS1(int nE, int n_isolate, PetscInt *i1, PetscInt *i2, PetscInt *j2, PetscInt *j3, PetscInt *j4, PetscInt *j5, PetscInt *j6, PetscInt *j7, PetscInt *i3, PetscInt *i4, PetscInt *i5, PetscInt *i6);

  void DME_init_AIS1(int nE, vector<double> E_eigen_array, PetscInt *i1, PetscInt *i2, PetscInt *j2, PetscInt *j3, PetscInt *j4, PetscInt *j5, PetscInt *j6, PetscInt *j7, PetscInt *i3, PetscInt *i4, PetscInt *i5, PetscInt *i6, Mat & H_temp);

  void DME_IS_AIS1(int nE, PetscInt *i1, PetscInt *i2, PetscInt *j2, PetscInt *j3, PetscInt *j4, PetscInt *j5, PetscInt *j6, PetscInt *j7, PetscInt *i3, PetscInt *i4, PetscInt *i5, PetscInt *i6, IS &isi1, IS &isi2, IS &isi3, IS &isi4, IS &isi5, IS &isi6, IS &isj2, IS &isj3, IS &isj4, IS &isj5, IS &isj6, IS &isj7);

  void DME_block_AIS1(Mat Ham, IS isi1, IS isi2, IS isi3, IS isi4, IS isi5, IS isi6, IS isj2, IS isj3, IS isj4, IS isj5, IS isj6, IS isj7, Mat &mat_i1j2, Mat &mat_i1j3, Mat &mat_i1j5, Mat &mat_i1j7, Mat &mat_i2j4, Mat &mat_i2j6, Mat &mat_i3j4, Mat &mat_i4j5, Mat &mat_i5j6, Mat &mat_i6j7, Mat &mat_j2i1, Mat &mat_j3i1, Mat &mat_j5i1, Mat &mat_j7i1, Mat &mat_j4i2, Mat &mat_j6i2, Mat &mat_j4i3, Mat &mat_j5i4, Mat &mat_j6i5, Mat &mat_j7i6);

  void DME_time_AIS1(Mat &H_temp, PetscScalar EXt, PetscScalar ELt, int nE, PetscInt *i1, PetscInt *i2, PetscInt *j2, PetscInt *j3, PetscInt *j4, PetscInt *j5, PetscInt *j6, PetscInt *j7, PetscInt *i3, PetscInt *i4, PetscInt *i5, PetscInt *i6, Mat mat_i1j2, Mat mat_i1j3, Mat mat_i1j5, Mat mat_i1j7, Mat mat_i2j4, Mat mat_i2j6, Mat mat_i3j4, Mat mat_i4j5, Mat mat_i5j6, Mat mat_i6j7, Mat mat_j2i1, Mat mat_j3i1, Mat mat_j5i1, Mat mat_j7i1, Mat mat_j4i2, Mat mat_j6i2, Mat mat_j4i3, Mat mat_j5i4, Mat mat_j6i5, Mat mat_j7i6);

  //void Propagate_APT_AIS1 ( PetscScalar cn_factor, Mat Ham, KSP &ksp, PC &pc, Vec &psi );
  void Propagate_APT_AIS1 ( PetscScalar cn_factor, Mat H_temp, KSP &ksp, PC &pc, Vec &psi );

// --- array set mat
  void DME_array_index_AIS1(int nE, int n_isolate, PetscInt i1[], PetscInt i2[], PetscInt j2[], PetscInt j3[], PetscInt j4[], PetscInt j5[], PetscInt j6[], PetscInt j7[], PetscInt i3[], PetscInt i4[], PetscInt i5[], PetscInt i6[]);
  void DME_array_static_AIS1(int nE, vector<double> E_eigen_array, PetscInt *i2, PetscInt *j5, vector<complexd> ai2j5, Mat &H_temp);
  void DME_array_values_AIS1(int nE, PetscScalar vi1j2[], PetscScalar vi1j3[], PetscScalar vi1j5[], PetscScalar vi1j7[], PetscScalar vi2j4[], PetscScalar vi2j6[], PetscScalar vi2j5[], PetscScalar vi3j4[], PetscScalar vi4j5[], PetscScalar vi5j6[], PetscScalar vi6j7[]);
 void DME_pointer_init_AIS1(int nE, vector<double> E_eigen_array, PetscInt *i1, PetscInt *i2, PetscInt *j2, PetscInt *j3, PetscInt *j4, PetscInt *j5, PetscInt *j6, PetscInt *j7, PetscInt *i3, PetscInt *i4, PetscInt *i5, PetscInt *i6, vector<complexd> ai1j2, vector<complexd> ai1j3, vector<complexd> ai1j5, vector<complexd> ai1j7, vector<complexd> ai2j4,  vector<complexd> ai2j6,  vector<complexd> ai2j5, vector<complexd> ai3j4, vector<complexd> ai4j5, vector<complexd> ai5j6, vector<complexd> ai6j7, Mat &H_temp);
 void DME_array_init_AIS1(int nE, vector<double> E_eigen_array, PetscInt i1[], PetscInt i2[], PetscInt j2[], PetscInt j3[], PetscInt j4[], PetscInt j5[], PetscInt j6[], PetscInt j7[], PetscInt i3[], PetscInt i4[], PetscInt i5[], PetscInt i6[], PetscScalar vi1j2[], PetscScalar vi1j3[], PetscScalar vi1j5[], PetscScalar vi1j7[], PetscScalar vi2j4[], PetscScalar vi2j6[], PetscScalar vi2j5[], PetscScalar vi3j4[], PetscScalar vi4j5[], PetscScalar vi5j6[], PetscScalar vi6j7[], Mat &H_temp);
 void DME_pointer_init_AIS1(int nE, vector<double> E_eigen_array, PetscInt *i1, PetscInt *i2, PetscInt *j2, PetscInt *j3, PetscInt *j4, PetscInt *j5, PetscInt *j6, PetscInt *j7, PetscInt *i3, PetscInt *i4, PetscInt *i5, PetscInt *i6, vector<complexd> ai1j2, vector<complexd> ai1j3, vector<complexd> ai1j5, vector<complexd> ai1j7, vector<complexd> ai2j4,  vector<complexd> ai2j6,  vector<complexd> ai2j5, vector<complexd> ai3j4, vector<complexd> ai4j5, vector<complexd> ai5j6, vector<complexd> ai6j7, Mat &H_temp);
 void DME_array_time_AIS1(PetscScalar EXt, PetscScalar ELt, int nE, PetscInt i1[], PetscInt i2[], PetscInt j2[], PetscInt j3[], PetscInt j4[], PetscInt j5[], PetscInt j6[], PetscInt j7[], PetscInt i3[], PetscInt i4[], PetscInt i5[], PetscInt i6[], PetscScalar vi1j2[], PetscScalar vi1j3[], PetscScalar vi1j5[], PetscScalar vi1j7[], PetscScalar vi2j4[], PetscScalar vi2j6[], PetscScalar vi3j4[], PetscScalar vi4j5[], PetscScalar vi5j6[], PetscScalar vi6j7[], Mat &H_temp);
 void DME_value_time_AIS1(PetscScalar EXt, PetscScalar ELt, int nE, PetscInt i1[], PetscInt i2[], PetscInt j2[], PetscInt j3[], PetscInt j4[], PetscInt j5[], PetscInt j6[], PetscInt j7[], PetscInt i3[], PetscInt i4[], PetscInt i5[], PetscInt i6[], PetscScalar vi1j2[], PetscScalar vi1j3[], PetscScalar vi1j5[], PetscScalar vi1j7[], PetscScalar vi2j4[], PetscScalar vi2j6[], PetscScalar vi3j4[], PetscScalar vi4j5[], PetscScalar vi5j6[], PetscScalar vi6j7[], Mat &H_temp);
 void DME_pointer_time_AIS1(PetscScalar EXt, PetscScalar ELt, int nE, PetscInt *i1, PetscInt *i2, PetscInt *j2, PetscInt *j3, PetscInt *j4, PetscInt *j5, PetscInt *j6, PetscInt *j7, PetscInt *i3, PetscInt *i4, PetscInt *i5, PetscInt *i6, vector<complexd> ai1j2, vector<complexd> ai1j3, vector<complexd> ai1j5, vector<complexd> ai1j7, vector<complexd> ai2j4,  vector<complexd> ai2j6, vector<complexd> ai3j4, vector<complexd> ai4j5, vector<complexd> ai5j6, vector<complexd> ai6j7, Mat &H_temp);
 void DME_full_time_AIS1(PetscScalar EXt, PetscScalar ELt, int nE, vector<double> E_eigen_array, PetscInt *i1, PetscInt *i2, PetscInt *j2, PetscInt *j3, PetscInt *j4, PetscInt *j5, PetscInt *j6, PetscInt *j7, PetscInt *i3, PetscInt *i4, PetscInt *i5, PetscInt *i6, vector<complexd> ai1j2, vector<complexd> ai1j3, vector<complexd> ai1j5, vector<complexd> ai1j7, vector<complexd> ai2j4,  vector<complexd> ai2j6, vector<complexd> ai2j5, vector<complexd> ai3j4, vector<complexd> ai4j5, vector<complexd> ai5j6, vector<complexd> ai6j7, Mat &H_temp);

}

#endif

