// #include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <sys/stat.h>
 
#include "BasisPropagate.h"
#include "BasisLaser.h"
#include "BasisConstant.h"
#include "petscksp.h"

using namespace std;
using namespace Basis_Propagate;
using namespace Basis_Laser;


int main( int argc , char *argv[] )
{
  
// **************  Initialize Petsc and MPI  *************  //
  static char help[] = "Parallel propagation for 3D TDSE on an eigen-energy basis by Crank-Nicoleson method.\n\n";
  PetscInitialize(&argc, &argv, (char *)0, help);   // MPI_Initialize is also called
  PetscMPIInt    rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

//  declare
  int n_basis=1, n_dip=1;
  double lambda1=800., i1=3.e12, n_c1=11., cep1_pi=0.;
  double w2=20., i2=1.e10, n_c2=5., cep2_pi=0., delay=0., dt=0.05, t_wait=1.;
  int nenv=2, iend_pulses;
  double Ipot=21.5645, E_adj=1.5;  // Ne
  vector<double> time, E1, E2, Efield, Afield;
  string output_folder, eng_file, dipole_file;

// input data
  if(rank==0) cout<<"--- start input data --- "<<endl;
  ifstream input_data;
  Open_inFile(input_data, "inp-tdse");
  input_data >> lambda1 >> i1 >> n_c1 >> cep1_pi >> nenv;
  input_data >> w2 >> i2 >> n_c2 >> cep2_pi;
  input_data >> delay >> dt >> t_wait;
  input_data >> n_basis >> n_dip;
  input_data >> E_adj >> Ipot;
  input_data >> eng_file;
  input_data >> dipole_file;
  input_data >> output_folder;

  const char* dir_name = output_folder.c_str();
  int dir_exist = mkdir (dir_name, S_IRWXU);  // create output folder
  if(rank==0){
   cout<< "number of basis: "<<n_basis <<endl;
   cout<< "created output folder: " << output_folder<<endl;
  }

// read energy and diploe from grasp
  ifstream input_eng;
  Open_inFile(input_eng, eng_file);
  vector<int> idx(n_basis);
  vector<double> E_levels(n_basis);
  for(int i=0; i<n_basis; i++) { 
    input_eng >> idx[i] >> E_levels[i]; 
    //if(rank==0) cout<<"energy"<<' '<<i<<' '<<idx[i]<<' '<<E_levels[i]<<endl;
  }

  ifstream input_dip;
  Open_inFile(input_dip, dipole_file);
  vector<int> idx_high(n_dip), idx_low(n_dip);
  vector<double> E_tran(n_dip), oscl(n_dip), mat_ele(n_dip), dipole(n_dip);
  for(int i=0; i<n_dip; i++) { 
    input_dip >> idx_high[i] >> idx_low[i] >> E_tran[i] >> oscl[i] >> mat_ele[i];
    if(idx_low[i]==1){  // ground - excited, adjust energy
      dipole[i]=sqrt( oscl[i]/(2.*(E_tran[i]/tran_au_wn+E_adj/tran_au_eV)) ); // in a.u.
    }
    else{  // between excited states
      dipole[i]=sqrt(oscl[i]/(2.*E_tran[i]/tran_au_wn)); // in a.u.
    }
    if(mat_ele[i]<0.) dipole[i] *= -1.;
    //if(rank==0) cout<<"dipole"<<' '<<idx_high[i]<<' '<<idx_low[i]<<' '<<dipole[i]<<endl;
  }
  

// create vec for wavefunction
  if(rank==0) cout<<"--- Create Vec and Mat ---"<<endl;
  Vec psi;
  Create_Vec(psi, n_basis);

// get local info for every processors
  vector<int> nlocal;
  Get_Local_Info (psi, nlocal);
   cout<<"nlocal="<<nlocal[0]<<'\t'<<"rstart="<<nlocal[1]<<'\t'<<"rend="<<nlocal[2]<<endl;

// vec for diag eigen energies
  Vec E_eigen;
  VecDuplicate(psi, &E_eigen);  // duplicate value=0 and local info

// set initial wave function as gound state 1s: the first element (indexed 0) is 1, the others are 0.
  VecSetValue(psi, 0, 1., INSERT_VALUES);  
  Assembly_vector(psi);  // assemble after set value before used

// create Mat for Hamiltonian
  Mat H_dipole;
  Create_Mat ( H_dipole, nlocal, n_basis, n_basis);  // associated with vec psi by nlocal
  MatSetUp(H_dipole);
//  Mat H_temp;

// read dme from disk and set values
    if(rank==0) cout<<"--- setting dme ---"<<endl;
  Set_DME_GRASP (n_basis, n_dip, idx_high, idx_low, dipole, H_dipole);
    if(rank==0) cout<<"--- end calculate dipole moments ---"<<endl;
  //MatSetOption(H_dipole, MAT_SYMMETRY_ETERNAL, PETSC_TRUE);  // set symmetric after setting vales before assembling
  Assembly_matrix_final ( H_dipole ); // finally assemble after setting vales before using
   PetscViewer viewer_Hdip;  // vec view by multi processors
   Create_Viewer ( viewer_Hdip, output_folder + "/h_dipole.dat" );
   MatView(H_dipole, viewer_Hdip);


// set eigen energy 
  PetscScalar  value=complexd(-Ipot/tran_au_eV,0.);  // initio complex
  VecSetValue(E_eigen, 0, value, INSERT_VALUES);
  for(int i=1; i<n_basis; i++){
     value = (E_levels[i]-E_levels[0]-E_adj-Ipot)/tran_au_eV ;     // set eigen energy, in a.u.
     VecSetValue(E_eigen, i, value, INSERT_VALUES);
  }  
  Assembly_vector(E_eigen);  // assemble after set value before used

// observe initial energy and state
  double E_init = Obs_Energy(psi, E_eigen);
  if(rank==0) { 
    cout<<"initial energy = "<< E_init <<endl;  // E = <psi| H |psi>
    //VecView(psi, PETSC_VIEWER_STDOUT_WORLD);
  }

// Laser field
  if(rank==0) cout<<"--- Laser pulse ---"<<endl;
  //One_pulse ( wavelength, intensity, n_cycle, dt, cep, t_wait, time, Efield, Afield, period);
  IR_SAP ( lambda1, i1, n_c1, cep1_pi, w2, i2, n_c2, cep2_pi, nenv, delay, dt, t_wait, time, E1, E2, Efield, Afield, iend_pulses);

  ofstream  output_field;
  Open_outFile( output_field, output_folder + "/laser.dat");

// time propagation
  KSP  ksp;
  PC   pc;
  KSPCreate(PETSC_COMM_WORLD, &ksp);

  if(rank==0) cout<<"--- start time propagation ---"<<endl;
  PetscScalar cn_factor = -0.5*I*dt;  //  -i*dt/2 , for dt[0]=dt[1]
    if(rank==0) cout<<"dt="<<dt<<",  scalar="<<cn_factor*value<<endl;
    PetscViewer viewer_t2;  // vec view by multi processors
    Create_Viewer ( viewer_t2, output_folder + "/h_test.dat" );
    MatScale(H_dipole, cn_factor);
    MatView(H_dipole, viewer_t2);

  PetscScalar Emid, E1mid, E2mid;  // complex by default
  double Ef,E1f,E2f;
  int nt = time.size();

  PetscViewer viewer_test;    
  Create_Viewer ( viewer_test, output_folder + "/psi_vec_t.dat" );  // view psi

  for(int i=0; i<nt; i++){
    
    if(i==nt-1){ 
      Ef = Efield[i];
      E1f = E1[i];
      E2f = E2[i];
    }
    else{
      Ef = (Efield[i+1]+Efield[i])/2.;  // average value of two adjacent points
      E1f = (E1[i+1]+E1[i])/2.;  
      E2f = (E2[i+1]+E2[i])/2.;  
    }
    Emid = complexd(Ef,0.);
    E1mid = complexd(E1f,0.);
    E2mid = complexd(E2f,0.);
   
    PropagateCN_linear_PETSc (cn_factor, Emid, E_eigen, H_dipole, ksp, pc, psi, viewer_test, rank);

    if(i%100==0){
      double norm_t = Obs_Norm(psi);
      if(rank==0) cout<< time[i] <<' '<< norm_t <<endl;
    }
    if(i%2==0){
      if(rank==0) output_field<<time[i] <<' '<<E1[i]<<' '<<E2[i]<<' '<<Efield[i]<<' '<<Afield[i]<<endl;
    }
  }

  if(rank==0) cout<<"--- end time propagation ---"<<endl;
  KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);   // view information for KSP

// test output
//   PetscViewer viewer_htemp;  // vec view by multi processors
//   Create_Viewer ( viewer_htemp, output_folder + "/h_temp.dat" );
//   MatView(H_temp, viewer_htemp);

// output final spectrum
   PetscScalar  *E_out, *cn2_out;  // complex by defaut
   vector<double> array_E(n_basis), array_cn2(n_basis);

   if(rank==0)  cout<<"--- output energy ---"<<endl;
   PetscViewer viewer_E;  // vec view by multi processors
   Create_Viewer ( viewer_E, output_folder + "/E_eigen_vec.dat" );
   VecView(E_eigen, viewer_E);

/*  VecGetArray(E_eigen, &E_out);   // get array 
   for(int i=0; i<n_basis; i++){
      if( nlocal[1] <= i && i < nlocal[2] ){ array_E[i] = real(E_out[i]); }   // store in each processor
   }
   VecRestoreArray(E_eigen, &E_out); //no longer access to array value, MUST restore it to vec, otherwise vec is changed
*/

   if(rank==0) cout<<"--- output spectrum ---"<<endl;
   PetscViewer viewer_psi;    
   Create_Viewer ( viewer_psi, output_folder + "/psi_vec.dat" );  // view psi
   VecView(psi, viewer_psi);

   Vec vec_cn2;    // view |psi|^2
   VecDuplicate(psi, &vec_cn2);  // duplicate value=0 and local info
   VecAbs(psi);   // psi = |c_n|
   VecPointwiseMult(vec_cn2, psi, psi);  // vec_cn2 = |c_n|^2
   PetscViewer viewer_cn2;    // vec view by multi processors
   Create_Viewer ( viewer_cn2, output_folder + "/cn2_vec.dat" );
   VecView(vec_cn2, viewer_cn2);

/*   VecGetArray(vec_cn2, &cn2_out);   // get array
   for(int i=0; i<n_basis; i++){
      if( nlocal[1] <= i && i < nlocal[2] ){ array_cn2[i] = real(cn2_out[i]); }  // store in each processor
   }
   VecRestoreArray(vec_cn2, &cn2_out); //no longer access to array value, MUST restore it to vec, otherwise vec is changed
*/

// finalize
  KSPDestroy(&ksp);
  MatDestroy(&H_dipole);
  VecDestroy(&psi);  
  VecDestroy(&E_eigen);  
  PetscViewerDestroy(&viewer_Hdip);
  PetscViewerDestroy(&viewer_test);
  PetscViewerDestroy(&viewer_E);
  PetscViewerDestroy(&viewer_psi);
  PetscViewerDestroy(&viewer_cn2);

  PetscFinalize();

  if(rank==0) cout << "---program end---" << endl;

}  // end of main

