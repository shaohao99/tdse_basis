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
  int n_basis;
  //double r_max_new, dr;
  double wavelength, intensity, n_cycle, cep, t_wait, period;
  vector<double> dt(2), time, Efield, Afield;
  //int n_max=16;
  //double E_min=0., E_max = 0.2, dE = 0.001;
  double n_front, n_back;
  string output_folder, wf_folder, dme_folder;

// input data
  if(rank==0) cout<<"--- start input data --- "<<endl;
  ifstream input_data;
  Open_inFile(input_data, "inp-tdse");
  input_data >> n_basis;
  //input_data >> r_max_new >> dr;
  input_data >> dt[0] >> dt[1];
  input_data >> wavelength>> intensity>> n_cycle>> cep>> t_wait;
  input_data >> n_front >> n_back;
  //input_data >> n_max;
  //input_data >> E_min >> E_max >> dE;
  input_data >> wf_folder;
  input_data >> dme_folder;
  input_data >> output_folder;

  const char* dir_name = output_folder.c_str();
  int dir_exist = mkdir (dir_name, S_IRWXU);  // create output folder
  if(rank==0){
   cout<< "number of basis: "<<n_basis <<endl;
   //cout<< "r_max, dr of new unifrom grids: " << r_max_new <<'\t'<< dr <<endl;
   cout<< "time steps: " << dt[0] <<'\t'<< dt[1]<<endl;
   cout<< "laser parameters: " << wavelength <<'\t'<< intensity <<'\t'<<n_cycle<<'\t'<<cep<<'\t'<<t_wait<<endl;
   cout<< "shape parameters " << n_front <<'\t'<< n_back <<endl;
   //cout<< "maximum of quantum number: "<<n_max <<endl;
   //cout<< "E_min, E_max, dE for energy spectrum: " << E_min <<'\t'<< E_max <<'\t'<< dE <<endl;
   cout<< "created input folder for wf: " << wf_folder <<endl;
   cout<< "created input folder fpr dme: " << dme_folder <<endl;
   cout<< "created output folder: " << output_folder<<endl;
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


    if(rank==0) cout<<"--- read orbit info ---"<<endl;
  ifstream input_obt;
  Open_inFile(input_obt, wf_folder + "/inp-orbit");
  vector<int> nwf(n_basis), l(n_basis), n(n_basis);
  vector<double> E_obt(n_basis), r_max(n_basis);
  for(int i=0; i<n_basis; i++){  // iteration of orbital
     input_obt >> n[i] >> l[i] >> r_max[i] >> nwf[i] >> E_obt[i];
     // cout<< n[i] <<'\t'<< l[i] <<'\t'<< nwf[i] <<'\t'<< E_obt[i] <<endl;
  }

/*
    if(rank==0) cout<<"--- read wf ---"<<endl;
  ifstream input_wf;
  Open_inFile(input_wf, wf_folder + "/wavefunction");
  int n_wf = 0;
  for(int i=0; i<n_basis; i++)
    for (int ir=0; ir<nwf[i]; ir++){ n_wf += 1;  }  // count total number
  vector<double> wf(n_wf), r(n_wf);
  for(int i=0; i<n_basis; i++){  // iteration of orbital
     for (int ir=0; ir<nwf[i]; ir++){
        int itotal = index_wf(i, nwf, ir);
        input_wf >> r[itotal] >> wf[itotal];    // wf = r*psi_orbital
        // cout<< itotal <<'\t'<< r[itotal] <<'\t'<< wf[itotal] <<endl;
     }
  }

    if(rank==0) cout<<"--- start calculating dipole moments ---"<<endl;
  const char* dir_name1 = dme_folder.c_str();
  int dir_exist1 = mkdir (dir_name1, S_IRWXU);  // create output folder
    // DME_linear (n_basis, l, nwf, r_max, r, wf, r_max_new, dr, H_dipole);
  DME_linear (n_basis, l, nwf, r_max, r, wf, r_max_new, dr, dme_folder);
  wf.resize(0);  r.resize(0);    // release memory occupied by old basis wavefunctions

 return 0;
*/

// read dme from disk and set values
    if(rank==0) cout<<"--- setting dme ---"<<endl;
  Set_DME (n_basis, l, H_dipole, dme_folder);
    if(rank==0) cout<<"--- end calculate dipole moments ---"<<endl;
    // MatSetOption(H_dipole, MAT_SYMMETRIC, PETSC_TRUE);  // set symmetric before assembling, not working
  Assembly_matrix_final ( H_dipole ); // finally assemble after setting vales before using

// set eigen energy 
  for(int i=0; i<n_basis; i++){
     PetscScalar  value = E_obt[i];     // set eigen energy
     VecSetValue(E_eigen, i, value, INSERT_VALUES);
  }  
  Assembly_vector(E_eigen);  // assemble after set value before used

/*
  if(rank==0){
    cout<<"--- view eigen energy ---"<<endl;
    VecView(E_eigen, PETSC_VIEWER_STDOUT_WORLD);
    cout<<"--- view dipole matrix ---"<<endl;
    MatView(H_dipole, PETSC_VIEWER_STDOUT_WORLD);
  }
  return 0;
*/

// observe initial energy and state
  double E_init = Obs_Energy(psi, E_eigen);
  if(rank==0) { 
    cout<<"initial energy = "<< E_init <<endl;  // E = <psi| H |psi>
    //VecView(psi, PETSC_VIEWER_STDOUT_WORLD);
  }

// Laser field
  if(rank==0) cout<<"--- Laser pulse ---"<<endl;
  One_pulse ( wavelength, intensity, n_cycle, dt, cep, t_wait, time, Efield, Afield, period);
  // SINn_SINm_pulse(n_front, n_back, wavelength, intensity, n_cycle, dt, cep, t_wait, time, Efield, Afield, period);

  ofstream  output_field;
  Open_outFile( output_field, output_folder + "/laser.dat");

// time propagation
  KSP  ksp;
  PC   pc;
  KSPCreate(PETSC_COMM_WORLD, &ksp);

  if(rank==0) cout<<"--- start time propagation ---"<<endl;
  PetscScalar cn_factor = -0.5*I*dt[0];  //  -i*dt/2 , for dt[0]=dt[1]
  PetscScalar Et;
  double Ef;
  int nt = time.size();

  for(int i=0; i<nt; i++){
    
    if(i==nt-1) Ef = Efield[i];
    else Ef = (Efield[i+1]+Efield[i])/2.;  // average value of two adjacent points
    Et = complexd(Ef,0.);
   
    PropagateCN_linear_PETSc (cn_factor, Et, E_eigen, H_dipole, ksp, pc, psi );

    if(i%100==0){
      double norm_t = Obs_Norm(psi);
      if(rank==0) cout<< time[i] <<'\t'<< norm_t <<endl;
    }
    if(i%10==0){
      if(rank==0) output_field<< time[i] <<'\t'<< Efield[i] <<'\t'<< Afield[i] <<endl;
    }
  }

  if(rank==0) cout<<"--- end time propagation ---"<<endl;
  KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);   // view information for KSP

// output final spectrum
   PetscScalar  *E_out, *cn2_out;
   vector<double> array_E(n_basis), array_cn2(n_basis);

   if(rank==0)  cout<<"--- output energy ---"<<endl;
   PetscViewer viewer_E;  // vec view by multi processors
   Create_Viewer ( viewer_E, output_folder + "/E_eigen_vec.dat" );
   VecView(E_eigen, viewer_E);

   VecGetArray(E_eigen, &E_out);   // get array 
   for(int i=0; i<n_basis; i++){
      if( nlocal[1] <= i && i < nlocal[2] ){ array_E[i] = real(E_out[i]); }   // store in each processor
   }
   VecRestoreArray(E_eigen, &E_out); //no longer access to array value, MUST restore it to vec, otherwise vec is changed

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

   VecGetArray(vec_cn2, &cn2_out);   // get array
   for(int i=0; i<n_basis; i++){
      if( nlocal[1] <= i && i < nlocal[2] ){ array_cn2[i] = real(cn2_out[i]); }  // store in each processor
   }
   VecRestoreArray(vec_cn2, &cn2_out); //no longer access to array value, MUST restore it to vec, otherwise vec is changed


// finalize
  KSPDestroy(ksp);
  MatDestroy(H_dipole);
  VecDestroy(psi);  
  VecDestroy(E_eigen);  

  PetscFinalize();

  if(rank==0) cout << "---program end---" << endl;

}  // end of main

