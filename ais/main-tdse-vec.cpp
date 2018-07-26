// #include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <sys/stat.h>
 
#include "BasisPropagate.h"
#include "BasisLaser.h"
#include "BasisConstant.h"

using namespace std;
using namespace Basis_Propagate;
using namespace Basis_Laser;


int main( int argc , char *argv[] )
{
  
// **************  Initialize Petsc and MPI  *************  //
  static char help[] = "Parallel propagation for AIS of He.\n\n";
  PetscInitialize(&argc, &argv, (char *)0, help);   // MPI_Initialize is also called
  PetscMPIInt    rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

//  declare
  int n_isolate,n_band;
  double Eg,Ea,Eband,dE;
  string output_folder, wf_folder, dme_folder;

  //vector<double> dt(2);
  //double wavelength1, intensity1, n_cycle1, cep1, period1;
  //double intensity2, n_cycle2, cep2, delay, t_wait=0.;
  int nenv, nhhg, iend_pulses;
  double lambda1, int1, n_c1, cep1_pi, int2, n_c2, cep2_pi, ohhg0, delay, dt, t_wait;

// input data
  if(rank==0) cout<<"--- start input data --- "<<endl;
  ifstream input_data;
  Open_inFile(input_data, "inp-tdse");
  input_data >> dt;
  input_data >> lambda1>> int1>> n_c1>> cep1_pi >> nenv;
  input_data >> int2 >> n_c2 >> cep2_pi >> ohhg0 >> nhhg;
  input_data >> delay >> t_wait;
  input_data >> n_isolate>> Eg>> Ea;
  input_data >> n_band>> Eband>> dE;  // eV
  input_data >> wf_folder;
  input_data >> dme_folder;
  input_data >> output_folder;

  vector<double> ohhg(nhhg), ahhg(nhhg), cephhg_pi(nhhg), time, E1, E2, Efield, Afield;

  Eg=Eg/tran_au_eV;  // a.u.
  Ea=Ea/tran_au_eV;
  dE=dE/tran_au_eV; 
  Eband=Eband/tran_au_eV;

  int nE=int(Eband/dE);
  int n_contin=nE*n_band;
  int n_basis=n_isolate + n_contin;
  int n_band_mid=int(n_band/2)+1;
  int nE_mid=int(nE/2)+1;
  double omegaL=45.5896/lambda1;   // IR frequency in a.u.
  vector<double> E_eigen_array(n_basis);

  const char* dir_name = output_folder.c_str();
  int dir_exist = mkdir (dir_name, S_IRWXU);  // create output folder
  if(rank==0){
   cout<<"=== "<<nE<<' '<<n_contin<<' '<<n_isolate<<' '<<n_band<<' '<<nE_mid<<endl; 
   cout<< "number of basis: "<<n_basis <<endl;
   cout<< "created input folder for wf: " << wf_folder <<endl;
   cout<< "created input folder for dme: " << dme_folder <<endl;
   cout<< "created output folder: " << output_folder<<endl;
  }

/*
// test pointer
  if(rank==0) cout<<"--- test pointer ---"<<endl;
  PetscScalar *test;
  PetscMalloc(nE*sizeof(PetscScalar),&test);
    for(int i=0; i<nE; i++){
       test[i]=1.;  
    }
  PetscFree(test);
*/

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
  if(rank==0) cout<<"--- set initio psi ---"<<endl;
  VecSetValue(psi, 0, 1., INSERT_VALUES);  
  Assembly_vector(psi);  // assemble after set value before used

// set eigen energy 
    if(rank==0) cout<<"--- setting E_eigen ---"<<endl;
  E_eigen_array[0]=Eg;
  E_eigen_array[1]=Ea;
  VecSetValue(E_eigen, 0, Eg, INSERT_VALUES);
  VecSetValue(E_eigen, 1, Ea, INSERT_VALUES);
  for(int j=0; j<n_band; j++){
    for(int i=0; i<nE; i++){
     int index=n_isolate+j*nE+i;
     E_eigen_array[index]=Ea+(j-n_band_mid)*omegaL+(i-nE_mid)*dE; 
     //PetscScalar  value = E_eigen_array[index];    // set eigen energy
     VecSetValue(E_eigen, index, E_eigen_array[index], INSERT_VALUES);
     //if(rank==0) cout <<"E_eigen="<<E_eigen_array[index]<<endl;
    }  
  }  
  Assembly_vector(E_eigen);  // assemble after set value before used
    if(rank==0) cout<<"--- end E_eigen ---"<<endl;



// create Mat for Hamiltonian
  if(rank==0) cout<<"--- Create mat of Hamiltonian ---"<<endl;
  Mat Ham;
  Create_Mat ( Ham, nlocal, n_basis, n_basis);  // associated with vec psi by nlocal
  MatSetUp(Ham); // setup before using MatSetValues

// set DME by hand
// to set nonzero structure of a matrix (using INSERT_VALUES) consumes time, it is important to do it before time evolution.
  if(rank==0) cout<<"--- setting dme index ---"<<endl;

  PetscInt *i1, *i2, *j2, *j3, *j4, *j5, *j6, *j7, *i3, *i4, *i5, *i6;
  PetscMalloc(1*sizeof(PetscInt),&i1);
  PetscMalloc(1*sizeof(PetscInt),&i2);
  PetscMalloc(1*sizeof(PetscInt),&j2);
  PetscMalloc(nE*sizeof(PetscInt),&j3);
  PetscMalloc(nE*sizeof(PetscInt),&j4);
  PetscMalloc(nE*sizeof(PetscInt),&j5);
  PetscMalloc(nE*sizeof(PetscInt),&j6);
  PetscMalloc(nE*sizeof(PetscInt),&j7);
  PetscMalloc(nE*sizeof(PetscInt),&i3);
  PetscMalloc(nE*sizeof(PetscInt),&i4);
  PetscMalloc(nE*sizeof(PetscInt),&i5);
  PetscMalloc(nE*sizeof(PetscInt),&i6);

  DME_index_AIS1(nE, n_isolate, i1, i2, j2, j3, j4, j5, j6, j7, i3, i4, i5, i6);  // define nonzero index
  ofstream  output_str;
  Open_outFile( output_str, output_folder + "/index.dat");
  if(rank==0){
    output_str<<"i1: "<<*i1<<endl;
    output_str<<"i2: "<<*i2<<endl;
    output_str<<"i3: "<<*i3<<endl;
    output_str<<"i4: "<<*i4<<endl;
    output_str<<"i5: "<<*i5<<endl;
    output_str<<"i6: "<<*i6<<endl;
    output_str<<"j2: "<<*j2<<endl;
    output_str<<"j3: "<<*j3<<endl;
    output_str<<"j4: "<<*j4<<endl;
    output_str<<"j5: "<<*j5<<endl;
    output_str<<"j6: "<<*j6<<endl;
    output_str<<"j7: "<<*j7<<endl;
  }


// ----------- set dipole moments, energies and nonzero structures of Ham
  if(rank==0) cout<<"--- calculate and store DME in arrays ---"<<endl;
  vector<complexd> ai1j2(1), ai1j3(nE), ai1j5(nE), ai1j7(nE), ai2j4(nE), ai2j6(nE), ai2j5(nE), ai3j4(nE*nE), ai4j5(nE*nE), ai5j6(nE*nE), ai6j7(nE*nE);
  DME_values_AIS1(nE, ai1j2, ai1j3, ai1j5, ai1j7, ai2j4, ai2j6, ai2j5, ai3j4, ai4j5, ai5j6, ai6j7);
  
  if(rank==0) cout<<"--- set initial values and nonzero structure of Ham ---"<<endl;
  DME_array_init_AIS1(nE, E_eigen_array, i1, i2, j2, j3, j4, j5, j6, j7, i3, i4, i5, i6, ai1j2, ai1j3, ai1j5, ai1j7, ai2j4, ai2j6,ai2j5, ai3j4, ai4j5, ai5j6, ai6j7, Ham);
  //MatSetOption(Ham,MAT_HERMITIAN,PETSC_TRUE);  // before assembly begin, set Hermitian only for this SetValues
  //MatSetOption(Ham,MAT_HERMITIAN_ETERNAL,PETSC_TRUE);  // before assembly begin, set Hermitian for all uses below
  Assembly_matrix_final ( Ham ); // finally assemble after setting vales before using
  MatSetOption(Ham,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE);  // after assembly end, keep the same nonzero structure in the following

/*  PetscViewer viewmat;  // Mat view by multi processors
  Create_Viewer ( viewmat, output_folder + "/H_init.dat" );
  MatView(Ham, viewmat);
  if(rank==0) cout<<"--- end view H_init ---"<<endl;*/

// observe initial energy and state
  double E_init = Obs_Energy(psi, E_eigen);
  if(rank==0) { 
    cout<<"Initial energy = "<< E_init <<endl;  // E = <psi| H |psi>
  }

// Laser field
  if(rank==0) cout<<"--- Laser pulse ---"<<endl;
  //One_pulse ( wavelength1, intensity1, n_cycle1, dt, cep1, t_wait, time, Efield, Afield, period1);
   for(int j=0; j<nhhg; j++){
     ohhg[j]=ohhg0+2.*double(j); // odd
     ahhg[j]=1./sqrt(double(nhhg));
     cephhg_pi[j]=cep2_pi;
   }
  IR_APT ( lambda1, int1, n_c1, cep1_pi, nenv, int2, n_c2, nhhg, ohhg, ahhg,  cephhg_pi, delay, dt, t_wait, time, E1, E2, Efield, Afield, iend_pulses);

  ofstream  output_field;
  Open_outFile( output_field, output_folder + "/laser.dat");

  PetscViewer viewmattemp;  // Mat view by multi processors
  Create_Viewer ( viewmattemp, output_folder + "/H_temp.dat" );

// =========== time propagation ===============
  KSP  ksp;
  PC   pc;
  KSPCreate(PETSC_COMM_WORLD, &ksp);

  PetscScalar cn_factor = -0.5*I*dt;  //  -i*dt/2 , for dt[0]=dt[1]
  PetscScalar EXt,ELt;
  double Ef1,Ef2;
  int nt = time.size();

  if(rank==0) cout<<"--- start time propagation ---"<<endl;
  for(int i=0; i<nt; i++){  
    
    if(i==nt-1) { 
      Ef1 = E1[i];  
      Ef2 = E2[i];  
    }
    else {
      Ef1 = (E1[i+1]+E1[i])/2.;  // average value of two adjacent points
      Ef2 = (E2[i+1]+E2[i])/2.;  // average value of two adjacent points
    }
    EXt = complexd(Ef2,0.);
    ELt = complexd(Ef1,0.);
   
    if(rank==0) cout<<"&&& time="<<time[i]<<endl;
    //if(rank==0) cout<<"### update Ham ###"<<endl;

    DME_array_time_AIS1(EXt, ELt, nE, i1, i2, j2, j3, j4, j5, j6, j7, i3, i4, i5, i6, ai1j2, ai1j3, ai1j5, ai1j7, ai2j4, ai2j6, ai3j4, ai4j5, ai5j6, ai6j7, Ham);  // update values of part of Ham in time, not changing the nz structure
    //MatSetOption(Ham,MAT_HERMITIAN,PETSC_TRUE);  // before assembly begin, set Hermitian
    Assembly_matrix_final(Ham);  // assemble after set values before use

    /*if(i%2000==0){
      MatView(Ham, viewmattemp);
    }*/

    //if(rank==0) cout<<"### propagate ###"<<endl;
    Propagate_APT_AIS1 ( cn_factor, Ham, ksp, pc, psi );  // propagate 1 time step

    //if(rank==0) cout<<"### output in time ###"<<endl;
    /*if(i%100==0){
      double norm_t = Obs_Norm(psi);
      if(rank==0) cout<<"time="<<' '<< time[i] <<' '<< norm_t <<endl;
    }*/
    if(i%2==0){
      if(rank==0) output_field<<time[i]<<' '<<Ef1<<' '<<Ef2<<' '<<Efield[i]<<' '<<Afield[i]<<endl;
    }

  } // end time propagation

  if(rank==0) cout<<"--- end time propagation ---"<<endl;
  KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);   // view information for KSP

// output final spectrum
/*   PetscScalar  *E_out, *cn2_out;
   vector<double> array_E(n_basis), array_cn2(n_basis);
*/

   if(rank==0)  cout<<"--- output energy ---"<<endl;
   PetscViewer viewer_E;  // vec view by multi processors
   Create_Viewer ( viewer_E, output_folder + "/E_eigen_vec.dat" );
   VecView(E_eigen, viewer_E);

/*
   VecGetArray(E_eigen, &E_out);   // get array 
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

/*
   VecGetArray(vec_cn2, &cn2_out);   // get array
   for(int i=0; i<n_basis; i++){
      if( nlocal[1] <= i && i < nlocal[2] ){ array_cn2[i] = real(cn2_out[i]); }  // store in each processor
   }
   VecRestoreArray(vec_cn2, &cn2_out); //no longer access to array value, MUST restore it to vec, otherwise vec is changed
*/

// finalize
  KSPDestroy(&ksp);
  MatDestroy(&Ham);
  VecDestroy(&psi);  
  VecDestroy(&E_eigen);  
  PetscFree(i1);
  PetscFree(i2);
  PetscFree(i3);
  PetscFree(i4);
  PetscFree(i5);
  PetscFree(i6);
  PetscFree(j2);
  PetscFree(j3);
  PetscFree(j4);
  PetscFree(j5);
  PetscFree(j6);
  PetscFree(j7);

  PetscFinalize();

  if(rank==0) cout << "---program end---" << endl;

}  // end of main

