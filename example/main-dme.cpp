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
  double r_max_new, dr;
  string wf_folder, dme_folder;

// input data
  if(rank==0) cout<<"--- start input data --- "<<endl;
  ifstream input_data;
  Open_inFile(input_data, "inp-dme");
  input_data >> n_basis;
  input_data >> r_max_new >> dr;
  input_data >> wf_folder;
  input_data >> dme_folder;

  if(rank==0){
   cout<< "number of basis: "<<n_basis <<endl;
   cout<< "r_max, dr of new unifrom grids: " << r_max_new <<'\t'<< dr <<endl;
   cout<< "created input folder for wf: " << wf_folder <<endl;
   cout<< "created input folder fpr dme: " << dme_folder <<endl;
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
    if(rank==0) cout<<"--- finish calculating dipole moments ---"<<endl;
  wf.resize(0);  r.resize(0);    // release memory occupied by old basis wavefunctions

  //return 0;

  MatDestroy(H_dipole);
  VecDestroy(psi);  
  VecDestroy(E_eigen);  

  PetscFinalize();

  if(rank==0) cout << "---program end---" << endl;

}  // end of main

