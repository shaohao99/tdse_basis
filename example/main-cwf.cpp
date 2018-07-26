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
  

  double k_max, dk;
  int n_theta;
  double r_max, dr;
  int n0, num_basis, l_max;
  string output_folder, cn_folder, wf_folder, dme_folder;

  cout<<"--- input data --- "<<endl;
  ifstream input_data;
  Open_inFile(input_data, "inp-cwf");
  input_data >> num_basis >> l_max >> n0;
  input_data >> r_max >> dr;
  input_data >> k_max >> dk;
  input_data >> n_theta;
  input_data >> wf_folder;
  input_data >> dme_folder;
  input_data >> cn_folder;
  input_data >> output_folder;

  const char* dir_name = output_folder.c_str();
  int dir_exist = mkdir (dir_name, S_IRWXU);  // create output folder
   cout<< "number of basis, maximum angular momentum, maximum number of real atomic states: "<<num_basis <<'\t'<< l_max <<'\t'<< n0 <<endl;
   cout<< "r_max, dr of new unifrom grids: " << r_max <<'\t'<< dr <<endl;
   cout<< "k_max, dk of radial momentum grids: " << k_max <<'\t'<< dk <<endl;
   cout<< "created input folder for wf: " << wf_folder <<endl;
   cout<< "created input folder for dme and wf_new: " << dme_folder <<endl;
   cout<< "created input folder for expasion coefficients: " << cn_folder <<endl;
   cout<< "created output folder: " << output_folder<<endl;

  cout<<"--- input basis info --- "<<endl;
  ifstream input_obt;
  Open_inFile(input_obt, wf_folder + "/inp-orbit");
  vector<int> l_basis(num_basis), n_basis(num_basis), nwf(num_basis);
  vector<double> E_obt(num_basis), r_max_obt(num_basis);
  for(int i=0; i<num_basis; i++){  // iteration of orbital
     input_obt >> n_basis[i] >> l_basis[i] >> r_max_obt[i] >> nwf[i] >> E_obt[i];
  }

  cout<<"--- input expansion coefficients  --- "<<endl;
  vector< vector<double> > c_nl(num_basis,2);
  ifstream input_coeff;
  Open_inFile(input_coeff, cn_folder + "/psi_vec.dat");
  for(int i=0; i<num_basis; i++){  // iteration of orbital
     input_coeff >> c_nl[i][0] >> c_nl[i][1];
     //cout << c_nl[i][0] <<'\t'<< c_nl[i][1] <<endl;
  }

  cout<<"--- input basis functions (uniform grids)  --- "<<endl;
  int nr=int(r_max/dr);
  int n_wf = nr*num_basis;
  vector<double> wf_new(n_wf);
  ifstream input_wf;
  Open_inFile(input_wf, dme_folder + "/wf_new.dat");
  //ofstream output_wf;
  //Open_outFile(output_wf, output_folder+"/wf_check.dat");
  for(int i=0; i<num_basis; i++)  // iteration of orbital
    for(int ir=0; ir<nr; ir++){  // iteration of grids
     int itotal = index2(i,ir,nr);
     input_wf >> wf_new[itotal];
     // if(wf_new[i]!=0.) cout << wf_new[itotal]<<"  ";
     //if(i==0 || i==1 || i==10 || i==100 || i== 600 || i==5000) output_wf << wf_new[itotal]<<endl;
  }

  cout<<"--- project to cwf  --- "<<endl;
  vector<double> k, theta, P_mom, dP_dk_dkrho;
  int nk=int(k_max/dk);
  Mom2D_cwf (k_max, dk, n_theta, wf_new, r_max, dr, n0, num_basis, l_max, l_basis, n_basis, c_nl, k, theta, P_mom, dP_dk_dkrho, output_folder);
     ofstream output_k;
     Open_outFile(output_k, output_folder+"/k_end.dat");
     ofstream output_theta;
     Open_outFile(output_theta, output_folder+"/thetak_end.dat");
     ofstream output_dP_dk;
     Open_outFile(output_dP_dk, output_folder+"/P_mom_end.dat");
     ofstream output_dP_dk_drho;
     Open_outFile(output_dP_dk_drho, output_folder+"/dP_dk_drho_end.dat");
     for(int i=0; i<nk; i++)  output_k << k[i] <<endl;  // in a.u.
     for(int i=0; i<n_theta; i++)  output_theta << theta[i] <<endl;  // in a.u.
     for(int i=0; i<nk*n_theta; i++)  output_dP_dk << P_mom[i] <<endl;  // dP_dk
     for(int i=0; i<nk*n_theta; i++)  output_dP_dk_drho << dP_dk_dkrho[i] <<endl;  // dP_dk

  cout<<"--- end of program  --- "<<endl;

}  // end of main

