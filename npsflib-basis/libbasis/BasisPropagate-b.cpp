#include "BasisPropagate.h"
#include "BasisConstant.h"

#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_coulomb.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include <opari_omp.h>

#include "petscksp.h"
#include "petscts.h"

using namespace std;
//using namespace PETSc_functions;


namespace Basis_Propagate{

 void Open_outFile( ofstream &filename, string name )
 {
  filename.open(name.c_str());
  filename.setf(ios::showpoint|ios::scientific);
 }

 void Open_inFile( ifstream &filename, string name )
 {
  filename.open(name.c_str());
  filename.setf(ios::showpoint|ios::scientific);
 }

  int index3( int k , int j , int i, int nj, int ni )
  {

    int index=( k*nj + j )*ni + i;
    return index;

  }

  int index2( int a , int b, int nb )
  {
    //return a*nb + b;
    int index=a;
    index*=nb;
    index+=b;
    return index;
  }

  int index_wf( int a, vector<int> na, int b )
  {
    int index=0;
    for(int i=0; i<a; i++){
      index += na[i];
    }
    index += b;
    return index;
  }

 void Spherical_grids ( double R0, double dr, double dtheta, double dphi, vector<double> &r, vector<double> &theta, vector<double> &phi)
 {

   int nr = int(R0/dr);
   int ntheta = int(pi/dtheta)+1;
   int nphi = int(2.*pi/dphi)+1;
   cout<<"number of spacial grids: "<< nr<<'\t'<<ntheta<<'\t'<<nphi<<endl;

   r.resize(nr);
   theta.resize(ntheta);
   phi.resize(nphi);

   for(int i=0; i<nr; i++) r[i] = (i+0.1)*dr;  // avoid r=0
   for(int i=0; i<ntheta; i++) theta[i] = i*dtheta;
   for(int i=0; i<nphi; i++) phi[i] = i*dphi;

 }


 void Basis1e_bound (double r, double theta, double phi, int n, int l_max, double Z, vector<complexd> &psi)
 {

   double x = cos(theta);
   int l_number = l_max+1;
   int m_number = l_number*l_number;
   psi.resize(l_number*m_number);
   int index_m=0;

   for(int l=0; l< l_number; l++){

     double R_nl = gsl_sf_hydrogenicR(n, l, Z, r);  // stable and correct. Z=1 for hygrogen atom
 
     for(int m=-l; m<= l; m++){

       double P_lm = gsl_sf_legendre_sphPlm (l, m, x );
       complexd Y_lm = P_lm * exp(I*double(m)*phi);

       psi[index2(index_m,l,l_number)] = R_nl * Y_lm;
       index_m += 1;
     }
   }

 }

 void Basis1e_continuum (double r, double theta, double phi, double k, int l_max, double Z, vector<complexd> &psi)
 {
 
   double xx = r*k;
   double eta = -1./k;
   int l_number = l_max+1;
   double fc_array[l_number];
   double F_exponent[l_number];
   int cwf_F = gsl_sf_coulomb_wave_F_array (0., l_max, eta, xx, fc_array, F_exponent );

   double x = cos(theta);
   int m_number = l_number*l_number;
   psi.resize(l_number*m_number);
   int index_m=0;

   for(int l=0; l< l_number; l++){

     for(int m=-l; m<= l; m++){
     
       double P_lm = gsl_sf_legendre_sphPlm (l, m, x );
       complexd Y_lm = P_lm * exp(I*double(m)*phi);
       
       psi[index2(index_m,l,l_number)] = fc_array[l] * Y_lm;
       index_m += 1;

     }
   }

 }


 void Basis1e_bound_m0 (double r, double theta, int n, int l_max, double Z, vector<double> &psi)
 {

   double x = cos(theta);
   int l_number = l_max+1;
   psi.resize(l_number);

   for(int l=0; l< l_number; l++){
     double R_nl = gsl_sf_hydrogenicR(n, l, Z, r);  // stable and correct. Z=1 for hygrogen atom
     double Y_lm = gsl_sf_legendre_sphPlm (l, 0, x );
     psi[l] = R_nl * Y_lm;
   }

 }

 void Basis1e_continuum_m0 (double r, double theta, double k, int l_max, double Z, vector<double> &psi)
 {

   double xx = r*k;
   double eta = -1./k;
   int l_number = l_max+1;
   double fc_array[l_number];
   double F_exponent[l_number];
   int cwf_F = gsl_sf_coulomb_wave_F_array (0., l_max, eta, xx, fc_array, F_exponent );

   double x = cos(theta);
   psi.resize(l_number);

   for(int l=0; l< l_number; l++){
       double Y_lm = gsl_sf_legendre_sphPlm (l, 0, x );
       psi[l] = fc_array[l] * Y_lm;
   }

 }

// set k_min = k_max if not using continuum states
 void Basis_struct ( int n_max, double k_min, double k_max, double dk, vector<double> &k, vector<int> &n_basis)
 {
   n_basis.resize(6);
   n_basis[1] = (n_max*(n_max+1))/2;  // number of bound states
   n_basis[3] = int( (k_max - k_min)/dk );  // number of k
   n_basis[4] = n_max;  // maximum l of continuum states
   n_basis[5] = n_basis[4]+1;  // number of l for every k
   n_basis[2] = n_basis[3]*n_basis[5];  // number of continuum states 
   n_basis[0] = n_basis[1] + n_basis[2];  // total number of basis
 
   k.resize(n_basis[3]);
   for(int n=0; n<n_basis[3]; n++){
     k[n]= k_min + n*dk;  // momuntum grids
   }
 }


// dipole matrix elements for linearly polarized laser: only uptriangle nonzero elements are calculated2.*(), H_dipole will be set symmetric after "set value" and "assemble".
 void DME_analytical_linear (vector<double> r, int n_max, vector<int> n_basis, vector<double> k, Vec &diag, Mat &H_dipole, int normalize)
 {

   int nr = r.size(); 
   double R0=r[nr-1];
   double dr = r[2]-r[1];
   double dk = k[2]-k[1];
   double Z=1.;

   int n_obt = n_basis[0];
   int n_bound = n_basis[1];
   int n_k = n_basis[3];
   int nr_obt = nr*n_obt;
   vector<double> R_obt(nr_obt);  // Allocate memory. Atomatically deallocated at the end of this function.

   gsl_set_error_handler_off();  // do not invoke due to gsl errors, and do not report the gsl errors.

   //cout<<"*** start calculating eigen functions and values of bound states ***"<<endl;
// calculate bound radial functions and save in memory
   int index_obt = 0;
   for(int n=1; n<=n_max; n++)
     for(int l=0; l<n; l++){

       PetscScalar  value = -0.5*Z*Z/double(n*n);     // set eigen energy
       VecSetValue(diag, index_obt, value, INSERT_VALUES);
      
       for(int i=0; i<nr; i++){
         int in2 = index2(index_obt, i, nr);
           //cout<<"gsl hydrogenicR: "<<gsl_sf_hydrogenicR(n, l, Z, r[i])<<endl;
           //R_obt[in2] = r[i]*gsl_sf_hydrogenicR(n, l, Z, r[i]); 
         gsl_sf_result  result;
         int status_bound = gsl_sf_hydrogenicR_e ( n, l, Z, r[i], &result );
         if( status_bound == GSL_SUCCESS)  R_obt[in2] = r[i]*result.val;  
         else if( status_bound == GSL_EUNDRFLW)  R_obt[in2] = 0.; 
         else if( status_bound == GSL_EOVRFLW)  R_obt[in2] = 1.e3; 
         else R_obt[in2] = 0.;
           //if(n_k==0)  R_obt[index_obt][i] = r[i]*gsl_sf_hydrogenicR(n, l, Z, r[i]);  // r*R_nl,  state normalized, used if there is no continuum states
           //else R_obt[index_obt][i] = r[i]*gsl_sf_hydrogenicR(n, l, Z, r[i])*pow(double(n),1.5);  // r*R_nl*n^(3/2), energy normalized, consistent with cotinuum states
       }
       index_obt++;
   }

   // if(n_k!=0) cout<<"*** eigen calculating functions and values of continuum states ***"<<endl;
// calcualte energy and number of states for continuum by WKB method applying boundary R0
   int lk_max = n_basis[4];
   int l_number = n_basis[5];
   vector<double> Ek(n_k);
   vector< vector<double> > NOS(n_k,l_number);

   for(int n=0; n<n_k; n++){

     Ek[n] = k[n]*k[n]/2.;   // energy in a.u.
        //cout<<"Ek"<<n<<"="<<Ek[n]<<endl;

     for(int l=0; l<l_number; l++){
        double rt = ( sqrt(Z*Z + double(l*(l+1))*Ek[n]) - Z) / Ek[n];  // position r corresponding to classical turing point as a function of Ek
          //cout<<"rt"<<n<<l<<"="<<rt<<endl;
        double rr=0.;
        NOS[n][l]=0.;
        do{            // calculate number of states at Ek as a function of Ek, for Ek > classical turning point, by WKB method.
          int i=0;
          rr=rt+(i+0.001)*dr;  // avoid rr=0 
          NOS[n][l] += sqrt( 2* (Ek[n] + 2.*Z/rr - double(l*(l+1))/(rr*rr) ) );
          i++;
        }while(rr>R0);
        NOS[n][l] *= dr;
        // cout<<"NOS"<<n<<l<<"="<<NOS[n][l]<<endl;
     }  // end of l

   }  // end of k

// calculate DOS and radial wave functions for continuum, save in memory --> file. Need to be discretized.
   double fc_array[l_number];
   double F_exponent[l_number];
   vector< vector<double> > DOS(n_k,l_number);

   for(int n=0; n<n_k; n++){

     for(int l=0; l<l_number; l++){  // calulate density of states at Ek as a function of Ek
        if(n==0) DOS[n][l] = (NOS[n+1][l] - NOS[n][l]) / (Ek[n+1]-Ek[n]);
        else if(n==n_k-1) DOS[n][l] = (NOS[n][l] - NOS[n-1][l]) / (Ek[n]-Ek[n-1]);
        else  DOS[n][l] = ( (NOS[n+1][l] - NOS[n][l]) / (Ek[n+1]-Ek[n]) + (NOS[n][l] - NOS[n-1][l] ) / (Ek[n]-Ek[n-1]) ) / 2.;
        // cout<<"DOS"<<n<<l<<"="<<DOS[n][l]<<endl;
     }  // end of l

     double eta = -Z/k[n];    // for input of Coulomb functions
     for(int i=0; i<nr; i++){   // calculate Coulomb functions and continuum states
       double xx = r[i]*k[n];
       int status_cwf = gsl_sf_coulomb_wave_F_array (0., lk_max, eta, xx, fc_array, F_exponent );
       for(int l=0; l<l_number; l++){
           //R_obt[index_obt+l][i] = fc_array[l];  // fc_array[l] = k*r*R_kl, corresponding to E=k^2 in unit of Ryd, save in the same array with bound states
         int in2 = index2(index_obt+l, i, nr);
         //if( status_cwf == GSL_SUCCESS) R_obt[in2] = fc_array[l] / k[n] ;  // r*R_kl
         if( status_cwf == GSL_SUCCESS) R_obt[in2] = fc_array[l]; 
         else if( status_cwf == GSL_EUNDRFLW) R_obt[in2] = 0.; 
         else if( status_cwf == GSL_EOVRFLW)  R_obt[in2] = 1.e3; 
         else R_obt[in2] = 0.;
           //R_obt[index_obt+l][i] = fc_array[l] / k[n] *sqrt(DOS[n][l])  ;  // energy normalized
       }
     }

     PetscScalar  value = Ek[n];   // set eigen energy
     for(int l=0; l<l_number; l++){
       VecSetValue(diag, index_obt+l, value, INSERT_VALUES);
     }
 
     index_obt += l_number;  // accumulated from bound states
   }  // end number of k

// calculate matrix elements
// initial data 
   int index_obt1=0;   // initial 0 before iteration of n1, l1

   cout<<"*** start iteration of orbitals ***"<<endl;
   for(int n1=1; n1<=n_max; n1++){   // n1, l1 for bound states
       for(int l1=0; l1<n1; l1++){

         cout<<"---obt1="<<index_obt1<<'\t'<<n1<<'\t'<<l1<<endl;
            double sum1=0.;
       int index_obt2=0;   // initial 0 before iteration of n2, l2
       for(int n2=1; n2<=n_max; n2++){   // 1. n2, l2 for bound states, bound-bound, including up and down triangle
         for(int l2=0; l2<n2; l2++){

            if(index_obt2==index_obt1){  // diagnal
                MatSetValue(H_dipole, index_obt1, index_obt2, 0., INSERT_VALUES);
            }

            if(index_obt2>index_obt1){  // only for up-triangle elements
               // if( n1!=n2 && (l2==l1+1 || l2==l1-1) ){    // off-diagonal, dipole selection rule
               // if( l2==l1 ){  // test orthogonal
              if( l2==l1+1 || l2==l1-1 ){    // off-diagonal, dipole selection rule
                 double fac = sqrt( (2.*double(l1)+1.) * (2.*double(l2)+1.) );
                 double c3j = gsl_sf_coupling_3j(2*l2, 2, 2*l1, 0, 0, 0); // integration of theta, input 2l for the first 3 arguments
                 double fac_l = fac*c3j*c3j;
                 double fac_r=0.; 
                 for(int i=0; i<nr; i++){    // integration of r
                   int in1 = index2(index_obt1, i, nr);
                   int in2 = index2(index_obt2, i, nr);
                   fac_r += R_obt[in2] * r[i] * R_obt[in1] ;
                    // fac_r += R_obt[in2] * R_obt[in1] ;
                 }
                 fac_r *= dr; 

                 PetscScalar offdiag;
                 double DOS_bound = pow(double(n2),3.);
                 if(normalize==1) offdiag = fac_l*fac_r*sqrt(DOS_bound);  // energy nomalized, multiply square root of density of final states n^3, consistent with cotinnum states
                 else offdiag = fac_l*fac_r;  // state normalized

                   double fac_r2 = fac_r*fac_r;
                   sum1 += fac_r2;
                   double f_os = (1./double(n1*n1)-1./double(n2*n2))*norm(offdiag);
                   cout<<"obt2="<<index_obt2<<'\t'<<n2<<'\t'<<l2<<'\t'<<fac_l<<'\t'<<fac_r2<<'\t'<<DOS_bound<<'\t'<<f_os<<endl;
                 MatSetValue(H_dipole, index_obt1, index_obt2, offdiag, INSERT_VALUES);
                 MatSetValue(H_dipole, index_obt2, index_obt1, offdiag, INSERT_VALUES);  // symmetric
              }   // end of election rule
            }  // end of up-triangle elements

            index_obt2++;
         } // end of l2
       }  // end of n2 for bound states
            cout <<"---sum of b-b = " << sum1 <<endl;

       sum1 = 0.;
       for(int n2=0; n2<n_k; n2++){   // 2. n2, l2 for continuum states, bound-free
         for(int l2=0; l2<l_number; l2++){

               // if( l2==l1 ){  // test orthogonal
            if(l2==l1+1 || l2==l1-1){    // off-diagonal, dipole selection rule
               double fac = sqrt( (2.*double(l1)+1.) * (2.*double(l2)+1.) );
               double c3j = gsl_sf_coupling_3j(2*l2, 2, 2*l1, 0, 0, 0); // integration of theta
               double fac_l = fac*c3j*c3j;
               double fac_r=0.; 
               for(int i=0; i<nr; i++){    // integration of r
                   int in1 = index2(index_obt1, i, nr);
                   int in2 = index2(index_obt2, i, nr);
                   fac_r += R_obt[in2] * r[i] * R_obt[in1] ;
                     // fac_r += R_obt[in2] * R_obt[in1] ;
               }
               fac_r *= dr;  // energy nomalized, multiply square root of density of final states

               PetscScalar offdiag;
               if(normalize==1) offdiag = fac_l*fac_r*sqrt(DOS[n2][l2]);  // energy normalized
               else offdiag = fac_l*fac_r;    // state normalized
                 double fac_r2 = fac_r*fac_r;
                 sum1 += fac_r2;
                 double f_os = 2.*( Ek[n2] + 0.5/double(n1*n1))*norm(offdiag);
                 cout<<"obt2="<<index_obt2<<'\t'<<k[n2]<<'\t'<<l2<<'\t'<<fac_l<<'\t'<<fac_r2<<'\t'<<DOS[n2][l2]<<'\t'<<f_os<<endl;
               MatSetValue(H_dipole, index_obt1, index_obt2, offdiag, INSERT_VALUES);
               MatSetValue(H_dipole, index_obt2, index_obt1, offdiag, INSERT_VALUES);  // symmetric
            }   // end of election rule

            index_obt2++;    // accumulated from bound states, hence index_obt2>index_obt1, all are up-triangle
         } // end of l2
       }  // end of n2 for continuum states
            cout <<"---sum of b-f = " << sum1 <<endl;

       index_obt1++;
     }  // end of l1
   }  // end of n1 for bound states


   for(int n1=0; n1<n_k; n1++){   // 3. n1, l1 for continuum states, free-free
     for(int l1=0; l1<l_number; l1++){

       cout<<"---obt1="<<index_obt1<<'\t'<<k[n1]<<'\t'<<l1<<endl;
       int index_obt2=n_bound;   // initial value = n_bound  before iteration of n2, l2 for continuum
          double sum2=0;
       for(int n2=0; n2<n_k; n2++){   // n2, l2 for continuum states, free-free
         for(int l2=0; l2<l_number; l2++){

            if(index_obt2==index_obt1){  // diagnal
                 MatSetValue(H_dipole, index_obt1, index_obt2, 0., INSERT_VALUES);
            }

            if(index_obt2>index_obt1){  // only for up-triangle elements
               // if(l2==l1){  
               if(l2==l1+1 || l2==l1-1){    // off-diagonal, dipole selection rule
                 double fac = sqrt( (2.*double(l1)+1.) * (2.*double(l2)+1.) );
                 double c3j = gsl_sf_coupling_3j(2*l2, 2, 2*l1, 0, 0, 0); // integration of theta
                 double fac_l = fac*c3j*c3j;
                 double fac_r=0.; 
                 for(int i=0; i<nr; i++){    // integration of r
                    int in1 = index2(index_obt1, i, nr);
                    int in2 = index2(index_obt2, i, nr);
                    fac_r += R_obt[in2] * r[i] * R_obt[in1] ;
                     //fac_r += R_obt[in2] * R_obt[in1] ;
                 }
                 fac_r *= dr;  // energy nomalized, multiply square root of density of final states

                 PetscScalar offdiag;
                 if(normalize==1) offdiag = fac_l*fac_r*sqrt(DOS[n2][l2]);  // energy normalized
                 else offdiag = fac_l*fac_r;    // state normalized
                   double fac_r2 = fac_r*fac_r;
                   sum2 += fac_r2;
                   double f_os = 2.*(Ek[n2]-Ek[n1])*norm(offdiag);
                   cout<<"obt2="<<index_obt2<<'\t'<<k[n2]<<'\t'<<l2<<'\t'<<fac_l<<'\t'<<fac_r2<<'\t'<<DOS[n2][l2]<<'\t'<<f_os<<endl;
                 MatSetValue(H_dipole, index_obt1, index_obt2, offdiag, INSERT_VALUES);
                 MatSetValue(H_dipole, index_obt2, index_obt1, offdiag, INSERT_VALUES);  // symmetric
              }   // end of election rule
            }  // end of up-triangle elements

            index_obt2++;   // 
         } // end of l2
       }  // end of n2 for continuum states
            cout <<"---sum of f-f = " << sum2 <<endl;

       index_obt1++;   // accumulated from bound states
     }  // end of l1
   }  // end of n1 for continuum states

 }  // end of matrix elements


// calculate basis wavefunctions, input: basis functions
 void DME_linear (int n_basis, vector<int> l, vector<int> nwf, vector<double> r_max, vector<double> r, vector<double> wf, double r_max_new, double dr, Mat &H_dipole)
 {

   ofstream  output_test1;
   Open_outFile( output_test1, "state_old.dat");
   ofstream  output_test2;
   Open_outFile( output_test2, "state_new.dat");

   for(int i=0; i<n_basis; i++){

     MatSetValue(H_dipole, i, i, 0., INSERT_VALUES);  // diagnal

     int nr1=nwf[i];  // pick up initial wf
     double r1[nr1], wf1[nr1];
     for(int ir = 0; ir<nr1; ir++){
       int itotal = index_wf(i, nwf, ir);
       r1[ir] = r[itotal];
       wf1[ir] = wf[itotal];
       if(i==399) output_test1 <<r1[ir]<<'\t'<<wf1[ir]<<endl;
     }

     cout<<"&&& initial &&&  "<<i<<endl;
     gsl_interp_accel *acc1 = gsl_interp_accel_alloc ();  // interpolation: cubic spline
     gsl_spline *spline1  = gsl_spline_alloc (gsl_interp_cspline, nr1);
     gsl_spline_init (spline1, r1, wf1, nr1);

     //cout<<"&&& cspline 2 &&&"<<endl;
     int n_new = int(r_max_new/dr);   // new orbit
     vector<double> wf1_new(n_new);
     double sum1=0., sum2=0.;
     for (int inew=0; inew<n_new; inew++){  
        double r_new = (inew+0.1) * dr ;
        if(r_new > r_max[i]) wf1_new[inew] = 0.;
        else wf1_new[inew] = gsl_spline_eval (spline1, r_new, acc1);
         sum1 += wf1_new[inew]*wf1_new[inew];
         // sum2 += wf1_new[inew]*wf1_new[inew]*r_new*r_new;
         if(i==399) output_test2 <<r_new<<'\t'<<wf1_new[inew]<<endl;
     }
     cout<<i<<'\t'<<sum1*dr<<'\t'<<sum2*dr<<endl;

     //cout<<"&&& cspline 3 &&&"<<endl;
     gsl_interp_accel_free (acc1);
     gsl_spline_free (spline1);
     //cout<<"&&& cspline 4 &&&"<<endl;

     for(int j=i+1; j<n_basis; j++){   // only for up-triangle elements

         if( l[i]==l[j]+1 || l[i]==l[j]-1 ){    // off-diagonal, dipole selection rule
            cout<<"### final ###  "<<j<<endl;

            double fac = sqrt( (2.*double(l[i])+1.) * (2.*double(l[j])+1.) );  // angular part
            double c3j = gsl_sf_coupling_3j(2*l[j], 2, 2*l[i], 0, 0, 0); // integration of theta, input 2l for the first 3 arguments
            double fac_l = fac*c3j*c3j;
 
            int nr2=nwf[j];  // pick up wf2
            double r2[nr2], wf2[nr2];
            for (int ir = 0; ir<nr2; ir++){
               int itotal = index_wf(j, nwf, ir);
               r2[ir] = r[itotal];
               wf2[ir] = wf[itotal];
                  //cout<<r2[ir]<<'\t'<<wf2[ir]<<endl;
            }

            //cout<<"*** cspline 1 ***"<<endl;
            gsl_interp_accel *acc2 = gsl_interp_accel_alloc ();  // interpolation: cubic spline
            gsl_spline *spline2  = gsl_spline_alloc (gsl_interp_cspline, nr2);
            gsl_spline_init (spline2, r2, wf2, nr2);

            cout<<"*** cspline 2 ***"<<endl;
            double fac_r = 0., wf2_new=0.;
            for (int inew=0; inew<n_new; inew++){  // intergration of r on new grirds
               double r_new = (inew+0.1) * dr ;
               if(r_new > r_max[j]) wf2_new = 0.;
               else wf2_new = gsl_spline_eval (spline2, r_new, acc2);
               fac_r += wf1_new[inew] * r_new * wf2_new ;
               //cout<<r_new<<'\t'<<wf1_new[inew]<<'\t'<<wf2_new<<'\t'<<fac_r<<endl;
            }
            fac_r *= dr;
            cout<<"*** cspline 3 ***"<<endl;

            gsl_interp_accel_free (acc2);
            gsl_spline_free (spline2);
            //cout<<"*** cspline 4 ***"<<endl;

            PetscScalar offdiag;   // set matrix elements
            offdiag = fac_l*fac_r;  
            MatSetValue(H_dipole, i, j, offdiag, INSERT_VALUES);
            MatSetValue(H_dipole, j, i, offdiag, INSERT_VALUES);  // symmetric
         }   // end of election rule

      }  // end of j

    }  // end of i

 }  // end of function


// PETSc functions:
  void Create_Vec (Vec &vec, int nsize)
  {
   VecCreate(PETSC_COMM_WORLD, &vec);
   VecSetSizes(vec, PETSC_DECIDE, nsize);
   VecSetFromOptions(vec);   // used by setting "-vec_type mpi" at run time
  }

  void Get_Local_Info ( Vec vec, vector<int> &nlocal)
  {
      PetscInt rstart, rend;
      VecGetOwnershipRange(vec, &rstart, &rend);

      PetscInt local_number;
      VecGetLocalSize(vec, &local_number);

      nlocal.resize(3);
      nlocal[1] = rstart;   // start index
      nlocal[2] = rend;     // end index
      nlocal[0] = local_number;  // number of local data

        // cout<<"nlocal="<<nlocal[0]<<'\t'<<"rstart="<<nlocal[1]<<'\t'<<"rend="<<nlocal[2]<<endl;
  }


  void Create_Mat ( Mat &matrix, vector<int> nlocal, PetscInt n, PetscInt m )
  {
      MatCreate(PETSC_COMM_WORLD, &matrix);
      PetscInt local_number=nlocal[0];
      MatSetSizes(matrix, local_number, local_number, n, m);
      MatSetFromOptions(matrix);   // used by set "-mat_type mpiaij" at run time, default mpiaij
      // MatSetType(matrix, MATMPIDENSE);  // MATMPIBAIJ, MATMPIDENSE
      // MatSetBlockSize(matrix, wf.n1);
  }


  void Assembly_vector ( Vec &vec )
  {
     VecAssemblyBegin(vec);
     VecAssemblyEnd(vec);
  }

  void Assembly_matrix_final ( Mat &matrix )
  {
     MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY);
     MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY);
  }

  void Assembly_matrix_flush ( Mat &matrix )
  {
     MatAssemblyBegin(matrix, MAT_FLUSH_ASSEMBLY);
     MatAssemblyEnd(matrix, MAT_FLUSH_ASSEMBLY);
  }

  void Create_Viewer (PetscViewer &viewer, string filename )
  {
      // PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);   
      PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
      PetscViewerSetType(viewer, PETSC_VIEWER_ASCII);
      PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_SYMMODU);   // PETSC_VIEWER_ + DEFAULT, ASCII_MATLAB, ASCII_COMMON, ASCII_INDEX, ASCII_SYMMODU

      char name[filename.size()];   // convert string to char
      for(int i=0; i<filename.size(); i++){
         name[i] = filename[i];
      }
        // cout<<name<<endl;
      PetscViewerFileSetName(viewer, name);
  }


// norm of vec
  double Obs_Norm ( Vec psi  )
  {
    PetscScalar  vnorm;
    VecDot (psi, psi, &vnorm);  // vnorm = <psi|psi> 
    double obs = abs(vnorm);
    return obs;
  } // end of Obs_Norm


// sqrt norm of vec
  PetscReal Obs_Norm_sqrt ( Vec psi  )
  {
    PetscReal  vnorm_sqrt;
    VecNorm(psi, NORM_2, &vnorm_sqrt);  // for type NROM_2, wf_norm = sqrt( a1^2 + a2^2 + ... ) is calculated
    return vnorm_sqrt;
  }

//  Normalize of a vec, phi' = phi/sqrt( <phi|phi> )
  void Normalize( Vec &psi )
  {
    PetscReal  norm_fac = 1./Obs_Norm_sqrt(psi);
    VecScale(psi, norm_fac);   
  } // end of Normalize


// calculate the mean value of any physical quantity that can be expressed by a matrix, e.g. the total/kinetic energy, 
  double Obs_Energy ( Vec psi, Vec E_eigen )
  {
     PetscScalar  obs_cplx;
     Vec temp_vec;
     VecDuplicate(psi, &temp_vec);  // get values of wavefucntion of the last step, as well as local info

     VecAbs(psi);  // |c_n|, psi is temporally changed, but not changed outside this function
     VecPointwiseMult(temp_vec, psi, psi);  // |c_n|^2
     VecDot(temp_vec, E_eigen, &obs_cplx);  // Sum_n { |c_n|^2*E_n }

     double obs = real(obs_cplx);
     return obs;
  }


// propagation by Cranck-Nicholson, implemented by PETSc. Define cn_factor=-i*dt/2.
  void PropagateCN_linear_PETSc ( PetscScalar cn_factor, PetscScalar Et, Vec diag, Mat H_dipole, KSP &ksp, PC &pc, Vec &psi )
  {

    PetscScalar  P_real1=complexd(1.,0.);
    PetscScalar  P_realm1=complexd(-1.,0.);
    PetscScalar  P_real2=complexd(2.,0.);

    //cout<< "*** vec duplicate ***"<<endl;
    Vec vec_temp;
    VecDuplicate(psi, &vec_temp);  // get values of wavefucntion of the last step, as well as local info
    Assembly_vector(vec_temp);

    //cout<< "*** Mat duplicate ***"<<endl;
    Mat H_temp;
    MatDuplicate(H_dipole, MAT_COPY_VALUES, &H_temp);  // H_temp = H_dipole, copy values and nonzero sturcture
    // Assembly_matrix_flush(H_temp);   // flushly assemble between setting values
    Assembly_matrix_final(H_temp);   // finally assemble after setting values brgore use
 
    //cout<< "*** Mat scale ***"<<endl;
    MatScale(H_temp, Et); // H_temp = Et*H_dipole
    MatDiagonalSet(H_temp, diag, ADD_VALUES);  // H_temp = Et*H_dipole + diag = H
    Assembly_matrix_final(H_temp);   // finally assemble after the last setting values
      // MatView(H_temp, PETSC_VIEWER_STDOUT_WORLD);

    //cout<< "*** Mat vec calc ***"<<endl;
    MatScale(H_temp, cn_factor); // H_temp = -i*dt/2*H
    MatShift(H_temp, P_real1);   // H_temp = I - i*dt/2*H = H_right
    MatMult(H_temp, psi, vec_temp);   // vec_temp = (I - i*dt/2*H) * psi(t) = b
    MatScale(H_temp, P_realm1);   // H_temp = -H_right 
    MatShift(H_temp, P_real2);   // H_temp = 2I - H_right = I + i*dt/2*H = H_left
      // MatView(H_temp, PETSC_VIEWER_STDOUT_WORLD);
    
    //cout<< "*** KSP setting ***"<<endl;
    KSPSetOperators(ksp, H_temp, H_temp, DIFFERENT_NONZERO_PATTERN);   // Set operators. Here the matrix that defines the linear system also serves as the preconditioning matrix.
    KSPSetType(ksp, KSPGMRES);   // slect the ksp method, difault is GMRES. Good: KSPGMRES > KSPCGS; Bad: KSPBICG, KSPCHEBYCHEV, KSPCG
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCJACOBI);    // Good: PCBJACOBI > PCJACOBI; Bad: PCILU, PCLU, PCASM
    KSPSetTolerances(ksp, 1.e-10, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT); // default: 1.e-7; diverge: 1e-4
    KSPSetFromOptions(ksp);   // Set runtime options, e.g., -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>. These options will override those specified above as long as KSPSetFromOptions() is called _after_ any other customization routines.

    //cout<< "*** KSP solve ***"<<endl;
    KSPSolve(ksp, vec_temp, psi);   // solve H_left * psi(t+dt) = b, and save the solution to the new psi
     // KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD); 

    //cout<< "*** destroy ***"<<endl;
    MatDestroy(H_temp);  // release memory
    VecDestroy(vec_temp);

  }


}  // end of namepsace



