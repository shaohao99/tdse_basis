#include "BasisPropagate.h"
#include "BasisConstant.h"

#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_coulomb.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

//#include <opari_omp.h>

#include "petscts.h"
#include "petscksp.h"
#include "petscvec.h"
#include "petscmat.h"
#include "petscis.h"


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

// transfer an integer to a character through ASCI code
  void int_to_string(int i, string &a)
  {

   int j1=48, j2=48;  // ascii code: 48 --> 0, for ouput file names
   if(i<10) { j1=48; j2=i+48; }
   if(10<=i && i<20) { j1=1+48; j2=i-10+48; }
   if(20<=i && i<30) { j1=2+48; j2=i-20+48; }
   if(30<=i && i<40) { j1=3+48; j2=i-30+48; }
   if(40<=i && i<50) { j1=4+48; j2=i-40+48; }
   if(50<=i && i<60) { j1=5+48; j2=i-50+48; }
 
   a=char(j1); 
   a += char(j2);

  }

// calculate Gamma fucntion by Lanczos approximation: http://en.wikipedia.org/wiki/Lanczos_approximation . good for small and large Re(z)
  inline complexd Gamma_Lanczos (complexd z)
  {

    complexd x,t;
    double g = 7.;
    vector<double> p(9);
    p[0] = 0.99999999999980993;
    p[1] = 676.5203681218851;
    p[2] = -1259.1392167224028;
    p[3] = 771.32342877765313;
    p[4] = -176.61502916214059;
    p[5] = 12.507343278686905;
    p[6] = -0.13857109526572012;
    p[7] = 9.9843695780195716e-6;
    p[8] = 1.5056327351493116e-7;

    if( real(z) < 0.5)  // reflection
        return pi / (sin(pi*z) * Gamma_Lanczos(1.-z));
    else
    {
        z -= 1.;
        x = p[0];
        for( int i=1; i<9; i++) x += p[i]/(z+double(i));
        t = z + g + 0.5;
        return sqrt(2.*pi) * pow(t,z+0.5) * exp(-t) * x;
    }

  }

// factorial function
  // inline int Factorial(int n)
  inline double Factorial(int n)
  {
      return (n<2) ? 1. : n*Factorial(n-1);
  }

// spherical harmonics by numerical recipes. Formula 6.8.2. 
 complexd Y_lm_plgndr ( double theta, double phi, int l, int m)
 {

   double x = cos(theta);
   double plm = plgndr(l, m, x);
     // double mm = m;
   double fac = sqrt( (2.*l+1.)/(4.*pi) * Factorial(l-m)/Factorial(l+m) );
   complexd Y_lm = fac * plm * exp(I*double(m)*phi);
   return Y_lm;
     // cout<<Y_lm<<"\t"<<fac<<"\t"<< plm<<"\t"<< exp(I*m*phi)<<"\t"<<I*m*phi<<"\t"<<m<<"\t"<<phi;

 }

// associated Legendre polynomial by numerical recipes. x=cos(theta)
  double plgndr(int l, int m, double x)
  {
        double fact,pll,pmm,pmmp1,somx2;
        int i,ll;

        if (m < 0 || m > l || fabs(x) > 1.0)
           cout<<"Bad arguments in routine plgndr"<<endl;
           pmm=1.0;
        if (m > 0) {
                somx2=sqrt((1.0-x)*(1.0+x));
                fact=1.0;
                for (i=1;i<=m;i++) {
                        pmm *= -fact*somx2;
                        fact += 2.0;
                }
        }
        if (l == m)
                return pmm;
        else {
                pmmp1=x*(2*m+1)*pmm;
                if (l == (m+1))
                        return pmmp1;
                else {
                        for (ll=m+2;ll<=l;ll++) {
                                pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
                                pmm=pmmp1;
                                pmmp1=pll;
                        }
                        return pll;
                }
        }
  }


// for analytical basis
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

/*
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
*/

// calculate basis wavefunctions, input: basis functions
// void DME_linear (int n_basis, vector<int> l, vector<int> nwf, vector<double> r_max, vector<double> &r, vector<double> &wf, double r_max_new, double dr, Mat &H_dipole, string folder)
 void DME_linear (int n_basis, vector<int> l, vector<int> nwf, vector<double> r_max, vector<double> &r, vector<double> &wf, double r_max_new, double dr, string folder)
 {

   ofstream  output_test1;
   Open_outFile( output_test1, folder + "/state_old.dat");
   ofstream  output_test2;
   Open_outFile( output_test2, folder + "/state_new.dat");
   ofstream  output_wf_new;
   Open_outFile( output_wf_new, folder + "/wf_new.dat");
   ofstream  output_dme;
   Open_outFile( output_dme, folder + "/dme.dat");

   int n_new = int(r_max_new/dr);   // new grids
   vector<double> wf_new(n_basis*n_new);

   for(int i=0; i<n_basis; i++){  // from old to new grids by spline

     cout<<"--- i="<<i<<endl;
     int ni=nwf[i];  
     double ri[ni], wf_i[ni];
     for(int iold = 0; iold<ni; iold++){
       int itotal = index_wf(i, nwf, iold);
       ri[iold] = r[itotal];
       wf_i[iold] = wf[itotal];
       if(i==6139) output_test1 <<ri[iold]<<'\t'<<wf_i[iold]<<endl;
     }
     cout<<"### read old wf ###"<<endl;

     //cout<<"&&& basis &&&  "<<i<<endl;
     gsl_interp_accel *acc = gsl_interp_accel_alloc ();  // interpolation: cubic spline
     gsl_spline *spline  = gsl_spline_alloc (gsl_interp_cspline, ni);
     gsl_spline_init (spline, ri, wf_i, ni);
     cout<<"### gsl inerpolate ###"<<endl;

     // double sum1=0.;
     for (int inew=0; inew<n_new; inew++){
        double r_new = (inew+0.1) * dr ;
        int itotal_new = index2(i,inew,n_new);
        if(r_new > r_max[i]) wf_new[itotal_new] = 0.;
        else wf_new[itotal_new] = gsl_spline_eval (spline, r_new, acc);
          // sum1 += wf_new[itotal_new]*wf_new[itotal_new];
        if(i==6139) output_test2 <<r_new<<'\t'<<wf_new[itotal_new]<<endl;
        output_wf_new <<wf_new[itotal_new]<<" ";
     }
     // cout<<i<<'\t'<<sum1*dr<<endl;
     cout<<"### spline for new wf ###"<<endl;

     gsl_interp_accel_free (acc);
     gsl_spline_free (spline);
     //cout<<"&&& enf of cspline &&&"<<endl;
  }  // end of spline

  r.resize(0);  wf.resize(0);  // release memory occupied by old basis wavefunctions

   cout<<"%%% start dme %%%"<<endl;
   for(int i=0; i<n_basis; i++){

     // MatSetValue(H_dipole, i, i, 0., INSERT_VALUES);  // diagnal

     for(int j=i+1; j<n_basis; j++){   // only for up-triangle elements

         if(l[i]==l[j]){    // test orthogonal
            double orthog = 0.;
            for (int inew=0; inew<n_new; inew++){ 
               int itotal_new = index2(i,inew,n_new);
               int jtotal_new = index2(j,inew,n_new);
               orthog += wf_new[itotal_new] * wf_new[jtotal_new];  // test orthogonal
            }
            orthog *= dr;
            cout <<"---     "<< i <<'\t'<< j <<'\t'<< orthog <<endl;
         }

         if( l[i]==l[j]+1 || l[i]==l[j]-1 ){    // off-diagonal, dipole selection rule

             // cout<<i<<'\t'<<j<<endl;
            double fac = sqrt( (2.*double(l[i])+1.) * (2.*double(l[j])+1.) );  // angular part
            double c3j = gsl_sf_coupling_3j(2*l[j], 2, 2*l[i], 0, 0, 0); // integration of theta, input 2l for the first 3 arguments
            double fac_l = fac*c3j*c3j;
 
            double fac_r = 0.;
            // int inew;
//#pragma omp parallel shared (fac_r,dr)
//{
//     #pragma omp for private(inew)
            for (int inew=0; inew<n_new; inew++){  // intergration of r on new grirds
                // int ip=omp_get_thread_num();
               double r_new = (inew+0.1) * dr ;
               int itotal_new = index2(i,inew,n_new);
               int jtotal_new = index2(j,inew,n_new);
               fac_r += wf_new[itotal_new] * r_new * wf_new[jtotal_new];  // dipole moment
            }
// }  // end of parallel
            fac_r *= dr;

            double dme = fac_l*fac_r;  
            output_dme<<dme<<" ";
            /*PetscScalar offdiag;   // set matrix elements
            offdiag = dme;  
            MatSetValue(H_dipole, i, j, offdiag, INSERT_VALUES);
            MatSetValue(H_dipole, j, i, offdiag, INSERT_VALUES);  // symmetric */
        }   // end of election rule

      }  // end of j
    }  // end of i
   cout<<"%%% end dme %%%"<<endl;

  wf_new.resize(0);  // release memory occupied by new basis wavefunctions

 }  // end of function

 void Set_DME (int n_basis, vector<int> l, Mat &H_dipole, string folder)
 {

   ifstream  input_dme;
   Open_inFile( input_dme, folder + "/dme.dat");

   for(int i=0; i<n_basis; i++){

     MatSetValue(H_dipole, i, i, 0., INSERT_VALUES);  // diagnal

     for(int j=i+1; j<n_basis; j++){   // only for up-triangle elements

         if( l[i]==l[j]+1 || l[i]==l[j]-1 ){    // off-diagonal, dipole selection rule

            double dme;
            input_dme>>dme;
            PetscScalar offdiag;   // set matrix elements
            offdiag = dme;
            MatSetValue(H_dipole, i, j, offdiag, INSERT_VALUES);
            MatSetValue(H_dipole, j, i, offdiag, INSERT_VALUES);  // symmetric
        }   // end of election rule

      }  // end of j
    }  // end of i

 }

// wrong !!
 void Set_DME_select (int n_basis, vector<int> l, vector<int> n, vector<int> nb, Mat &H_dipole, string folder)
 {

   ifstream  input_dme;
   Open_inFile( input_dme, folder + "/dme.dat");
   int i_select=0;

   for(int i=0; i<n_basis; i++){

     int j_select = i_select + 1;  // index of array, wrong!!!

     if(n[i]<=nb[0] || n[i]>nb[1]){  // select bound states
       MatSetValue(H_dipole, i_select, i_select, 0., INSERT_VALUES);  // diagnal
       i_select++;  // index of row
     }

     for(int j=i+1; j<n_basis; j++){   // only for up-triangle elements

         if( l[i]==l[j]+1 || l[i]==l[j]-1 ){    // off-diagonal, dipole selection rule

            double dme;
            input_dme>>dme;

            if(n[i]<=nb[0] || n[i]>nb[1]){  // select bound states
              PetscScalar offdiag;   // set matrix elements
              offdiag = dme;
              MatSetValue(H_dipole, i_select, j_select, offdiag, INSERT_VALUES);
              MatSetValue(H_dipole, j_select, i_select, offdiag, INSERT_VALUES);  // symmetric
              j_select++;
            }
        }   // end of election rule

      }  // end of j
    }  // end of i

 }


// elliptical polarized: Salpeter book, p254, Eq(60.11)
 void DME_elliptical (int n_nl, int &n_basis, vector<int> l, vector<int> nwf, vector<double> r_max, vector<double> &r, vector<double> &wf, double r_max_new, double dr, string folder)
 {

   ofstream  output_test1;
   Open_outFile( output_test1, folder + "/state_old.dat");
   ofstream  output_test2;
   Open_outFile( output_test2, folder + "/state_new.dat");
   ofstream  output_dme;
   Open_outFile( output_dme, folder + "/dme.dat");

   int n_new = int(r_max_new/dr);   // new grids
   vector<double> wf_new(n_nl*n_new);

   for(int i=0; i<n_nl; i++){  // from old to new grids by spline

     int ni=nwf[i];
     double ri[ni], wf_i[ni];
     for(int iold = 0; iold<ni; iold++){
       int itotal = index_wf(i, nwf, iold);
       ri[iold] = r[itotal];
       wf_i[iold] = wf[itotal];
       if(i==6139) output_test1 <<ri[iold]<<'\t'<<wf_i[iold]<<endl;
     }

     //cout<<"&&& basis &&&  "<<i<<endl;
     gsl_interp_accel *acc = gsl_interp_accel_alloc ();  // interpolation: cubic spline
     gsl_spline *spline  = gsl_spline_alloc (gsl_interp_cspline, ni);
     gsl_spline_init (spline, ri, wf_i, ni);

     // double sum1=0.;
     for (int inew=0; inew<n_new; inew++){
        double r_new = (inew+0.1) * dr ;
        int itotal_new = index2(i,inew,n_new);
        if(r_new > r_max[i]) wf_new[itotal_new] = 0.;
        else wf_new[itotal_new] = gsl_spline_eval (spline, r_new, acc);
          // sum1 += wf_new[itotal_new]*wf_new[itotal_new];
        if(i==6139) output_test2 <<r_new<<'\t'<<wf_new[itotal_new]<<endl;
     }
     // cout<<i<<'\t'<<sum1*dr<<endl;

     gsl_interp_accel_free (acc);
     gsl_spline_free (spline);
     //cout<<"&&& enf of cspline &&&"<<endl;
  }  // end of spline

  r.resize(0);  wf.resize(0);  // release memory occupied by old basis wavefunctions

   for(int i=0; i<n_nl; i++)
      for(int mi=-l[i]; mi<=l[i]; mi++ ){  // iteration of m
         n_basis++;  // count number of basis
   }

   for(int i=0; i<n_nl; i++){
     for(int j=i+1; j<n_nl; j++){   // only for up-triangle elements

         if( l[i]==l[j]+1 || l[i]==l[j]-1 ){    // off-diagonal, dipole selection rule

            // radial
            double fac_r = 0.;
            int inew;
            for (inew=0; inew<n_new; inew++){  // intergration of r on new grirds
                // int ip=omp_get_thread_num();
               double r_new = (inew+0.1) * dr ;
               int itotal_new = index2(i,inew,n_new);
               int jtotal_new = index2(j,inew,n_new);
               fac_r += wf_new[itotal_new] * r_new * wf_new[jtotal_new];  // dipole moment
            }
            fac_r *= dr;

            // expand subspace of m
            double x_lm=0., xy_lm=0.;
            complexd y_lm=complexd(0.,0.);
            for(int mi=-l[i]; mi<=l[i]; mi++ ){  // iteration of m
               for(int mj=-l[i]; mj<=l[i]; mj++ ) {  // iteration of m'

                 // angula
                 if(mi==mj-1 && l[i]==l[j]-1 ){  // m'=m+1, l'=l+1
                    xy_lm = sqrt( (l[i]+mi+2)*(l[i]+mi+1)/(2*l[i]+3)/(2*l[i]+1) );
                    y_lm = -I*xy_lm/2.;
                    x_lm = xy_lm/2.;
                    double dme_x = x_lm*fac_r;
                    complexd dme_y = y_lm*fac_r;
                    output_dme<<dme_x<<" "<<dme_y<<endl;
                 }
                 else if(mi==mj+1 && l[i]==l[j]-1 ){  // m'=m-1, l'=l+1
                    xy_lm = -sqrt( (l[i]-mi+2)*(l[i]-mi+1)/(2*l[i]+3)/(2*l[i]+1) );
                    y_lm = I*xy_lm/2.;
                    x_lm = xy_lm/2.;
                    double dme_x = x_lm*fac_r;
                    complexd dme_y = y_lm*fac_r;
                    output_dme<<dme_x<<" "<<dme_y<<endl;
                 }
                 else if(mi==mj-1 && l[i]==l[j]+1 ){  // m'=m+1, l'=l-1
                    xy_lm = -sqrt( (l[i]-mi)*(l[i]-mi-1)/(2*l[i]+1)/(2*l[i]-1) );
                    y_lm = -I*xy_lm/2.;
                    x_lm = xy_lm/2.;
                    double dme_x = x_lm*fac_r;
                    complexd dme_y = y_lm*fac_r;
                    output_dme<<dme_x<<" "<<dme_y<<endl;
                 }
                 else if(mi==mj+1 && l[i]==l[j]+1 ){  // m'=m-1, l'=l-1
                    xy_lm = -sqrt( (l[i]+mi)*(l[i]+mi-1)/(2*l[i]+1)/(2*l[i]-1) );
                    y_lm = I*xy_lm/2.;
                    x_lm = xy_lm/2.;
                    double dme_x = x_lm*fac_r;
                    complexd dme_y = y_lm*fac_r;
                    output_dme<<dme_x<<" "<<dme_y<<endl;
                 }

               }  // end of m'
            }  // end of m

        }   // end of election rule

      }  // end of j
    }  // end of i

  wf_new.resize(0);  // release memory occupied by new basis wavefunctions

 }  // end of function

/*
// set dme for circular pulse
 void Set_DME_circular (int n_basis, vector<int> l, Mat &H_dipole_x,  Mat &H_dipole_y string folder)
 {

   ifstream  input_dme;
   Open_inFile( input_dme, folder + "/dme.dat");

   for(int i=0; i<n_basis; i++){

     MatSetValue(H_dipole, i, i, 0., INSERT_VALUES);  // diagnal

     for(int j=i+1; j<n_basis; j++){   // only for up-triangle elements

         if( l[i]==l[j]+1 || l[i]==l[j]-1 ){    // off-diagonal, dipole selection rule

            int indexm=0; 
            for( int mi=-l[i]; mi<=l[i]; mi++ ){ 
               indexm++;   // count index of m subspace

               for( int mj=-l[i]; mj<=l[i]; mj++ ) {  // iteration of m

                 double dme_x, dme_y;
                 PetscScalar offdiag_x[], offdiag_y[];   // set matrix elements

                 if(mi==mj-1 && l[i]==l[j]-1 ){  // m'=m+1, l'=l+1
                   input_dme>>dme_x>>dme_y;
                   offdiag_x[indexm] = dme_x;
                   offdiag_y[indexm] = dme_y;
                 }

                 MatSetValue(H_dipole, i, j, offdiag, INSERT_VALUES);
                 MatSetValue(H_dipole, j, i, offdiag, INSERT_VALUES);  // symmetric

               }  // end of m'
            }  // end of m

        }   // end of election rule

      }  // end of j
    }  // end of i

 }
*/

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

  void Assembly_matrix_final ( Mat &matrix )  // use MAT_FINAL_ASSEMBLY for the final assembly before using the matrix.
  {
     MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY);
     MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY);
  }

  void Assembly_matrix_flush ( Mat &matrix ) // Use MAT_FLUSH_ASSEMBLY when switching between ADD_VALUES and INSERT_VALUES in MatSetValues()
  {
     MatAssemblyBegin(matrix, MAT_FLUSH_ASSEMBLY);
     MatAssemblyEnd(matrix, MAT_FLUSH_ASSEMBLY);
  }

  void Create_Viewer (PetscViewer &viewer, string filename )
  {
      // PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);   
      PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
      //PetscViewerSetType(viewer, PETSC_VIEWER_ASCII);
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
     VecDuplicate(psi, &temp_vec);  // duplicate local info

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
    VecDuplicate(psi, &vec_temp);  // duplicate local info
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

/*
// DME offdiagonal nonzero structure, before time iteration
  void DME_str_AIS1(int nE, int n_isolate, vector<int> &i1, vector<int> &i2,  vector<int>&j2, vector<int> &j3, vector<int> &j4, vector<int> &j5, vector<int> &j6, vector<int> &j7, vector<int> &i3, vector<int> &i4, vector<int> &i5, vector<int> &i6){  // globel index

    //i1.resize(1);

    i1[0]=0; 
    i2[0]=1; 
    j2[0]=1; 
    for(int j=0; j<nE; j++){ // column
       j3[j]=n_isolate+j;
       j4[j]=n_isolate+nE+j;
       j5[j]=n_isolate+2*nE+j;
       j6[j]=n_isolate+3*nE+j;
       j7[j]=n_isolate+4*nE+j;
    }
    for(int i=0; i<nE; i++){ // row
       i3[i]=n_isolate+i;
       i4[i]=n_isolate+nE+i;
       i5[i]=n_isolate+2*nE+i;
       i6[i]=n_isolate+3*nE+i;
    }

  }
*/

// DME offdiagonal nonzero structure, before time iteration
  void DME_str_AIS1(int nE, int n_isolate, PetscInt *i1, PetscInt *i2, PetscInt *j2, PetscInt *j3, PetscInt *j4, PetscInt *j5, PetscInt *j6, PetscInt *j7, PetscInt *i3, PetscInt *i4, PetscInt *i5, PetscInt *i6){  // globel index

    *i1=0; 
    *i2=1; 
    *j2=1; 
    for(int j=0; j<nE; j++){ // column
       *j3=n_isolate+j;
       *j4=n_isolate+nE+j;
       *j5=n_isolate+2*nE+j;
       *j6=n_isolate+3*nE+j;
       *j7=n_isolate+4*nE+j;
       j3++;
       j4++;
       j5++;
       j6++;
       j7++;
    }
    for(int i=0; i<nE; i++){ // row
       *i3=n_isolate+i;
       *i4=n_isolate+nE+i;
       *i5=n_isolate+2*nE+i;
       *i6=n_isolate+3*nE+i;
       i3++;
       i4++;
       i5++;
       i6++;
    }

  }  // end str


// set DME values, before time iteration
  void DME_init_AIS1(int nE, vector<double> E_eigen_array, PetscInt *i1, PetscInt *i2, PetscInt *j2, PetscInt *j3, PetscInt *j4, PetscInt *j5, PetscInt *j6, PetscInt *j7, PetscInt *i3, PetscInt *i4, PetscInt *i5, PetscInt *i6, Mat &H_temp){

    PetscScalar *vi1j2, *vi1j3, *vi1j5, *vi1j7, *vi2j4, *vi2j6, *vi2j5, *vi3j4, *vi4j5, *vi5j6, *vi6j7;
    PetscScalar *vj2i1, *vj3i1, *vj5i1, *vj7i1, *vj4i2, *vj6i2, *vj5i2, *vj4i3, *vj5i4, *vj6i5, *vj7i6;
    PetscMalloc(1*sizeof(PetscScalar),&vi1j2);
    PetscMalloc(1*sizeof(PetscScalar),&vj2i1);
    PetscMalloc(nE*sizeof(PetscScalar),&vi1j3);
    PetscMalloc(nE*sizeof(PetscScalar),&vj3i1);
    PetscMalloc(nE*sizeof(PetscScalar),&vi1j5);
    PetscMalloc(nE*sizeof(PetscScalar),&vj5i1);
    PetscMalloc(nE*sizeof(PetscScalar),&vi1j7);
    PetscMalloc(nE*sizeof(PetscScalar),&vj7i1);
    PetscMalloc(nE*sizeof(PetscScalar),&vi2j4);
    PetscMalloc(nE*sizeof(PetscScalar),&vj4i2);
    PetscMalloc(nE*sizeof(PetscScalar),&vi2j6);
    PetscMalloc(nE*sizeof(PetscScalar),&vj6i2);
    PetscMalloc(nE*sizeof(PetscScalar),&vi2j5);
    PetscMalloc(nE*sizeof(PetscScalar),&vj5i2);
    PetscMalloc(nE*nE*sizeof(PetscScalar),&vi3j4);
    PetscMalloc(nE*nE*sizeof(PetscScalar),&vj4i3);
    PetscMalloc(nE*nE*sizeof(PetscScalar),&vi4j5);
    PetscMalloc(nE*nE*sizeof(PetscScalar),&vj5i4);
    PetscMalloc(nE*nE*sizeof(PetscScalar),&vi5j6);
    PetscMalloc(nE*nE*sizeof(PetscScalar),&vj6i5);
    PetscMalloc(nE*nE*sizeof(PetscScalar),&vi6j7);
    PetscMalloc(nE*nE*sizeof(PetscScalar),&vj7i6);
    double factor=1.;  //1000.;

    *vi1j2=2.*factor;  // g,a, by XUV
    *vj2i1=conj(*vi1j2);
    for(int i=0; i<nE; i++){
       *vi1j3=1.*factor;  // g,cotinuum, by XUV
       *vj3i1=conj(*vi1j3);
       vi1j3++;
       vj3i1++;

       *vi1j5=1.*factor;
       *vj5i1=conj(*vi1j5);
       vi1j5++;
       vj5i1++;

       *vi1j7=1.*factor;
       *vj7i1=conj(*vi1j7);
       vi1j7++;
       vj7i1++;

       *vi2j4=4.*factor;  // a,cotinuum, by IR
       *vj4i2=conj(*vi2j4);
       vi2j4++;
       vj4i2++;

       *vi2j6=4.*factor;
       *vj6i2=conj(*vi2j6);
       vi2j6++;
       vj6i2++;

       *vi2j5=3.*factor;  // CI
       *vj5i2=conj(*vi2j5);
       vi2j5++;
       vj5i2++;

       for(int j=0; j<nE; j++){
         *vi3j4=5.*factor;  // continuum,cotinuum, by IR
         *vj4i3=conj(*vi3j4);
         vi3j4++;
         vj4i3++;

         *vi4j5=5.*factor; 
         *vj5i4=conj(*vi4j5);
         vi4j5++;
         vj5i4++;

         *vi5j6=5.*factor; 
         *vj6i5=conj(*vi5j6);
         vi5j6++;
         vj6i5++;

         *vi6j7=5.*factor; 
         *vj7i6=conj(*vi6j7);
         vi6j7++;
         vj7i6++;

       }
    }

// offdiagonal: up triangle
    MatSetValues(H_temp,1,i1,1,j2,vi1j2,INSERT_VALUES);
    MatSetValues(H_temp,1,i1,nE,j3,vi1j3,INSERT_VALUES);
    MatSetValues(H_temp,1,i1,nE,j5,vi1j5,INSERT_VALUES);
    MatSetValues(H_temp,1,i1,nE,j7,vi1j7,INSERT_VALUES);
    MatSetValues(H_temp,1,i2,nE,j4,vi2j4,INSERT_VALUES);
    MatSetValues(H_temp,1,i2,nE,j6,vi2j6,INSERT_VALUES);
    MatSetValues(H_temp,1,i2,nE,j5,vi2j5,INSERT_VALUES);
    MatSetValues(H_temp,nE,i3,nE,j4,vi3j4,INSERT_VALUES);
    MatSetValues(H_temp,nE,i4,nE,j5,vi4j5,INSERT_VALUES);
    MatSetValues(H_temp,nE,i5,nE,j6,vi5j6,INSERT_VALUES);
    MatSetValues(H_temp,nE,i6,nE,j7,vi6j7,INSERT_VALUES);
// offdiagonal: down triangle, Hermitian
    MatSetValues(H_temp,1,j2,1,i1,vj2i1,INSERT_VALUES);
    MatSetValues(H_temp,nE,j3,1,i1,vj3i1,INSERT_VALUES);
    MatSetValues(H_temp,nE,j5,1,i1,vj5i1,INSERT_VALUES);
    MatSetValues(H_temp,nE,j7,1,i1,vj7i1,INSERT_VALUES);
    MatSetValues(H_temp,nE,j4,1,i2,vj4i2,INSERT_VALUES);
    MatSetValues(H_temp,nE,j6,1,i2,vj6i2,INSERT_VALUES);
    MatSetValues(H_temp,nE,j5,1,i2,vj5i2,INSERT_VALUES);
    MatSetValues(H_temp,nE,j4,nE,i3,vj4i3,INSERT_VALUES);
    MatSetValues(H_temp,nE,j5,nE,i4,vj5i4,INSERT_VALUES);
    MatSetValues(H_temp,nE,j6,nE,i5,vj6i5,INSERT_VALUES);
    MatSetValues(H_temp,nE,j7,nE,i6,vj7i6,INSERT_VALUES);

// diagonal
    int n_b=E_eigen_array.size();
    for(int i=0; i<n_b; i++){
       MatSetValue(H_temp, i, i, E_eigen_array[i], INSERT_VALUES);  // diagnal = eigen energies
    }

  }  // end DME_static


// IS generate: index in parallel ???
  void DME_IS_AIS1(int nE, PetscInt *i1, PetscInt *i2, PetscInt *j2, PetscInt *j3, PetscInt *j4, PetscInt *j5, PetscInt *j6, PetscInt *j7, PetscInt *i3, PetscInt *i4, PetscInt *i5, PetscInt *i6, IS &isi1, IS &isi2, IS &isi3, IS &isi4, IS &isi5, IS &isi6, IS &isj2, IS &isj3, IS &isj4, IS &isj5, IS &isj6, IS &isj7){
//  void DME_IS_AIS1(int nE, PetscInt *i1, IS *isi1){

    /*ISCreateGeneral(PETSC_COMM_WORLD, 1, i1, PETSC_COPY_VALUES, &isi1);  //rows
    ISCreateGeneral(PETSC_COMM_WORLD, 1, i2, PETSC_COPY_VALUES, &isi2);
    ISCreateGeneral(PETSC_COMM_WORLD, nE, i3, PETSC_COPY_VALUES, &isi3);
    ISCreateGeneral(PETSC_COMM_WORLD, nE, i4, PETSC_COPY_VALUES, &isi4);
    ISCreateGeneral(PETSC_COMM_WORLD, nE, i5, PETSC_COPY_VALUES, &isi5);
    ISCreateGeneral(PETSC_COMM_WORLD, nE, i6, PETSC_COPY_VALUES, &isi6);
    ISCreateGeneral(PETSC_COMM_WORLD, 1, j2, PETSC_COPY_VALUES,  &isj2); // columns
    ISCreateGeneral(PETSC_COMM_WORLD, nE, j3, PETSC_COPY_VALUES, &isj3);
    ISCreateGeneral(PETSC_COMM_WORLD, nE, j4, PETSC_COPY_VALUES, &isj4);
    ISCreateGeneral(PETSC_COMM_WORLD, nE, j5, PETSC_COPY_VALUES, &isj5);
    ISCreateGeneral(PETSC_COMM_WORLD, nE, j6, PETSC_COPY_VALUES, &isj6);
    ISCreateGeneral(PETSC_COMM_WORLD, nE, j7, PETSC_COPY_VALUES, &isj7);*/

    ISCreateGeneral(PETSC_COMM_WORLD, 1, i1, &isi1);  //rows
    ISCreateGeneral(PETSC_COMM_WORLD, 1, i2, &isi2);
    ISCreateGeneral(PETSC_COMM_WORLD, nE, i3, &isi3);
    ISCreateGeneral(PETSC_COMM_WORLD, nE, i4, &isi4);
    ISCreateGeneral(PETSC_COMM_WORLD, nE, i5, &isi5);
    ISCreateGeneral(PETSC_COMM_WORLD, nE, i6, &isi6);
    ISCreateGeneral(PETSC_COMM_WORLD, 1, j2,  &isj2); // columns
    ISCreateGeneral(PETSC_COMM_WORLD, nE, j3, &isj3);
    ISCreateGeneral(PETSC_COMM_WORLD, nE, j4, &isj4);
    ISCreateGeneral(PETSC_COMM_WORLD, nE, j5, &isj5);
    ISCreateGeneral(PETSC_COMM_WORLD, nE, j6, &isj6);
    ISCreateGeneral(PETSC_COMM_WORLD, nE, j7, &isj7);

  }   // end of globel index

// get submatrix form H_static
  void DME_block_AIS1(Mat Ham, IS isi1, IS isi2, IS isi3, IS isi4, IS isi5, IS isi6, IS isj2, IS isj3, IS isj4, IS isj5, IS isj6, IS isj7, Mat &mat_i1j2, Mat &mat_i1j3, Mat &mat_i1j5, Mat &mat_i1j7, Mat &mat_i2j4, Mat &mat_i2j6, Mat &mat_i3j4, Mat &mat_i4j5, Mat &mat_i5j6, Mat &mat_i6j7, Mat &mat_j2i1, Mat &mat_j3i1, Mat &mat_j5i1, Mat &mat_j7i1, Mat &mat_j4i2, Mat &mat_j6i2, Mat &mat_j4i3, Mat &mat_j5i4, Mat &mat_j6i5, Mat &mat_j7i6){


    //cout<<"$$$ Get submatrix $$$"<<endl;
    MatGetSubMatrix(Ham,isi1,isj2,1,MAT_INITIAL_MATRIX,&mat_i1j2);  // up triangle
    MatGetSubMatrix(Ham,isi1,isj3,1,MAT_INITIAL_MATRIX,&mat_i1j3);
    MatGetSubMatrix(Ham,isi1,isj5,1,MAT_INITIAL_MATRIX,&mat_i1j5);
    MatGetSubMatrix(Ham,isi1,isj7,1,MAT_INITIAL_MATRIX,&mat_i1j7);
    MatGetSubMatrix(Ham,isi2,isj4,1,MAT_INITIAL_MATRIX,&mat_i2j4);
    MatGetSubMatrix(Ham,isi2,isj6,1,MAT_INITIAL_MATRIX,&mat_i2j6);
    MatGetSubMatrix(Ham,isi3,isj4,1,MAT_INITIAL_MATRIX,&mat_i3j4);
    MatGetSubMatrix(Ham,isi4,isj5,1,MAT_INITIAL_MATRIX,&mat_i4j5);
    MatGetSubMatrix(Ham,isi5,isj6,1,MAT_INITIAL_MATRIX,&mat_i5j6);
    MatGetSubMatrix(Ham,isi6,isj5,1,MAT_INITIAL_MATRIX,&mat_i6j7);
    MatGetSubMatrix(Ham,isj2,isi1,1,MAT_INITIAL_MATRIX,&mat_j2i1);  // low triangle
    MatGetSubMatrix(Ham,isj3,isi1,1,MAT_INITIAL_MATRIX,&mat_j3i1);
    MatGetSubMatrix(Ham,isj5,isi1,1,MAT_INITIAL_MATRIX,&mat_j5i1);
    MatGetSubMatrix(Ham,isj7,isi1,1,MAT_INITIAL_MATRIX,&mat_j7i1);
    MatGetSubMatrix(Ham,isj4,isi2,1,MAT_INITIAL_MATRIX,&mat_j4i2);
    MatGetSubMatrix(Ham,isj6,isi2,1,MAT_INITIAL_MATRIX,&mat_j6i2);
    MatGetSubMatrix(Ham,isj4,isi3,1,MAT_INITIAL_MATRIX,&mat_j4i3);
    MatGetSubMatrix(Ham,isj5,isi4,1,MAT_INITIAL_MATRIX,&mat_j5i4);
    MatGetSubMatrix(Ham,isj6,isi5,1,MAT_INITIAL_MATRIX,&mat_j6i5);
    MatGetSubMatrix(Ham,isj7,isi6,1,MAT_INITIAL_MATRIX,&mat_j7i6);

    /*MatGetSubMatrix(Ham,isi1,isj2,MAT_INITIAL_MATRIX,&mat_i1j2);  // up triangle
    MatGetSubMatrix(Ham,isi1,isj3,MAT_INITIAL_MATRIX,&mat_i1j3);
    MatGetSubMatrix(Ham,isi1,isj5,MAT_INITIAL_MATRIX,&mat_i1j5);
    MatGetSubMatrix(Ham,isi1,isj7,MAT_INITIAL_MATRIX,&mat_i1j7);
    MatGetSubMatrix(Ham,isi2,isj4,MAT_INITIAL_MATRIX,&mat_i2j4);
    MatGetSubMatrix(Ham,isi2,isj6,MAT_INITIAL_MATRIX,&mat_i2j6);
    MatGetSubMatrix(Ham,isi3,isj4,MAT_INITIAL_MATRIX,&mat_i3j4);
    MatGetSubMatrix(Ham,isi4,isj5,MAT_INITIAL_MATRIX,&mat_i4j5);
    MatGetSubMatrix(Ham,isi5,isj6,MAT_INITIAL_MATRIX,&mat_i5j6);
    MatGetSubMatrix(Ham,isi6,isj5,MAT_INITIAL_MATRIX,&mat_i6j7);
    MatGetSubMatrix(Ham,isj2,isi1,MAT_INITIAL_MATRIX,&mat_j2i1);  // low triangle
    MatGetSubMatrix(Ham,isj3,isi1,MAT_INITIAL_MATRIX,&mat_j3i1);
    MatGetSubMatrix(Ham,isj5,isi1,MAT_INITIAL_MATRIX,&mat_j5i1);
    MatGetSubMatrix(Ham,isj7,isi1,MAT_INITIAL_MATRIX,&mat_j7i1);
    MatGetSubMatrix(Ham,isj4,isi2,MAT_INITIAL_MATRIX,&mat_j4i2);
    MatGetSubMatrix(Ham,isj6,isi2,MAT_INITIAL_MATRIX,&mat_j6i2);
    MatGetSubMatrix(Ham,isj4,isi3,MAT_INITIAL_MATRIX,&mat_j4i3);
    MatGetSubMatrix(Ham,isj5,isi4,MAT_INITIAL_MATRIX,&mat_j5i4);
    MatGetSubMatrix(Ham,isj6,isi5,MAT_INITIAL_MATRIX,&mat_j6i5);
    MatGetSubMatrix(Ham,isj7,isi6,MAT_INITIAL_MATRIX,&mat_j7i6);*/

}


// set offdiagonal DME values, during time iteration
//  void DME_time_AIS1(Mat &Ham, PetscScalar EXt, PetscScalar ELt, IS isi1, IS isi2, IS isi3, IS isi4, IS isi5, IS isi6, IS isj2, IS isj3, IS isj4, IS isj5, IS isj6, IS isj7, Mat mat_i1j2, Mat mat_i1j3, Mat mat_i1j5, Mat mat_i1j7, Mat mat_i2j4, Mat mat_i2j6, Mat mat_i3j4, Mat mat_i4j5, Mat mat_i5j6, Mat mat_i6j7, Mat mat_j2i1, Mat mat_j3i1, Mat mat_j5i1, Mat mat_j7i1, Mat mat_j4i2, Mat mat_j6i2, Mat mat_j4i3, Mat mat_j5i4, Mat mat_j6i5, Mat mat_j7i6){
  void DME_time_AIS1(Mat &H_temp, PetscScalar EXt, PetscScalar ELt, int nE, PetscInt *i1, PetscInt *i2, PetscInt *j2, PetscInt *j3, PetscInt *j4, PetscInt *j5, PetscInt *j6, PetscInt *j7, PetscInt *i3, PetscInt *i4, PetscInt *i5, PetscInt *i6, Mat mat_i1j2, Mat mat_i1j3, Mat mat_i1j5, Mat mat_i1j7, Mat mat_i2j4, Mat mat_i2j6, Mat mat_i3j4, Mat mat_i4j5, Mat mat_i5j6, Mat mat_i6j7, Mat mat_j2i1, Mat mat_j3i1, Mat mat_j5i1, Mat mat_j7i1, Mat mat_j4i2, Mat mat_j6i2, Mat mat_j4i3, Mat mat_j5i4, Mat mat_j6i5, Mat mat_j7i6){

    //cout<<"$$$ Mat scale $$$"<<endl;
    MatScale(mat_i1j2,EXt);  // up triangle
    MatScale(mat_i1j3,EXt);
    MatScale(mat_i1j5,EXt);
    MatScale(mat_i1j7,EXt);
    MatScale(mat_i2j4,ELt);
    MatScale(mat_i2j6,ELt);
    MatScale(mat_i3j4,EXt);
    MatScale(mat_i4j5,EXt);
    MatScale(mat_i5j6,EXt);
    MatScale(mat_i6j7,EXt);
    MatScale(mat_j2i1,EXt);  // low triangle
    MatScale(mat_j3i1,EXt);
    MatScale(mat_j5i1,EXt);
    MatScale(mat_j7i1,EXt);
    MatScale(mat_j4i2,ELt);
    MatScale(mat_j6i2,ELt);
    MatScale(mat_j4i3,ELt);
    MatScale(mat_j5i4,ELt);
    MatScale(mat_j6i5,ELt);
    MatScale(mat_j7i6,ELt);
/*
// ouput mat for test
    PetscViewer viewtime;  // Mat view by multi processors
    Create_Viewer ( viewtime, "H_time.dat" );
    MatView(mat_i1j2, viewtime);
    MatView(mat_i1j3, viewtime);
    MatView(mat_i2j4, viewtime);
    MatView(mat_i4j5, viewtime);
    MatView(mat_j5i4, viewtime);
*/

// Get values from mat to array
    PetscScalar *vi1j2, *vi1j3, *vi1j5, *vi1j7, *vi2j4, *vi2j6, *vi3j4, *vi4j5, *vi5j6, *vi6j7;
    PetscScalar *vj2i1, *vj3i1, *vj5i1, *vj7i1, *vj4i2, *vj6i2, *vj4i3, *vj5i4, *vj6i5, *vj7i6;

    MatGetArray(mat_i1j2,&vi1j2);  // up triangle
    MatGetArray(mat_i1j3,&vi1j3);
    MatGetArray(mat_i1j5,&vi1j5);
    MatGetArray(mat_i1j7,&vi1j7);
    MatGetArray(mat_i2j4,&vi2j4);
    MatGetArray(mat_i2j6,&vi2j6);
    MatGetArray(mat_i3j4,&vi3j4);
    MatGetArray(mat_i4j5,&vi4j5);
    MatGetArray(mat_i5j6,&vi5j6);
    MatGetArray(mat_i6j7,&vi6j7);
    MatGetArray(mat_j2i1,&vj2i1);  // down triangle
    MatGetArray(mat_j3i1,&vj3i1); 
    MatGetArray(mat_j5i1,&vj5i1); 
    MatGetArray(mat_j7i1,&vj7i1); 
    MatGetArray(mat_j4i2,&vj4i2); 
    MatGetArray(mat_j6i2,&vj6i2); 
    MatGetArray(mat_j4i3,&vj4i3); 
    MatGetArray(mat_j5i4,&vj5i4); 
    MatGetArray(mat_j6i5,&vj6i5); 
    MatGetArray(mat_j7i6,&vj7i6); 

// Mat set values
// offdiagonal: up triangle
    MatSetValues(H_temp,1,i1,1,j2,vi1j2,INSERT_VALUES);
    MatSetValues(H_temp,1,i1,nE,j3,vi1j3,INSERT_VALUES);
    MatSetValues(H_temp,1,i1,nE,j5,vi1j5,INSERT_VALUES);
    MatSetValues(H_temp,1,i1,nE,j7,vi1j7,INSERT_VALUES);
    MatSetValues(H_temp,1,i2,nE,j4,vi2j4,INSERT_VALUES);
    MatSetValues(H_temp,1,i2,nE,j6,vi2j6,INSERT_VALUES);
    MatSetValues(H_temp,nE,i3,nE,j4,vi3j4,INSERT_VALUES);
    MatSetValues(H_temp,nE,i4,nE,j5,vi4j5,INSERT_VALUES);
    MatSetValues(H_temp,nE,i5,nE,j6,vi5j6,INSERT_VALUES);
    MatSetValues(H_temp,nE,i6,nE,j7,vi6j7,INSERT_VALUES);
// offdiagonal: down triangle, Hermitian
    MatSetValues(H_temp,1,j2,1,i1,vj2i1,INSERT_VALUES);
    MatSetValues(H_temp,nE,j3,1,i1,vj3i1,INSERT_VALUES);
    MatSetValues(H_temp,nE,j5,1,i1,vj5i1,INSERT_VALUES);
    MatSetValues(H_temp,nE,j7,1,i1,vj7i1,INSERT_VALUES);
    MatSetValues(H_temp,nE,j4,1,i2,vj4i2,INSERT_VALUES);
    MatSetValues(H_temp,nE,j6,1,i2,vj6i2,INSERT_VALUES);
    MatSetValues(H_temp,nE,j4,nE,i3,vj4i3,INSERT_VALUES);
    MatSetValues(H_temp,nE,j5,nE,i4,vj5i4,INSERT_VALUES);
    MatSetValues(H_temp,nE,j6,nE,i5,vj6i5,INSERT_VALUES);
    MatSetValues(H_temp,nE,j7,nE,i6,vj7i6,INSERT_VALUES);

// sub matrix restore
    MatRestoreArray(mat_i1j2,&vi1j2);  // up triangle
    MatRestoreArray(mat_i1j3,&vi1j3);
    MatRestoreArray(mat_i1j5,&vi1j5);
    MatRestoreArray(mat_i1j7,&vi1j7);
    MatRestoreArray(mat_i2j4,&vi2j4);
    MatRestoreArray(mat_i2j6,&vi2j6);
    MatRestoreArray(mat_i3j4,&vi3j4);
    MatRestoreArray(mat_i4j5,&vi4j5);
    MatRestoreArray(mat_i5j6,&vi5j6);
    MatRestoreArray(mat_i6j7,&vi6j7);
    MatRestoreArray(mat_j2i1,&vj2i1);  // down triangle
    MatRestoreArray(mat_j3i1,&vj3i1);
    MatRestoreArray(mat_j5i1,&vj5i1);
    MatRestoreArray(mat_j7i1,&vj7i1);
    MatRestoreArray(mat_j4i2,&vj4i2);
    MatRestoreArray(mat_j6i2,&vj6i2);
    MatRestoreArray(mat_j4i3,&vj4i3);
    MatRestoreArray(mat_j5i4,&vj5i4);
    MatRestoreArray(mat_j6i5,&vj6i5);
    MatRestoreArray(mat_j7i6,&vj7i6);

/*
    cout<<"$$$ Submatrix update $$$"<<endl;
    MatSubMatrixUpdate(mat_i1j2,Ham,isi1,isj2);  // up triangle
    MatSubMatrixUpdate(mat_i1j3,Ham,isi1,isj3);
    MatSubMatrixUpdate(mat_i1j5,Ham,isi1,isj5);
    MatSubMatrixUpdate(mat_i1j7,Ham,isi1,isj7);
    MatSubMatrixUpdate(mat_i2j4,Ham,isi2,isj4);
    MatSubMatrixUpdate(mat_i2j6,Ham,isi2,isj6);
    MatSubMatrixUpdate(mat_i3j4,Ham,isi3,isj4);
    MatSubMatrixUpdate(mat_i4j5,Ham,isi4,isj5);
    MatSubMatrixUpdate(mat_i5j6,Ham,isi5,isj6);
    MatSubMatrixUpdate(mat_i6j7,Ham,isi6,isj7);
    MatSubMatrixUpdate(mat_j2i1,Ham,isj2,isi1);  // low triangle
    MatSubMatrixUpdate(mat_j3i1,Ham,isj3,isi1);
    MatSubMatrixUpdate(mat_j5i1,Ham,isj5,isi1);
    MatSubMatrixUpdate(mat_j7i1,Ham,isj7,isi1);
    MatSubMatrixUpdate(mat_j4i2,Ham,isj4,isi2);
    MatSubMatrixUpdate(mat_j6i2,Ham,isj6,isi2);
    MatSubMatrixUpdate(mat_j4i3,Ham,isj4,isi3);
    MatSubMatrixUpdate(mat_j5i4,Ham,isj5,isi4);
    MatSubMatrixUpdate(mat_j6i5,Ham,isj6,isi5);
    MatSubMatrixUpdate(mat_j7i6,Ham,isj7,isi6);
*/
  }  // submatrix update


// propagation by Cranck-Nicholson of APT and AIS1, implemented by PETSc. Define cn_factor=-i*dt/2.
  void Propagate_APT_AIS1 ( PetscScalar cn_factor, Mat Ham, KSP &ksp, PC &pc, Vec &psi )
  {

    PetscScalar  P_real1=complexd(1.,0.);
    PetscScalar  P_realm1=complexd(-1.,0.);
    PetscScalar  P_real2=complexd(2.,0.);

    //cout<< "*** vec duplicate ***"<<endl;
    Vec vec_temp;
    VecDuplicate(psi, &vec_temp);  // duplicate local info
    Assembly_vector(vec_temp);

   //cout<< "*** Mat duplicate ***"<<endl;
    Mat H_temp;
    MatDuplicate(Ham, MAT_COPY_VALUES, &H_temp);  // H_temp = H_static, copy values and nonzero sturcture
    // Assembly_matrix_flush(H_temp);   // flushly assemble between setting values
    Assembly_matrix_final(H_temp);   // finally assemble after setting values brgore use

    //cout<< "*** Mat vec calc ***"<<endl;
    MatScale(H_temp, cn_factor); // H = -i*dt/2*H
    MatShift(H_temp, P_real1);   // H = I - i*dt/2*H = H_right
    MatMult(H_temp, psi, vec_temp);   // vec_temp = (I - i*dt/2*H) * psi(t) = b
    MatScale(H_temp, P_realm1);   // H = -H_right 
    MatShift(H_temp, P_real2);   // H = 2I - H_right = I + i*dt/2*H = H_left
      // MatView(H, PETSC_VIEWER_STDOUT_WORLD);
    
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

    VecScale(psi,2.);

    //cout<< "*** destroy ***"<<endl;
    MatDestroy(H_temp);  // release memory
    VecDestroy(vec_temp);

  } // end APT_AIS1


// project to Coulomb wave
  void Mom2D_cwf ( double k_max, double dk, int n_theta, vector<double> wf_new, double r_max, double dr, int n0, int num_basis, int l_max, vector<int> l_basis, vector<int> n_basis, vector< vector<double> > c_nl, vector<double> &k, vector<double> &theta, vector<double> &P_mom, vector<double> &dP_dk_dkrho, string output_folder)
 {

    ofstream output_pl;
    Open_outFile(output_pl, output_folder+"/Pl_mom_end.dat");

    int nk=int(k_max/dk); 
    int nr=int(r_max/dr); 
    double dtheta=pi/(n_theta-1);
    int l_num = l_max+1;
    vector<complexd> fac_phase(l_num);
    vector<complexd> Ak(l_num);

    k.resize(nk);  
    theta.resize(n_theta); 
    P_mom.resize(nk*n_theta);
    dP_dk_dkrho.resize(nk*n_theta);

    for(int i=0;i<nk;i++) {k[i] = (i+0.5)*dk; }   // k_Ryd[i] = k[i]/sqrt(2.); }
    for(int j=0;j<n_theta;j++) theta[j] = j*dtheta; 

    double fac = sqrt(2./pi);
    //vector< vector<complexd> > Y_lm(n_theta,l_num);
    vector< vector<complexd> > Y_lm;
    Y_lm.resize(n_theta);
    for(int i=0; i<n_theta; i++) Y_lm[i].resize(l_num);

    for(int l=0; l<l_num; l++)  
      for(int j=0;j<n_theta;j++)  // calculate spherical function
      {
         Y_lm[j][l] =  Y_lm_plgndr (theta[j], 0., l, 0);  // m=0
      }

// test output
   ofstream  output_cwf;
   Open_outFile( output_cwf, "cwf_test.dat");

// big iteration
    for(int ik=0; ik<nk; ik++){  // interation of k

      double eta = -1./k[ik];
      for(int l=0; l<l_num; l++){  // phase factor
         complexd zz = double(l+1) + I*eta;
         double sigma = arg(Gamma_Lanczos(zz));  // Coulomb phase shift
         fac_phase[l] = conj( pow(I,l) * exp(I*sigma) );
      }

      for(int l=0; l<l_num; l++){  // interation of l
         Ak[l] = complexd(0.,0.);  // initial value 0 for sum later
      }

      for(int ir=0; ir<nr; ir++){  // interation of r, radial integration

          double r = (ir+0.1) * dr ;  // new grids, must be the same as that in DME_linear

          double x = r*k[ik];
          double fc_array[l_num],  Rcwf[l_num];  // calculate Coulomb wave
          double F_exponent[l_num];
          int status_cwf = gsl_sf_coulomb_wave_F_array (0., l_num-1, eta, x, fc_array, F_exponent );
   
          if(k[ik]==0.005 || k[ik]==0.015 || k[ik]==0.025 || k[ik]==0.035 || k[ik]==0.045)
             output_cwf<< r <<" "; 

          for(int l=0; l<l_num; l++){  // interation of l
 
            if( status_cwf == GSL_SUCCESS) Rcwf[l] = fc_array[l];  // pick up good values of Coulomb wave and save in an array
            else if( status_cwf == GSL_EUNDRFLW) Rcwf[l] = 0.;
            else if( status_cwf == GSL_EOVRFLW)  Rcwf[l] = 1.; // 1.e5; // exp(F_exponent[l]);
            // else Rcwf[l] = 0.;
            else if( fc_array[l] != fc_array[l]) Rcwf[l] = 0.; // 1.; // fc_array[l]=nan
            else Rcwf[l] = fc_array[l];

            // test coulomb wave
            if(k[ik]==0.005 || k[ik]==0.015 || k[ik]==0.025 || k[ik]==0.035 || k[ik]==0.045)
                output_cwf<< Rcwf[l] <<" "; 


            for(int i=0; i<num_basis; i++){  // iteration of continuum basis
              if(l_basis[i]==l && n_basis[i]>n0){  // nonzero for l'=l and exclude bound states
                 int itotal = index2(i,ir,nr);  // total index of new wf
                 complexd value = (c_nl[i][0]+I*c_nl[i][1])*wf_new[itotal]*Rcwf[l]; //(r*r);  // R_nl*R_kl = wf_new/r * fc_array/r
                 // cout<<c_nl[i][0]<<'\t'<<c_nl[i][1]<<'\t'<<wf_new[itotal]<<'\t'<<Rcwf[l]<<endl;
                 Ak[l] += value;  // sum of n and integration of r, no sum of l
            // if(ik==0 && 382<ir && ir<=386) cout << k[i] <<'\t'<<r<<'\t'<<l<<'\t'<<c_nl[i][0]<<'\t'<<c_nl[i][1]<<'\t'<<wf_new[itotal]<<'\t'<<fc_array[l]<<'\t'<<Rcwf[l]<<'\t'<<value<<'\t'<<Ak[l]<<endl;
              }  // end of selection rule
            }  // end of basis

          }  // end of l

          if(k[ik]==0.005 || k[ik]==0.015 || k[ik]==0.025 || k[ik]==0.035 || k[ik]==0.045)
              output_cwf<<endl; 

      }  // end r

      for(int l=0; l<l_num; l++){  // interation of l
        Ak[l] *= dr;
        // if(ik%20==0) cout <<"---  "<<ik <<'\t'<< l <<'\t'<< Ak[l] <<endl;
      }

      for(int it=0; it<n_theta; it++){  // interation of theta

        complexd sum1 = complexd(0.,0.);
        for(int l=0; l<l_num; l++){  // sum of l
          complexd Pl_mom = fac * fac_phase[l]*Ak[l]*Y_lm[it][l] / k[ik];
          output_pl<<real(Pl_mom)<<" "<<imag(Pl_mom)<<" "; //output same k, theta, diff l in one line
          sum1 += fac_phase[l]*Ak[l]*Y_lm[it][l];
          // cout<<fac_phase[l]<<'\t'<<Ak[l]<<'\t'<<Y_lm[it][l]<<endl;
        }

        // cout <<"--- "<< sum1 <<'\t'<< k[ik] <<'\t'<< P_mom[index2(it,ik,nk)] <<endl;
        P_mom[index2(it,ik,nk)] = norm(fac * sum1 / k[ik]);  // 2D momentum distribution
        output_pl<<endl;  // new line for new k, theta

        double factor = k[ik]*sin(theta[it]);   // k_rho
        dP_dk_dkrho[index2(it,ik,nk)]  = abs( factor*P_mom[index2(it,ik,nk)] );  // differential distribution

      }  // end theta

    }  // end of k

 }


}  // end of namepsace



