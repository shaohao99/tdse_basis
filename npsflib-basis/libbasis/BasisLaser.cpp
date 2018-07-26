#include "BasisLaser.h"
#include "BasisConstant.h"

using namespace std;


namespace Basis_Laser{

 void One_pulse (double wavelength, double intensity, double n_cycle, vector<double> dt, double cep, double t_wait, vector<double> &time, vector<double> &Efield, vector<double> &Afield, double &period)
 {

    double amplitude = sqrt( intensity/3.5101e+16 );
    double w = 45.5896/ wavelength;
    period = 2*pi/w;
    double duration = n_cycle*period;

    int N1 = int( duration/dt[0] );
    int N2 = int( t_wait/dt[1] );
    int N = N1+N2;

    time.resize(N);    Efield.resize(N);    Afield.resize(N);
    time[0]=0.;    Efield[0]=0.;    Afield[0]=0.;
    double timestep;

    for(int i=1; i<N; i++){

          if(i<=N1) timestep = dt[0];
          else timestep = dt[1];

          time[i] = time[i-1] + timestep;

          if(time[i]<=duration)  Efield[i] = amplitude * pow(sin(pi*time[i]/duration), 2) * sin(w*time[i] + cep);
          else Efield[i] =0. ;

          Afield[i] = Afield[i-1] - (Efield[i-1] + Efield[i])/2. * timestep * lightC_au;
    }

 }  // end One_pulse


 void Two_pulses ( vector<double> wavelength, vector<double> intensity, vector<double> n_cycle, vector<double> dt, vector<double> cep, double delay, vector<double> &time, vector<double> &Efield, vector<double> &Afield, double &period1)
 {
    double amplitude1 = sqrt( intensity[0]/3.5101e+16 );
    double amplitude2 = sqrt( intensity[1]/3.5101e+16 );
    double w1 = 45.5896/ wavelength[0];
    double w2 = 45.5896/ wavelength[1];
    period1 = 2*pi/w1;
    double period2 = 2*pi/w2;
    double duration1 = n_cycle[0]*period1;
    double duration2 = n_cycle[1]*period2;
    double dt1 = dt[0];
    double dt2 = dt[1];
    double cep1 = cep[0];
    double cep2 = cep[1];

    int N1 = int(duration1/dt1) + 1;
    if( N1 % 2 == 0) N1 +=1;   // N1 must be odd number   
    int N2 = int(duration2/dt2) + 1;

    int N;   // total time grids number
    if(delay <= (duration1 - duration2)/2.) N = N1 + N2 - int(duration2/dt1) + 100;
    else if(delay >= (duration1 + duration2)/2.) N = N1 + N2 + int( (delay - duration1/2. - duration2/2. )/dt1 ) + 100;
    else N = N1 + N2 + 100;

    int istart = int((N1+1)/2) + int((delay-duration2/2.)/dt1);  // starting point of the 2nd pulse
    int iend = istart + N2;  // ending point of the 2nd pulse

    time.resize(N);
    Efield.resize(N);
    Afield.resize(N);
    time[0]=0.;
    Efield[0]=0.;
    Afield[0]=0.;

    double time2=0.;
    double E1=0., E2=0.;

    /*cout<<"pulse1:"<<'\t'<<w1<<'\t'<<period1<<'\t'<<duration1<<'\t'<<dt1<<'\t'<<N1<<endl;
    cout<<"pulse2:"<<'\t'<<w2<<'\t'<<period2<<'\t'<<duration2<<'\t'<<dt2<<'\t'<<N2<<endl;
    cout<<N<<'\t'<<istart<<'\t'<<iend<<endl;*/

    for(int i=1; i<N; i++){

       if(i>=istart && i<=iend){

          time[i] = time[i-1] + dt2;
          time2 += dt2;
          if(time[i]<=duration1)  E1 = amplitude1 * pow(sin(pi*time[i]/duration1), 2) * sin(w1*time[i] + cep1);
          else E1 =0. ;
          if(time2<=duration2)  E2 = amplitude2 * pow(sin(pi*time2/duration2), 2) * sin(w2*time2 + cep2);
          else E2 =0. ;

          Efield[i] = E1 + E2;
          Afield[i] = Afield[i-1] - (Efield[i-1] + Efield[i])/2. * dt2 * lightC_au;
       }
       else{
          time[i] = time[i-1] + dt1;
          // cout<<"---out of the 2d pulse"<<'\t'<<time[i]<<endl;
          if(time[i]<=duration1)  E1 = amplitude1 * pow(sin(pi*time[i]/duration1), 2) * sin(w1*time[i] + cep1);
          else E1 =0. ;

          Efield[i] = E1;
          Afield[i] = Afield[i-1] - (Efield[i-1] + Efield[i])/2. * dt1 * lightC_au;
       }
    }
 }  // end Two_pulses


 void SINn_SINm_pulse ( double n_front, double n_back, double wavelength, double intensity, double n_cycle, vector<double> dt, double cep, double t_wait, vector<double> &time, vector<double> &Efield, vector<double> &Afield, double &period)
 {

    double amplitude = sqrt( intensity/3.5101e+16 );
    double w = 45.5896/ wavelength;
    period = 2*pi/w;
    double duration = n_cycle*period;

    int N1 = int( duration/dt[0] );
    int N2 = int( t_wait/dt[1] );
    int N = N1+N2;

    time.resize(N);    Efield.resize(N);    Afield.resize(N);
    time[0]=0.;    Efield[0]=0.;    Afield[0]=0.;
    double timestep;

    for(int i=1; i<N; i++){

          if(i<=N1) timestep = dt[0];
          else timestep = dt[1];

          time[i] = time[i-1] + timestep;

          if(time[i]<=duration/2.)  Efield[i] = amplitude * pow(sin(pi*time[i]/duration), n_front) * sin(w*time[i] + cep);
          else if(time[i]>duration/2. && time[i]<=duration)  Efield[i] = amplitude * pow(sin(pi*time[i]/duration), n_back) * sin(w*time[i] + cep);
          else Efield[i] =0. ;

          Afield[i] = Afield[i-1] - (Efield[i-1] + Efield[i])/2. * timestep * lightC_au;
    }
 }  // end

// IR + APT
 void IR_APT ( double lambda1, double i1, double n_c1, double cep1_pi, int nenv, double i2, double n_c2, int nhhg,vector<double> ohhg, vector<double> ahhg, vector<double> cephhg_pi, double delay, double dt, double t_wait, vector<double> &time, vector<double> &E1, vector<double> &E2, vector<double> &Efield, vector<double> &Afield, int &iend_pulses)
 {

    double a1 = sqrt( i1/3.5101e+16 );
    double a2 = sqrt( i2/3.5101e+16 );
    double w1 = 45.5896/ lambda1;  // a.u.
    double period1 = 2*pi/w1;
    double duration1 = n_c1*period1;
    double duration2 = n_c2*period1;
    double cep1 = cep1_pi*pi;
    delay = delay*period1;  // delay in a.u.
    t_wait = t_wait*period1;  // delay in a.u.

    vector<double> cephhg(nhhg);
    for(int j=0; j<nhhg; j++){
      cephhg[j]=cephhg_pi[j]*pi;
    }

// count number of points
    int N1 = int(duration1/dt) + 1;
    if( N1 % 2 == 0) N1 +=1;   // N1 must be odd number   
    int N2 = int(duration2/dt) + 1;
    if( N2 % 2 == 0) N2 +=1;   // N2 must be odd number  
    int N3 = int(t_wait/dt) + 1;

    int N;   // total time grids number
    if(abs(delay) <= (duration1 - duration2)/2.) N = N1 + N3 ;
    else if(abs(delay) >= (duration1 + duration2)/2.) N = N1 + N2 + int( (abs(delay) - duration1/2. - duration2/2. )/dt ) + N3;
    else N = int( (duration1/2. + abs(delay) + duration2/2.)/dt ) + N3;
    //N_pulse=N-N3;  // number points during the pulse

    int istart1, iend1, istart2, iend2;   // start and end index of two pulses
    if(delay<-abs(duration1/2.-duration2/2.)) { 
       istart1=int( (duration2/2. + abs(delay) - duration1/2.)/dt );
       istart2=0;
    } else{
       istart1=0;
       istart2=int( (duration1/2. + delay -duration2/2.)/dt );
    }
    iend1 = istart1 + N1 + 1;
    iend2 = istart2 + N2 + 1;
    iend_pulses=max(iend1,iend2);

    time.resize(N);
    E1.resize(N);
    E2.resize(N);
    Efield.resize(N);
    Afield.resize(N);

// initial values
    for(int i=0; i<N; i++){
       time[i] = ( i-(istart1+int((N1+1)/2)) )*dt;  // golable time, 0 at center of IR
       E1[i]=0.;
       E2[i]=0.;
       Efield[i]=0.;
       Afield[i]=0.;
    }

    /*cout<<"pulse1:"<<'\t'<<w1<<'\t'<<period1<<'\t'<<duration1<<'\t'<<N1<<endl;
    cout<<"pulse2:"<<'\t'<<duration2<<'\t'<<N2<<endl;
    cout<<N<<'\t'<<time[istart1]<<'\t'<<time[iend1]<<'\t'<<time[istart2]<<'\t'<<time[iend2]<<'\t'<<endl;*/
  
    //ofstream output_field; 
    //Open_File(output_field, "field_test.dat");

// calculate fields
    for(int i=1; i<N; i++){

       if(i>=istart1 && i<=iend1){
         E1[i]=a1*pow(sin(pi*(time[i]-time[istart1])/duration1),nenv)*sin(w1*(time[i]-time[istart1])+cep1);
       }

       if(i>=istart2 && i<=iend2){
         double carrier=0.;
         for(int j=0; j<nhhg; j++){
           carrier += ahhg[j]*sin( ohhg[j]*w1*(time[i]-time[istart2]) + cephhg[j]); 
         }
         E2[i]=a2*pow(sin(pi*(time[i]-time[istart2])/duration2),2)*carrier;
       }

       Efield[i] = E1[i] + E2[i];
       Afield[i] = Afield[i-1] - (Efield[i-1] + Efield[i])/2. * dt * lightC_au; // A(t)=-c*Int(E(t),t)

       //output_field<<time[i]<<' '<<Efield[i]<<' '<<Afield[i]<<' '<<dt<<' '<<lightC_au<<endl;

    } // end of time   

 }  // end IR_APT

// IR + SAP
 void IR_SAP ( double lambda1, double i1, double n_c1, double cep1_pi, double w2, double i2, double n_c2, double cep2_pi, int nenv, double delay, double dt, double t_wait, vector<double> &time, vector<double> &E1, vector<double> &E2, vector<double> &Efield, vector<double> &Afield, int &iend_pulses)
 {

    double a1 = sqrt( i1/3.5101e+16 );
    double a2 = sqrt( i2/3.5101e+16 );
    double w1 = 45.5896/ lambda1;  // a.u.
    w2 = w2/tran_au_eV;  // a.u.
    double period1 = 2*pi/w1;
    double period2 = 2*pi/w2;
    double duration1 = n_c1*period1;
    double duration2 = n_c2*period2;
    double cep1 = cep1_pi*pi;
    double cep2 = cep2_pi*pi;
    delay = delay*period1;  // delay in a.u.
    t_wait = t_wait*period1;  // delay in a.u.

// count number of points
    int N1 = int(duration1/dt) + 1;
    if( N1 % 2 == 0) N1 +=1;   // N1 must be odd number   
    int N2 = int(duration2/dt) + 1;
    if( N2 % 2 == 0) N2 +=1;   // N2 must be odd number  
    int N3 = int(t_wait/dt) + 1;

    int N;   // total time grids number
    if(abs(delay) <= (duration1 - duration2)/2.) N = N1 + N3 ;
    else if(abs(delay) >= (duration1 + duration2)/2.) N = N1 + N2 + int( (abs(delay) - duration1/2. - duration2/2. )/dt ) + N3;
    else N = int( (duration1/2. + abs(delay) + duration2/2.)/dt ) + N3;
    //N_pulse=N-N3;  // number points during the pulse

    int istart1, iend1, istart2, iend2;   // start and end index of two pulses
    if(delay<-abs(duration1/2.-duration2/2.)) { 
       istart1=int( (duration2/2. + abs(delay) - duration1/2.)/dt );
       istart2=0;
    } else{
       istart1=0;
       istart2=int( (duration1/2. + delay -duration2/2.)/dt );
    }
    iend1 = istart1 + N1 + 1;
    iend2 = istart2 + N2 + 1;
    iend_pulses=max(iend1,iend2);

    time.resize(N);
    E1.resize(N);
    E2.resize(N);
    Efield.resize(N);
    Afield.resize(N);

// initial values
    for(int i=0; i<N; i++){
       time[i] = ( i-(istart1+int((N1+1)/2)) )*dt;  // golable time, 0 at center of IR
       E1[i]=0.;
       E2[i]=0.;
       Efield[i]=0.;
       Afield[i]=0.;
    }

    /*cout<<"pulse1:"<<'\t'<<w1<<'\t'<<period1<<'\t'<<duration1<<'\t'<<N1<<endl;
    cout<<"pulse2:"<<'\t'<<duration2<<'\t'<<N2<<endl;
    cout<<N<<'\t'<<time[istart1]<<'\t'<<time[iend1]<<'\t'<<time[istart2]<<'\t'<<time[iend2]<<'\t'<<endl;*/
  
    //ofstream output_field; 
    //Open_File(output_field, "field_test.dat");

// calculate fields
    for(int i=1; i<N; i++){

       if(i>=istart1 && i<=iend1){
         E1[i]=a1*pow(sin(pi*(time[i]-time[istart1])/duration1),nenv)*sin(w1*(time[i]-time[istart1])+cep1);
       }

       if(i>=istart2 && i<=iend2){
         E2[i]=a2*pow(sin(pi*(time[i]-time[istart2])/duration2),nenv)*sin(w2*(time[i]-time[istart2])+cep2);
       }

       Efield[i] = E1[i] + E2[i];
       Afield[i] = Afield[i-1] - (Efield[i-1] + Efield[i])/2. * dt * lightC_au; // A(t)=-c*Int(E(t),t)

       //output_field<<time[i]<<' '<<Efield[i]<<' '<<Afield[i]<<' '<<dt<<' '<<lightC_au<<endl;

    } // end of time   

 }  // end IR_SAP

}  // end namespace
