#ifndef BasisLaser_H
#define BasisLaser_H

#include <math.h>
#include <complex>
#define complexd complex<double>
#include <vector>
#include <iostream>

using namespace std;

namespace Basis_Laser{

 void One_pulse ( double wavelength, double intensity, double n_cycle, vector<double> dt, double cep, double t_wait, vector<double> &time, vector<double> &Efield, vector<double> &Afield, double &period);

 void Two_pulses ( vector<double> wavelength, vector<double> intensity, vector<double> n_cycle, vector<double> dt, vector<double> cep, double delay, vector<double> &time, vector<double> &Efield, vector<double> &Afield, double &period1);

 void SINn_SINm_pulse ( double n_front, double n_back, double wavelength, double intensity, double n_cycle, vector<double> dt, double cep, double t_wait, vector<double> &time, vector<double> &Efield, vector<double> &Afield, double &period);

 void IR_APT ( double lambda1, double i1, double n_c1, double cep1_pi, int nenv, double i2, double n_c2, int nhhg,vector<double> ohhg, vector<double> ahhg, vector<double> cephhg_pi, double delay, double dt, double t_wait, vector<double> &time, vector<double> &E1, vector<double> &E2, vector<double> &Efield, vector<double> &Afield, int &iend_pulses);

 void IR_SAP ( double lambda1, double i1, double n_c1, double cep1_pi, double lambda2, double i2, double n_c2, double cep2_pi, int nenv, double delay, double dt, double t_wait, vector<double> &time, vector<double> &E1, vector<double> &E2, vector<double> &Efield, vector<double> &Afield, int &iend_pulses);

}  // end namespace
#endif


