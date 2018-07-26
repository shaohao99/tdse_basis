//Constants
#ifndef BASISCONSTANT_H
#define BASISCONSTANT_H
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <vector>
#include <string>	
#include <complex>
#define complexd complex<double>

const std::complexd I=std::complexd(0.,1.); 
const double lightC_au = 137.036;
const double one_by_lightC_au = 1./lightC_au;
const double pi = 3.141592653589793238463383;
const double charge_e1_au = -1.;
const double charge_e2_au = -1.;
const double mass_e1_au = 1.;
const double mass_e2_au = 1.;
const double charge_n1_au=1.;
const double charge_n2_au=1.;
const double mass_n1_au=1836.;
const double mass_n2_au=1836.;
const double tran_au_eV=27.2113845;
const double tran_au_wn=219474.63;
const double tran_au_fs=0.024188843265;

const double e_Euler = 2.71828247254;

enum gauge_t { lengthgauge, velocitygauge, othergauge };

#endif	/* BASISCONSTANT_H */
