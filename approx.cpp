#include <stdio.h>
#include <math.h>

int main() {
	
	double j = 0.1;				// dimensionless angular momentum j = J / M^2
	double q = 3*j*j;			// dimensionless quadrupole moment q = Q / M^3
	double M = 1.4;				// mass in solar mass units  
	
	double b = 1;				// beta to beta_cusp ratio - from 0 to 1
	
    const double c = 2.99792458e10;			// speed of light in [cm/s]
    const double G = 6.67428e-8;			// gravitational constant in [cm^3/(g*s^2)]
    const double Msun = 1.98892e33;			// solar mass in [g]
    double R = pow(c, 3) / (G * M * Msun);		// unit conversion
    
    double nu0 = 1 / (2 * M_PI) * R / pow(6, 1.5);

	double nuk;					// Keplerian frequency = variable

	double r0;
	
		r0 = pow((pow(6, 1.5) - j * nuk / nu0 + 0.5 * (q - pow(j, 2)) * pow(nuk / nu0, 2)) / (nuk / nu0), 2.0 / 3.0);
	
	double RR,BB,KK,TT;
	
		RR = 1 - b*b * (1.0 / 33.0 + 10.0 / nu0 * pow(pow(r0, 1.5) - 18.4 + 15 * j - 4 * q,2) / (1 - 1.5 * j + 0.6 * q));
		BB = 1 - 0.2 * b*b * (1 + j);
		KK = 1 - b*b * (14.0 / nu0 * pow(pow(r0, 1.5) - 40.0 / 3.0 + 10 * j - 2.5 * q,2) / (1 - 1.7 * j + 0.7 * q));
		TT = 1 - b*b * (5 * j / nu0 * (1 - 2 / (1 - 2 * j - 2 * pow(j, 2) + 3 * q) * pow(r0 - 6.1 + 2.5 * j, 2)) + 9 * (q - pow(j, 2)) / nu0 * (r0 - 7.2));

	double radaxisym, radprec, vertaxisym, vertprec;
	
		radaxisym = nuk * RR * sqrt(1 - 6/r0 + 8 * j/pow(r0, 1.5) + 57 * j*j/pow(r0, 2.5) - 3 * q/pow(r0, 2) * (1 + 19/pow(r0, 1.5)));
		radprec = nuk * (1 - BB * sqrt(1 - 6/r0 + 8 * j/pow(r0, 1.5) + 57 * j*j/pow(r0, 2.5) - 3 * q/pow(r0, 2) * (1 + 19/pow(r0, 1.5))));
		vertaxisym = nuk * KK * sqrt(1 - 4 * j/pow(r0, 1.5) - 24 * j*j/pow(r0, 2.5) + 3 * q /pow(r0, 2) * (1 + 8/pow(r0, 1.5)));
		vertprec = nuk * (1 - TT * sqrt(1 - 4 * j/pow(r0, 1.5) - 24 * j*j/pow(r0, 2.5) + 3 * q /pow(r0, 2) * (1 + 8/pow(r0, 1.5))));
			
  return 0;
} 
