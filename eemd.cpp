// eemd.cpp 

#include <cstdlib>
#include <cmath>
#include <ctime>
#include "emd_core.cpp"
#include <iostream>

void randn(double *W, int sz) {
	srand((int)time(NULL));
	//srand(0);
	double u = 0;
	double v = 0;
	for (int i = 0; i < sz; i++) {
		u = (rand() *1.0) / RAND_MAX;
		v = (rand() *1.0) / RAND_MAX; 
		W[i] = sqrt( -2 * log(u) ) * cos(2*3.1415926*v);
	}
}

double Std(double *Y, int sz) {
	double mean = 0;
	for (int i = 0; i < sz; i++)
		mean = mean + Y[i];
	mean = mean / sz;
	double sigma = 0;
	for (int i = 0; i < sz; i++) {
		sigma = sigma + pow( (Y[i] - mean) , 2);
		// if (sigma < 0) {
		// 	std::cout << "Warning! sigma = " << sigma << " at position " << i << endl;
		// 	system("PAUSE");
		// }
	}
	if (sigma > 0)
		sigma = sqrt( sigma / sz );
	else {
		// std::cout << "Error! sigma^2 = " << sigma << " before taking square root!" << endl;
		// system("PAUSE");
		sigma = 0;
	}
	//std::cout << "sigma = " << sigma << endl; 
	return sigma;
}

void eemd(double *modes,
		double *Y, int sz, int goal, int ens, double nstd) {

	int m, i, c, k, t;
	

	/* parse input arguments */

	int goal1 = goal + 1;
	int MAX = sz;
	int GMX = goal1;

	/* create output arguments */

	/* Core function */
	double *m1 = new double[GMX*MAX];
	double *m2 = new double[GMX*MAX];
	double *tmp = new double[GMX*MAX];
	double *wn = new double[MAX];
	double *Y1 = new double[MAX];
	double *Y2 = new double[MAX];
	double sigma = Std(Y, sz);
    //printf("Initializing... \n");
    for (t = 0; t < sz*goal1; t++)
        tmp[t] = 0;
	for (k = 0; k < ens; k++) {
		randn(wn, sz);
		for (i = 0; i < sz; i++) {
			if (sigma > 0)
				Y1[i] = Y[i]/sigma + wn[i]*nstd;
			else
				Y1[i] = 0 + wn[i]*nstd;
			if (nstd > 0)
				if (sigma > 0)
					Y2[i] = Y[i]/sigma - wn[i]*nstd;
				else
					Y2[i] = 0 - wn[i]*nstd;
		}

		emd_core(m1, Y1, sz, goal);
		if (nstd > 0)
			emd_core(m2, Y2, sz, goal);

		for (t = 0; t < sz*goal1; t++) {
            //printf("%f <= %f \n", tmp[t], m1[t] + m2[t]);
			tmp[t] = tmp[t] + m1[t];
			if (nstd > 0)
				tmp[t] = tmp[t] + m2[t];
        } // END of for-t
	} // end of for-k
	if (nstd > 0)
		for (t = 0; t < sz*goal1; t++)
			modes[t] = tmp[t]*sigma / (ens * 2);
	else
		for (t = 0; t < sz*goal1; t++)
			modes[t] = tmp[t]*sigma / ens;
	delete[] m1, m2, tmp, wn, Y1, Y2;

} // end of eemd
	
