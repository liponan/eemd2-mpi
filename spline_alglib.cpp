// spline.cpp : 定義 DLL 應用程式的匯出函式。
//

#include "interpolation.h"

using namespace alglib;

void spline(double *YY, 
		int *X, double *Y, int m1, int m2) {

	// m1: length of discrete extremas
	// m2: length of original data points

	/* Core function */

	double* dx = new double[m1];
	spline1dinterpolant s;
	real_1d_array rx;
	real_1d_array ry;

	for (int j = 0; j < m1; j++)
		dx[j] = X[j];

	rx.setcontent(m1, dx);
	ry.setcontent(m1, Y);


	double m;
	
	if (m1 > 2 ) {

		spline1dbuildcubic(rx, ry, s);

		for (int j = 0; j < m2; j++) {
			YY[j] = spline1dcalc(s, j);
		} // end of for-j

	} // end of if

	else {

		m = (Y[1] - Y[0]) / (m2 - 1);

		for (int j = 0; j < m2; j++) {
			YY[j] = Y[0] + m * j;		
		} // end of for-j
		
	} //end of else

}
