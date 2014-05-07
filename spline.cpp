// spline.cpp : 定義 DLL 應用程式的匯出函式。
//

#include "spline.hpp"
#include <iostream>
#include <fstream>
	//using namespace std;

void spline(double *YY, 
		int *X, double *Y, int m1, int m2) {
	// m1: length of discrete extremas
	// m2: length of original data points
	

	/* Core function */
	//fstream fout("spline_log.txt", ios::out | ios::app);

	using namespace magnet::math;

	Spline spline;
	double m;
	
	if (m1 > 2 ) {
		//fout << ">";
		for (int i = 0; i < m1; i++) {
			spline.addPoint(X[i],Y[i]);
			//fout << "\t" << X[i] << ", " << Y[i];
		}
		//fout << endl << "yy = [";
		for (int j = 0; j < m2; j++) {
			//cout << ">>\t" << j << ", ";
			YY[j] = spline(j);
			//fout << YY[j] << " ";
			//cout << YY[j] << endl;
		}
	}
	else {
		m = (Y[1] - Y[0]) / (m2 - 1);
		//fout << "m = " << Y[1] << " - " << Y[0] << " / " << (m2-1) << " = " << m<< endl; 
		//fout << "yy' = [";
		for (int j = 0; j < m2; j++) {
			//cout << ">>\t" << j << ", ";
			YY[j] = Y[0] + m * j;
			//fout << YY[j] << " ";
			//cout << YY[j] << endl;
		}
	}
	//fout << "];" << endl;
} // END of spline
