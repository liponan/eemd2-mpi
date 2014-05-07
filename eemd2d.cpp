
/*

#include <math.h>
#include <cstdio>
#include <cstdlib>
*/

#include <iostream>
#include <ctime>
#include "eemd.cpp"

	using namespace std;

void eemd2d(double *modes,
		double *Img, int U, int V, int goal, int ens, double nstd) {
	int m, i, j, c, k;

	int goal1 = goal + 1;
	int sz = U * V;

	int GMX = goal1;

	/* Core function */

	//ofstream fout("spline_log.txt");

	// Solve through dim#1
	double *tmp1 = new double[V];
	double *tmp2 = new double[U];
	double *tm1 = new double[GMX*V];
	double *tm2 = new double[GMX*U];
	double *modes1 = new double[GMX*U*V];
	double *modes3 = new double[GMX*GMX*U*V];

    
	int time0, time1, time2;
	double eta;
	time0 = (int)time(NULL);
	for (i = 0; i < U; i++) {
		//printf("loading... ");
        for (j = 0; j < V; j++) 
			tmp1[j] = Img[i + j * U]; // end of for-j
		time1 = (int)time(NULL);
		//fout << "Solving row " << i+1 << "/" << U << "..." << endl;
        cout << "Solving row " << i+1 << "/" << U << "..." << endl;
		eemd(tm1, tmp1, V, goal, ens, nstd);
		/*for (int t = 0; t < V*goal1; t++)
			tm1[t] = 0;*/
		//emd_core(tm1, tmp1, V, goal);
		time2 = (int)time(NULL);
		eta = (time2 - time0) * (U - i - 1.0) / (i + 1.0);
		cout << "Done in " << time2 - time1 << "s. ";
		cout << (i+1)*100.0/U << "% solved. " << eta << "s to go... " << endl;
		for (j = 0; j < V; j++) {
			for (m = 0; m < goal1; m++) 
				modes1[i + U * (j + m * V )] = tm1[j + m * V];
		}
	} // end of for-i
	delete [] tmp1;
	delete [] tm1;

	int time00 = (int)time(NULL);
	// Solve through dim#2
	for (m = 0; m < goal1; m++) {
		time1 = (int)time(NULL);
		cout << "Solving mode " << m+1 << "/" << goal1 << "..." << endl;
		for (j = 0; j < V; j++) {
			for (i = 0; i < U; i++)
				tmp2[i] = modes1[i + U * (j + m * V)];
			
			//fout << "Solving col " << j+1 << "/" << V 
			//	<< " in mode " << m+1 << "/" << goal1 << "..." << endl;
			
			/*for (int t = 0; t < U*goal1; t++)
				tm2[t] = 0;*/
			//emd_core(tm2, tmp2, U, goal);
			eemd(tm2, tmp2, U, goal, ens, nstd);
			
			for (i = 0; i < U; i++)
				for (k = 0; k < goal1; k++)
					modes3[i + U * (j + V * (m + k * goal1))] = tm2[i + k*U];
		} // end of for-j
		time2 = (int)time(NULL);
		eta = (time2 - time00) * (goal1 - m - 1.0) / (m + 1.0);
		cout << "Done in " << time2 - time1 << "s. ";
		cout << (m+1)*100.0/(goal1) << "% solved. " << eta << "s to go... " << endl;
	} // end of for-m
	delete[] tm2, tmp2, modes1;

	// Combine modes
	for (int t = 0; t < U*V*goal1; t++)
		modes[t] = 0;
	for (i = 0; i < U; i++) {
		for (j = 0; j < V; j++) {
			for (m = 0; m < goal1; m++) {
				for (k = m; k < goal1; k++) {
					modes[i + U * (j + m * V )] += modes3[i + U * (j + V * (m + k * goal1))];
					modes[i + U * (j + m * V )] += modes3[i + U * (j + V * (k + m * goal1))];
				} // end of for-k
				modes[i + U * (j + m * V )] -= modes3[i + U * (j + V * (m + m * goal1))];
			} // end of for-m
		} // end of for-j
	} //end of for-i
	delete[] modes3;

	cout << "Elasped time: " << (int)time(NULL) - time0 << "s. " << endl;
	} // end of eemd2d
	
