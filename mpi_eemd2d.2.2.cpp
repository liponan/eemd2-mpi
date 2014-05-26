// mpi_eemd2d.cpp
// written by Neil Po-Nan Li 2012/11/23 @ IPAS
// modified 2013/12/20 : For better stablity, only the root node reads the input txt file.
// v.2.0    2014/04/30 : Interpolation engine changed to GSL
// v.2.1    2014/05/05 : Better memory ultilization
// v.2.2 	2014/05/07 : use binary file for export
// v.2.3 	2014/05/26 : read binary file as input



#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <string>
#include <fstream>
#include <mpi.h>

#include "eemd.cpp"

	using namespace std;

	int readBinaryHeader(int* dim, int* lg, string filename);
	int readBinaryImage(double *Y, string filename);
	void writeBinary(string filename, int dim, int* lg, double *Y);
	int toDo(int, int, int);

int main(int argc, char *argv[])
{
	// initialize MPI
	MPI_Init(&argc, &argv);
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	int t1 = MPI_Wtime();
	int t2, t3, t4, t5;
	int dt, eta_time;
	int dim = 0;
	int d = 0;
	int lg[3] = {0};
	bool flag = true;
	int H, W, SZ;

	int bin_flag1 = 1;
	int bin_flag2 = 1;

	if (world_rank == 0) {
		
		// read the input file		
		if (argc < 2) {
			cout << "No input file!!" << endl;
			return 0;
		} else {
			cout << "Loading " << argv[1] << endl;
			bin_flag1 = readBinaryHeader(&dim, lg, argv[1]);
			if (bin_flag1 == 0)
				flag = false;
			else {
				cout << "Can't open " << argv[1] << ". Please try again... " << endl;
				return 1;
			}
		}
			
		// read the header
		W = lg[1];
		H = lg[0];
		SZ = H * W;
	} // end of if (world_rank == 0)

	MPI_Bcast(&SZ,  1, MPI_INT, 0, MPI_COMM_WORLD); // broadcast the data length SZ
	MPI_Bcast(&dim, 1, MPI_INT, 0, MPI_COMM_WORLD); // broadcast the dimension dim
	MPI_Bcast(&H,   1, MPI_INT, 0, MPI_COMM_WORLD); // broadcast the height H
	MPI_Bcast(&W,   1, MPI_INT, 0, MPI_COMM_WORLD); // broadcast the width W

	// print the dimensions
	if (world_rank == 0) {
		cout << "# of dimensions: " << dim << endl;
		cout << "Size: [ ";
		for (d = 0; d < dim; d++)
			cout << lg[d] << " ";
		cout << "]" << endl;
	} // end of if (world_rank == 0)

	int i = 0, j = 0, k = 0;
	int t = 0, m = 0;

	// Now we know the exact data size, so let's declare the memory block for the data
	double *img = NULL;
	double tmpNum;

	// load the image into memory
	if (world_rank == 0) {
		img = new double[SZ];
		bin_flag2 = readBinaryImage(img, argv[1]);
	} // end of if (world_rank == 0)

	//MPI_Bcast(img, SZ, MPI_DOUBLE, 0, MPI_COMM_WORLD); // broadcast the IMG data to all nodes

	// find the parameters or use the default values
	int goal = 3;
	int ens = 1;
	double nstd = 0;
	

	int timecode =  (int)time(NULL) % 10000;

	// You all know the argument values, so I don't need to broadcast them
	switch (argc) {
		case 6: 
			timecode = atoi(argv[5]);
		case 5:
			nstd = atof(argv[4]);
		case 4:
			ens = atoi(argv[3]);	
		case 3:
			goal = atoi(argv[2]);
		default:
			break;
	}
	int goalt = goal + 1;

	char timecode_str[5];
	sprintf(timecode_str, "%d", timecode);
	
	double *modes1 = NULL; 
	double *rootBuff1 = NULL;
	if (world_rank == 0) {
		cout << "# of modes: " << goal << endl;
		cout << "# of ensembles: " << ens << endl;
		cout << "Amplitude of white noise:" << nstd << endl;
		rootBuff1 = new double[SZ*goalt];
		modes1    = new double[SZ*goalt]; 
		cout << "EEMD for dimension #1 starts!" << endl;
	}


	/* now establish the plan for parallel programing */
	int *mCnts = new int[world_size];
	int *nCnts = new int[world_size];
	int *mDisps = new int[world_size];
	int *nDisps = new int[world_size];
	mCnts[0] = toDo(H, 0, world_size) * W;
	nCnts[0] = toDo(W, 0, world_size) * H;
	mDisps[0] = 0;
	nDisps[0] = 0;
	// cout << "world_size=" << world_size << endl;
	for (t = 1; t < world_size; t++) {
		mCnts[t] = toDo(H, t, world_size) * W;
		nCnts[t] = toDo(W, t, world_size) * H;
		mDisps[t] = mDisps[t-1] + mCnts[t-1];
		nDisps[t] = nDisps[t-1] + nCnts[t-1];
		// cout << mCnts[t] << " " << nCnts[t] 
		// 	<< " " << mDisps[t] << " " << nDisps[t] << endl;
	}

	// new plan for gathering post-eemd data
	int *uCnts  = new int[world_size];
	int *vCnts  = new int[world_size];
	int *uDisps = new int[world_size];
	int *vDisps = new int[world_size];
	uCnts[0] = toDo(H, 0, world_size) * W * goalt;
	vCnts[0] = toDo(W, 0, world_size) * H * goalt;
	uDisps[0] = 0;
	vDisps[0] = 0;

	for (t = 1; t < world_size; t++) {
		uCnts[t] = toDo(H, t, world_size) * W * goalt;
		vCnts[t] = toDo(W, t, world_size) * H * goalt;
		uDisps[t] = uDisps[t-1] + uCnts[t-1];
		vDisps[t] = vDisps[t-1] + vCnts[t-1];
	}

	/* 2D-EEMD core function */
	
	// for #1 dimension: solve for each ROW
	int myTodo = toDo(W, world_rank, world_size);
	// cout << "Node " << world_rank << " gets " << myTodo << " jobs to do..." << endl;
	double *myBuff1  = new double[nCnts[world_rank]];
	double *myModes1 = new double[nCnts[world_rank]*goalt];
	double *inTmp1  = new double[H];
	double *outTmp1 = new double[H * goalt];
	MPI_Scatterv(img, nCnts, nDisps, MPI_DOUBLE, 
		myBuff1, nCnts[world_rank], MPI_DOUBLE,
		0, MPI_COMM_WORLD);
	for (t = 0; t < myTodo; t++) {
		for(i = 0; i < H; i++) 
			inTmp1[i] = myBuff1[i + t*H];
		eemd(outTmp1, inTmp1, H, goal, ens, nstd);
		for (i = 0; i < H*goalt; i++)
			myModes1[i + t * H*goalt] = outTmp1[i];
	}
	// char metaFilename[10];
	// sprintf(metaFilename,"meta%d.m", world_rank);
	// string metaStr(metaFilename);
	// int lg0[3] = {myTodo, H, goalt};
	// printArray(metaFilename, myModes1, 3, lg0);
	delete[] img, myBuff1, inTmp1, outTmp1;

	// gather and re-sort the data: rootBuff1 => modes1
	MPI_Gatherv(myModes1, vCnts[world_rank], MPI_DOUBLE, 
		rootBuff1, vCnts, vDisps, MPI_DOUBLE,
		0, MPI_COMM_WORLD);
	if (world_rank == 0) {
		for (k = 0; k < goalt; k++)
			for (i = 0; i < H; i++)
				for (j = 0; j < W; j++)
					modes1[i + H*(j + k*W)] = rootBuff1[i + H*(k + j*goalt)];
	//string metaFilename = "meta1.m";
	//int lg0[3] = {H, W, goalt};
	//printArray(metaFilename, modes1, 3, lg0);
	}
	delete[] myModes1, rootBuff1;
	t2 = MPI_Wtime();
	dt = t2 - t1;
	if (world_rank == 0)
		cout << "EEMD stage 1 done in " << dt << "s" << endl;
	

	// for #2 dimension: solve for each COL in each MODE
	//MPI_Bcast(modes1, SZ*goalt, MPI_DOUBLE, 0, MPI_COMM_WORLD); // seed for 2nd EEMD phase
	myTodo = toDo(H, world_rank, world_size);
	// cout << "Node " << world_rank << " gets " << myTodo << " jobs to do..." << endl;
	// buffers for each mode only
	double *modeBuff2in = new double[SZ]; // buffer for each mode from modes1
	double *myBuff2 = new double[mCnts[world_rank]];
	double *myModes2 = new double[mCnts[world_rank]*goalt];
	double *inTmp2 = new double[W];
	double *outTmp2 = new double[W * goalt];
	// buffer for inter-mode data collecting
	double *modeBuff2out = NULL;
	// buffers for post-EEMD data collecting
	double *modes2 = NULL;
	// memory allocation for root-only arrays`
	if (world_rank == 0) {
		modeBuff2out = new double[SZ * goalt];
		modes2 = new double[SZ * goalt*goalt];
	}
	// parallel in each mode
	for (m = 0; m < goalt; m++) {
		t3 = MPI_Wtime();
		if (world_rank == 0)
			for (i = 0; i < H; i++)
				for (j = 0; j < W; j++)
					modeBuff2in[j + i*W] = modes1[i + H*(j + m*W)];
		MPI_Scatterv(modeBuff2in, mCnts, mDisps, MPI_DOUBLE, 
			myBuff2, mCnts[world_rank], MPI_DOUBLE,
			0, MPI_COMM_WORLD);
		for (t = 0; t < myTodo; t++) {
			for (j = 0; j < W; j++)
				inTmp2[j] = myBuff2[j + t*W];
			eemd(outTmp2, inTmp2, W, goal, ens, nstd);
			for (j = 0; j < W*goalt; j++)
				myModes2[j + t * W*goalt] = outTmp2[j];
		}
		MPI_Gatherv(myModes2, uCnts[world_rank], MPI_DOUBLE, 
			modeBuff2out, uCnts, uDisps, MPI_DOUBLE,
			0, MPI_COMM_WORLD);
		if (world_rank == 0) {
			for (i = 0; i < H; i++)
				for (j = 0; j < W; j++)
					for (k = 0; k < goalt; k++)
						modes2[i + H*(j + W*(m + k*goalt))]
						 = modeBuff2out[j + W*(k+ i*goalt)];
			
		t4 = MPI_Wtime();
		dt = t4 - t3;
		eta_time = (t4 - t2) * (goalt-m-1) / (m+1);
		cout << "Mode " << m+1 << "/" << goal << " solved in " << dt << "s.  ";
		cout << eta_time << "s to go..." << endl; 
		}
	}
	delete[] modes1;
	delete[] modeBuff2in, myBuff2, myModes2, inTmp2, outTmp2;
	delete[] modeBuff2out;
	if (world_rank == 0) {
		dt = t4 - t2;
		cout << "EEMD stage 2 done in " << dt << "s" << endl;
		cout << "Combining modes... " << endl;
	}

	// combine modes
	double *modes = NULL;
	if (world_rank == 0) {
		modes = new double[SZ * goalt];
		for (i = 0; i < H; i++)
			for (j = 0; j < W; j++) {
				for (m = 0; m < goalt; m++) {
					modes[i + H*(j + m*W)] = 0;
					for (k = m; k < goalt; k++) {
						modes[i + H*(j + m*W)] += modes2[i + H*(j + W*(k + m*goalt))];
						modes[i + H*(j + m*W)] += modes2[i + H*(j + W*(m + k*goalt))]; 
					}
					modes[i + H*(j + m*W)] -=  modes2[i + H*(j + W*(m + m*goalt))]; 
				}
			}
	}
	delete[] modes2;


	// export the result to file
	t4 = MPI_Wtime();
	if (world_rank == 0) {
		cout << "Done! Writing to file... " << endl;
		int dim2 = 3;
		int lg2[3] = {H, W, goalt};
		//int timecode = (int)time(NULL) % 10000;
		//char timecode_str[4];
		//sprintf(timecode_str, "%d", timecode);
		string filenameStr(argv[1]);
		string filename_export // v
		 = string(filenameStr,0,filenameStr.length()-4)+"_modes" + timecode_str + ".bin";
		string filename_log    // v
		 = string(filenameStr,0,filenameStr.length()-4)+"_log" + timecode_str + ".txt";
		writeBinary(filename_export, dim2, lg2, modes);
		cout << filename_export << " exported!" << endl;
		ofstream fout(filename_log.c_str());
		fout << "Input: " << filenameStr << endl;
		fout << "# of modes: " << goal << endl;
		fout << "# of ensembles: " << ens << endl;
		fout << "Amp. of white noise: " << nstd << endl;
		fout << "# of cores: " << world_size << endl;
		fout << "Elapsed time: " << t4 - t1 << "s" << endl;
		fout.close();
		cout << "Elapsed time: " << t4 - t1 << "s" << endl;
	}
	MPI_Finalize();
	delete[] modes;
	return 0;
}

int readBinaryHeader(int* dim, int* lg, string filename) {
	FILE *file;
	char filename_char[20];
	strcpy(filename_char, filename.c_str()); 
	file = fopen(filename_char , "rb");

	// first byte: dimension number
	fread(dim, sizeof(int), 1, file);

	// 2nd~4th bytes: size in each dimension
	for (int i = 0; i < *dim; i++) {
		fread(&lg[i], sizeof(int), 1, file);
	} // end of for
	return 0;
}

int readBinaryImage(double *Y, string filename) {
	FILE *file;
	char filename_char[20];
	strcpy(filename_char, filename.c_str()); 
	file = fopen(filename_char , "rb");

	int dim = 0;
	int lg[3] = {0};

	// first byte: dimension number
	fread(&dim, sizeof(int), 1, file);

	// 2nd~4th bytes: size in each dimension
	int sz = 1;
	for (int i = 0; i < dim; i++) {
		sz *= lg[i];
		fread(&lg[i], sizeof(int), 1, file);
	} // end of for

	// remaing bytes: data (double)
	fread(Y, sizeof(double), sz, file);
	return 0;
}

void writeBinary(string filename, int dim, int* lg, double *Y) {
	
	FILE *file;
	char filename_char[20];
	strcpy(filename_char, filename.c_str()); 
	file = fopen(filename_char , "wb");
	rewind (file);
	
	// first byte: dimension number
	fwrite(&dim, sizeof(int), 1, file);

	// 2nd~4th bytes: size in each dimension
	int sz = 1;
	for (int i = 0; i < dim; i++) {
		sz *= lg[i];
		fwrite(&lg[i], sizeof(int), 1, file);
	}

	// remaing bytes: data (double)
	fwrite(Y, sizeof(double), sz, file);	

	// close file
	fclose(file);
}

int toDo(int N, int myrank, int world_size) {
	int num = floor((N*1.0)/world_size);
	int remaining = N - num * world_size;
	if (myrank >= world_size - remaining)
		num++;
	return num;
}
