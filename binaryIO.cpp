#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string.h>
#include <fstream>

	using namespace std;

int readBinaryHeader(int* dim, int* lg, string filename);
int readBinaryImage(double *Y, string filename);
void writeBinary(string filename, int dim, int* lg, double *Y);

int readBinaryHeader(int *dim, int *lg, string filename) {
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
	fclose(file);

	return 0;
}

int readBinaryImage(double *Y, string filename) {
	FILE *file;
	char filename_char[20];
	strcpy(filename_char, filename.c_str()); 
	file = fopen(filename_char , "rb");
	rewind(file);

	int dim = 0;
	int lg[3] = {1};

	// first byte: dimension number
	fread(&dim, sizeof(int), 1, file);

	// 2nd~4th bytes: size in each dimension
	int sz = 1;
	for (int i = 0; i < dim; i++) {
		fread(&lg[i], sizeof(int), 1, file);
		sz *= lg[i];	
	} // end of for

	// remaing bytes: data (double)
	fread(Y, sizeof(double), sz, file);
	
	fclose(file);
	return 0;
}

void writeBinary(string filename, int dim, int *lg, double *Y) {
	
	FILE *file;
	char filename_char[20];
	strcpy(filename_char, filename.c_str()); 
	file = fopen(filename_char , "wb");
	
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