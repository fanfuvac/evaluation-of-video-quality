#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <math.h>  
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include "psnr.h"   

#include <omp.h>

using namespace std;

#define CHUNKSIZE   100

//#define LINUX
using namespace std;

#define header_size 14+40
/*

//gets array of data from BMP file, skipping the header
PictureData getData(string filename){
	ifstream frame;
	frame.open (filename.c_str(),ios::in|ios::binary);
	//frame.seekg ( 14+4 );
	PictureData data;
	frame.seekg (0, ios::end);
	data.size = frame.tellg();
	data.size-=header_size;
	data.data = new  char [data.size];
	frame.seekg (14+sizeof(int), ios::beg);
	char * tmp=new char[sizeof(int)];
	frame.read (tmp, sizeof(int));
	memcpy(&data.width,tmp,sizeof(int));
	
	frame.seekg (14+2*sizeof(int), ios::beg);
	frame.read (tmp, sizeof(int));
	memcpy(&data.height,tmp,sizeof(int));
	
	frame.seekg (header_size, ios::beg);
    frame.read (data.data, data.size-header_size);
    frame.close();
	return data;
	//out = new char[]
}*/

double countPSNR(unsigned char * data1, unsigned  char * data2, int size){
	double sum = 0;
	int pixel;
	int nthreads;
	omp_set_num_threads(1);
	double subSum[8] = { 0 };
	//#pragma omp parallel
	{
		
		//#pragma omp master
		nthreads =omp_get_num_threads();
		
		//#pragma omp  for 
		for (int n = 0; n < 8; ++n) {
			int start = size * 3 / 8 * n;
			for (int i = start; i < start + size * 3 / 8; i++) {
				pixel = ((unsigned char)data1[i] - (unsigned char)data2[i])*((unsigned char)data1[i] - (unsigned char)data2[i]);
				//pixel = pow((unsigned char)data1[i] - (unsigned char)data2[i],2);
				//cout<<int((unsigned char)frame->data[i])<<endl;
				//#pragma omp atomic
				subSum[n] += pixel;
			}
		}
	}

	for (int i = 0; i < 8; i++) sum += subSum[i];
	double max = pow(2, 8) - 1;
	double mse = sum / size / 3;
	double psnr = 10 * log10(max*max / mse);
	//cout << nthreads << endl;
	return psnr;
}
