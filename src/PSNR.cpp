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



//#define LINUX
using namespace std;



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
