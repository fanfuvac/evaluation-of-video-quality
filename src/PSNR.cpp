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


	for (int i = 0; i < size; i++) {
		pixel = ((unsigned char)data1[i] - (unsigned char)data2[i])*((unsigned char)data1[i] - (unsigned char)data2[i]);
		//pixel = pow((unsigned char)data1[i] - (unsigned char)data2[i],2);
		//cout<<int((unsigned char)frame->data[i])<<endl;
		//#pragma omp atomic
		sum += pixel;
	}



	//for (int i = 0; i < 8; i++) sum += subSum[i];
	double max = pow(2, 8) - 1;
	double mse = sum / size;
	double psnr = 10 * log10(max*max / mse);
	//cout << nthreads << endl;
	return psnr;
}
