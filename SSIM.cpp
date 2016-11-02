#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <math.h>  
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include "SSIM.h"   


using namespace std;
void getLuma(unsigned char *in, unsigned char *out,int size) {
	//#pragma omp parallel for schedule(static, 100)
	for (int i = 0; i < size; i++) {
		//out[i] = round(0.216*in[3 * i] + 0.7152*in[3 * i + 1] + 0.0722*in[3 * i + 2]); //get Luma from RGB picture
		//out[i] = round((double)1/3*in[3 * i] + (double)1/3*in[3 * i + 1] + (double)1/3*in[3 * i + 2]); //get Luma from RGB picture
		out[i] = round(0.299*in[3 * i] + 0.587*in[3 * i + 1] + 0.114*in[3 * i + 2]); //get Luma from RGB picture
	}
}


//count SSIM
double countSSIM(unsigned char * datain1, unsigned char * datain2, int size,int width) { 
	//unsigned char * data1 = (unsigned char*)datain1;
	//unsigned char * data2 = (unsigned char*)datain2;
	double * tmpRes=new double[size];
	unsigned char * data1 = new unsigned char[size];
	unsigned char * data2 = new unsigned char[size];
	getLuma(datain1, data1, size);
	getLuma(datain2, data2,size);
	unsigned char * rect1 = new unsigned char[RECT_SIZE];
	unsigned char * rect2 = new unsigned char[RECT_SIZE];
	int k = 0;
	

	//#pragma omp parallel
	//nthreads = omp_get_num_threads();

//	#pragma omp parallel for schedule(static, 20)
	for (int i = 0; i < size / width - RECT_SQRT; i += SKIP_SIZE) {
		
		for (int j = 0; j < width - RECT_SQRT; j += SKIP_SIZE,k++) {
			//for (int i = 0; i < size-(RECT_SQRT-1)*width; i+=SKIP_SIZE) {
			getRect(data1, i*width+j, width, rect1);
			getRect(data2, i*width + j, width, rect2);
			tmpRes[k] = countRectangle(rect1, rect2);
			//if (tmpRes[k] < 0) cout << "low result: " << i<< ": " << j<< " :" << tmpRes[k] << endl;
			/*delete[] rect1;
			delete[] rect1;*/
		}
	}
	double res= countRes(tmpRes, k);
	delete[] tmpRes;
	delete[] data1;
	delete[] data2;
	delete[] rect1;
	delete[] rect2;
	return res;
}

//count average SSIM value from SSIM values per rectangle
double countRes(double * tmpRes,int count){
	double sum=0;
	for (int i = 0; i < count; i += 1) {
		sum+=tmpRes[i];
		
	}
	return sum/(double)count;
	
}

//return one rectangle with RECT_SIZE pixels
void getRect(unsigned char* data,int start,int width, unsigned char * out ){
	
	for(int i=0;i<RECT_SQRT;i++){
		memcpy(out + i*RECT_SQRT,data+start+i*width,RECT_SQRT);
	}
	//return out;
}

//count ssim of one rectangle with RECT_SIZE pixels
double countRectangle(unsigned char * data1, unsigned char * data2) {
	
	double avg1=countAvg(data1);
	double avg2 = countAvg(data2);

	double var1 = countVariance(data1,avg1);
	double var2 = countVariance(data2,avg2);

	double cov = countCovariance(data1, data2, avg1, avg2);


	double ssim = ((2*avg1*avg2+C1)*(2*cov+C2)) / ((avg1*avg1+ avg2*avg2 +C1)*(var1+var2+C2));
	return ssim;
}
//count avg value of given rectangle 
double countAvg(unsigned char * data) {
	double avg=0;
	for (int i = 0; i < RECT_SIZE; i++) {
		avg += data[i];
	}
	avg = avg / (double)RECT_SIZE;
	return avg;
}

//count variance of given rectangle
double countVariance(unsigned char * data, double avg) {
	double var = 0;
	for (int i = 0; i < RECT_SIZE; i++) {
		var += (data[i]-avg)*(data[i] - avg);
	}
	var = var/ (double)RECT_SIZE ;
	return var;
}

//count covariance of given rectangle
double countCovariance(unsigned char * data1, unsigned char * data2, double avg1, double avg2) {
	double cov = 0;
	for (int i = 0; i < RECT_SIZE; i++) {
		cov += (data1[i] - avg1)*(data2[i] - avg2);
	}
	cov = cov/ (double)RECT_SIZE;
	//if (cov < 0) cout << "neg "<<cov << endl;
	return cov;
}