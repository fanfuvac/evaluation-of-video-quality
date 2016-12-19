#pragma once
double countSSIM(unsigned char * datain1, unsigned  char * datain2, int size,int width); 
double countRectangle(unsigned char * data1, unsigned  char * data2);
double countAvg(unsigned char * data);
double countVariance(unsigned char * data, double avg);
double countCovariance(unsigned char * data1, unsigned  char * data2, double avg1, double avg2);
void getRect(unsigned char* data, int start, int width, unsigned char * out);
double countRes(double * tmpRes, int count);


#define RECT_SIZE 64 //number of pixels in rectangle
#define RECT_SQRT 8 //size of side of rectangle
#define C1 6.5025 //(K1*L )**2
#define C2 58.5225 //(K2*L )**2
#define SKIP_SIZE 8 //how many pixels should be skipped when going to next window, should be pow of 2((have to be?)
//#define CHUNK_SIZE 4 //how many frames will be executed at once

extern  int CHUNK_SIZE;