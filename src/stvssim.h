#pragma once
#include "main.h"   

#define RECT_SIZE_ARPS 64 //number of pixels in rectangle
#define RECT_SQRT_ARPS 8 //size of side of rectangle
#define ZERO_MVMT 12 //level of SAD for zero movement prejudgement(2 bit per pixel)
#define FRAME_CNT 33 //count of frames being evaluated by SSIM-3D, need to be odd
#define ROOD_SIZE 2 //size of the rood for ARPS
#define RECT_SQRT_3D 11 //size side of rectangle for 3D-SSIM
#define RECT_SIZE_3D 11*11 //size of rectangle for 3D-SSIM
#define FRAME_SKIP 16 //how mane frames will be skipped between 2 frames (every how many frames is stVSSIM being counted)
#define INT_MAX 2147483647

using namespace std;
//vector used for ARPS
struct vector { 
	int x;
	int y;
};
double ** countMetricSTVSSIM(FILE ** streams, FILE * ref, int files_count, PictureData * frame, double ** results, int *& frames);

double countSTVSSIM(unsigned char ** datain1, unsigned char ** datain2, int size, int width);
void shiftData(unsigned char ** data, int size);
vector countARPS(unsigned char * block, unsigned char * framePrev, int x, int y, int width, int height, int T);
double countDelta(unsigned char*** filter, unsigned char*** cube1, unsigned char*** cube2, double muX, double muY);
double countDeltaSqr(unsigned char*** filter, unsigned char*** cube, double mu);
double countMu(unsigned char*** filter, unsigned char*** cube);
double countSSIM3D(unsigned char *** filter, unsigned char ***  cube1, unsigned char ***  cube2);
unsigned char *** generateCube();
void fillCube(unsigned char ** datain, int pos, unsigned char *** out, int width);

unsigned char **** generateFilters();

int countSAD(unsigned char * rect1, unsigned  char * rect2);