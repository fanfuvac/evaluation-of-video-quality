#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <math.h>  
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include "stvssim.h"
#include "SSIM.h"
#include <stdlib.h>     /* abs */
#include <map>

using namespace std;
double countSTVSSIM(unsigned char ** datain1, unsigned char ** datain2, int size, int width) {
	unsigned char * out = new unsigned char[RECT_SIZE];
	
	int T = ROOD_SIZE;
	vector vct;
	vct.x = 0;
	vct.y = 0;
	int filter;
	unsigned char **** filters = generateFilters();
	double * tmpRes = new double[size];
	unsigned char *** cube1 = generateCube();
	unsigned char *** cube2 = generateCube();
	int k = 0;
	for (int i = 0; i < size / width - RECT_SQRT; i += SKIP_SIZE) {
		for (int j = 0; j < width - RECT_SQRT; j += SKIP_SIZE, k++) {
			getRect(datain1[FRAME_CNT / 2], i, width, out);
			//if (abs(vct.x) > abs(vct.y)) T = abs(vct.x); FIXME
			//if (abs(vct.x) < abs(vct.y)) T = abs(vct.y);
			vct = countARPS(out, datain1[FRAME_CNT / 2 - 1], j, i, width, size / width, T);

			if ((vct.x > vct.y * 2 && vct.x*-1 < 2 * vct.y) || (vct.x < vct.y * 2 && vct.x*-1 > 2 * vct.y)) { //y=0
				filter = 0;
			}
			else if ((vct.y > vct.x * 2 && vct.y > -2 * vct.x) || (vct.y < vct.x * 2 && vct.y < -2 * vct.x)) { //x=0
				filter = 2;
			}
			else if ((vct.y > vct.x / 2 && vct.y < 2 * vct.x) || (vct.y < vct.x / 2 && vct.y > 2 * vct.x)) { //y=x
				filter = 3;
			}
			else if ((vct.y > vct.x / -2 && vct.y < -2 * vct.x) || (vct.y <vct.x / -2 && vct.y>-2 * vct.x)) { //y=-x
				filter = 1;
			}
			else if (vct.x == 0 && vct.y == 0) {
				filter = 8;
			}
			else if (vct.x == vct.y * 2) { //exactly between 2 axes
				filter = 4;
			}
			else if (vct.x == vct.y * -2) {
				filter = 5;
			}
			else if (-2 * vct.x == vct.y) {
				filter = 6;
			}
			else if (2 * vct.x == vct.y) {
				filter = 7;
			}
			else {
				cout << "WUT - nonsense vector " << vct.x << " " << vct.y << endl;
			}



			//3D-SSIM part
			fillCube(datain1, i*width + j, cube1,width);
			fillCube(datain2, i*width + j, cube2, width);
			if (filter < 4) {
				tmpRes[k] = countSSIM3D(filters[filter], cube1, cube2);
				//cout << tmpRes[k] << endl;
			}
			else if (filter < 8) {
				double a, b;
				switch (filter) {
				case 4:
					a = countSSIM3D(filters[0], cube1, cube2);
					b = countSSIM3D(filters[3], cube1, cube2);
					break;
				case 5:
					a = countSSIM3D(filters[0], cube1, cube2);
					b = countSSIM3D(filters[1], cube1, cube2);
					break;
				case 6:
					a = countSSIM3D(filters[1], cube1, cube2);
					b = countSSIM3D(filters[2], cube1, cube2);
					break;
				case 7:
					a = countSSIM3D(filters[2], cube1, cube2);
					b = countSSIM3D(filters[3], cube1, cube2);
					break;
				}
				tmpRes[k] = (a + b) / 2;
				//cout << tmpRes[k] << endl;
			}
			else {
				double a, b, c, d;
				a = countSSIM3D(filters[0], cube1, cube2);
				b = countSSIM3D(filters[1], cube1, cube2);
				c = countSSIM3D(filters[2], cube1, cube2);
				d = countSSIM3D(filters[3], cube1, cube2);
				tmpRes[k] = (a + b+c+d) / 4;
				//cout << tmpRes[k] << endl;
			}
			if (tmpRes[k] > 1) {
				cout<< tmpRes[k] <<endl;
			}
		}
	}
	double res = countRes(tmpRes, k);
	delete[] tmpRes;
	delete[] out;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < RECT_SQRT_3D; j++) {
			for (int k = 0; k < RECT_SQRT_3D; k++) {
				delete filters[i][j][k];
				
			}
			delete[] cube1[i][j];
			delete[] cube2[i][j];
			delete[] filters[i][j];
		}
	}

	return res;

}
double countSSIM3D(unsigned char *** filter, unsigned char ***  cube1, unsigned char ***  cube2) {
	double muX = countMu(filter,cube1);
	double muY = countMu(filter, cube2);

	double deltaSqrX = countDeltaSqr(filter, cube1, muX);
	double deltaSqrY = countDeltaSqr(filter, cube2, muY); 

	double delta = countDelta(filter, cube1,cube2, muX, muY);
	
	double ssim3D = ((2 * muX*muY + C1)*(2 * delta + C2)) / ((muX*muX + muY*muY + C1)*(deltaSqrX + deltaSqrY + C2));
	return ssim3D;
}
double countMu(unsigned char*** filter,unsigned char*** cube){
	double res = 0;
	int res2=0;
	for (int alpha = 0; alpha < RECT_SQRT_3D; alpha++) {
		for (int beta = 0; beta < RECT_SQRT_3D; beta++) {
			for (int gamma = 0; gamma < FRAME_CNT; gamma++) {
				res2 += filter[gamma][alpha][beta];
				res += filter[gamma][alpha][beta] * cube[gamma][alpha][beta];
			}
		}
	}
	return res/(RECT_SQRT_3D*FRAME_CNT);
}

double countDeltaSqr(unsigned char*** filter, unsigned char*** cube,double mu) {
	double res = 0;
	for (int alpha = 0; alpha < RECT_SQRT_3D; alpha++) {
		for (int beta = 0; beta < RECT_SQRT_3D; beta++) {
			for (int gamma = 0; gamma < FRAME_CNT; gamma++) {
				res += filter[gamma][alpha][beta] * (cube[gamma][alpha][beta]-mu)*(cube[gamma][alpha][beta] - mu);
			}
		}
	}
	return res/ (RECT_SQRT_3D*FRAME_CNT);
}

double countDelta(unsigned char*** filter, unsigned char*** cube1, unsigned char*** cube2, double muX, double muY) {
	double res = 0;
	for (int alpha = 0; alpha < RECT_SQRT_3D; alpha++) {
		for (int beta = 0; beta < RECT_SQRT_3D; beta++) {
			for (int gamma = 0; gamma < FRAME_CNT; gamma++) {
				res += filter[gamma][alpha][beta] * (cube1[gamma][alpha][beta] - muX)*(cube2[gamma][alpha][beta] - muX);
			}
		}
	}
	return res/(RECT_SQRT_3D*FRAME_CNT);
}

//Generates 3D array used for SSIM 3D
unsigned char *** generateCube(){
	unsigned char *** out = new unsigned char**[FRAME_CNT];
	for (int i = 0; i < FRAME_CNT; i++) {
		out[i] = new unsigned char*[RECT_SQRT_3D];
		for (int j = 0; j < RECT_SQRT_3D; j++) {
			out[i][j] = new unsigned char[RECT_SQRT_3D];
		}
	}

	return out;
}

//Fill 3D array with data of surroundings pixels
void fillCube(unsigned char ** datain, int pos, unsigned char *** out, int width) {
	for (int i = 0; i < FRAME_CNT; i++) {
		for (int j = 0; j < RECT_SQRT_3D; j++) {
			memcpy(out[i][j], datain[i] + pos + width*j, RECT_SQRT_3D);
		}
	}
}
unsigned char **** generateFilters() {
	unsigned char **** out = new unsigned char***[4];
	for (int i = 0; i < 4; i++) {
		out[i]=new unsigned char**[FRAME_CNT];
		for (int j = 0; j < FRAME_CNT; j++) {
			out[i][j] = new unsigned char*[RECT_SQRT_3D];
			for (int k = 0; k < RECT_SQRT_3D; k++) {
				out[i][j][k] = new unsigned char[RECT_SQRT_3D];
			}
		}
	}
	for (int i = 0; i < FRAME_CNT; i++) { //FIXME modify to create 2D array and memcpy
		for (int j = 0; j < RECT_SQRT_3D; j++) {
			for (int k = 0; k < RECT_SQRT_3D; k++) {
				if (j  == RECT_SQRT_3D / 2) { //horizontal filter
					out[0][i][j][k] = 1;
				}
				else {
					out[0][i][j][k] = 0;
				}
				if (k == RECT_SQRT_3D / 2) { //vertical filter
					out[1][i][j][k] = 1;
				}
				else {
					out[1][i][j][k] = 0;
				}


				if (j == k) { //x=y filter
					out[2][i][j][k] = 1;
				}
				else {
					out[2][i][j][k] = 0;
				}
				if (k+j == RECT_SQRT_3D ) { //x=-y filter
					out[3][i][j][k] = 1;
				}
				else {
					out[3][i][j][k] = 0;
				}
			}
		}
	}
	return out;
}

vector countARPS(unsigned char * block, unsigned char * framePrev, int x,int y,int width, int height, int T) {
	unsigned char * out = new unsigned char[RECT_SIZE];
	getRect(framePrev, x*y, width, out);
	int sad=countSAD(block, out);
	vector vOut;
	if (sad < ZERO_MVMT) {
		vOut.x = 0;
		vOut.y = 0;
		return vOut;
	}
	map<int , int > past;
	int xOrig = x;
	int yOrig = y;
	int res[5];
	while (1) {
		getRect(framePrev, x+y*width, width, out);
		res[0] = countSAD(block, out);
		if (x - T - RECT_SQRT / 2 > 0) {
			getRect(framePrev, (x-T)+y*width, width, out);
			res[3] = countSAD(block, out);
		}
		else res[3] = INT_MAX;
		if (x + T + RECT_SQRT / 2< width) {
			getRect(framePrev, (x + T)+y*width, width, out);
			res[1] = countSAD(block, out);
		}
		else res[1] = INT_MAX;
		if (y + T + RECT_SQRT/2< height) {
			getRect(framePrev, x+(y+T)*width, width, out);
			res[2] = countSAD(block, out);
		}
		else res[2] = INT_MAX;
		if (y - T - RECT_SQRT / 2 > 0) {
			getRect(framePrev, x+(y-T)*width, width, out);
			res[4] = countSAD(block, out);
		}
		else res[4] = INT_MAX;
		int min = INT_MAX;
		int minPos = 0;
		for (int i = 0; i < 5; i++) {
			if (res[i] < min) {
				min = res[i];
				minPos = i;
			}
		}
		if (minPos==0) {
			vOut.x = x - xOrig;
			vOut.y = y - yOrig;
			return vOut;
		}
		switch (minPos) {
		case 1: 
			x += T;
			break;
		case 2:	
			y += T;
			break;
		case 3: 
			x -= T;
			break;
		case 4:	
			y -= T;
			break;
		}	
	}
}

int countSAD(unsigned char * rect1, unsigned  char * rect2) {
	int sad = 0;
	for (int i = 0; i < RECT_SIZE; i++) {
		sad += abs(rect1[i] - rect2[i]);
	}
	return sad;
}