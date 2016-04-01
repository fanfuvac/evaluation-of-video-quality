#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <math.h>  
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include "png_decode.h"   
using namespace std;

#define header_size 14+40


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
}
int compare (const void * a, const void * b)
{
  return ( *(float*)a - *(float*)b );
}
int main(){
	/*BMPData frame=getData("frame.bmp");
	BMPData frame2=getData("frame_2.bmp");*/
	int frame_count=452;
	float * results=new float[frame_count];
	char buffer [50];
	for (int j=1;j<=frame_count;j++){
		sprintf (buffer, "../data/Pictures%d.png", j);
		PictureData * frame=read_png_file(buffer);
		sprintf (buffer, "../data2/Pictures%d.png", j);
		PictureData * frame2=read_png_file(buffer);
		float sum=0;
		int pixel;
		for(int i=0;i<frame->width*frame->height;i++){
			if(i%1==0){
				pixel=pow((unsigned char)frame->data[i]-(unsigned char)frame2->data[i],2);
				//cout<<int((unsigned char)frame->data[i])<<endl;
				sum+=pixel;
			}
		}
		float max=pow(2,8)-1;
		
		float mse=sum/frame->width/frame->height/3;
		float psnr=10*log10(max*max/mse);
		results[j]=psnr;
		cout<<psnr<<endl;
		delete frame->data;
		delete frame2->data;
		delete frame;
		delete frame2;
		//unsigned char a=frame->data[198];
		/*cout<<"SUM: "<<sum<<endl;
		cout<<"MSE: "<<mse<<endl;
		cout<<"MAX^2: "<<max*max<<endl;
		cout<<"PSNR: "<<10*log10(max*max/mse)<<endl;*/
	}
	float sum=0;
	int frames=frame_count;
	for (int i=0;i<frame_count;i++){
		if(isfinite(results[i]))
			sum+=results[i];
		else frames--;
	}
	cout<<"AVG: "<<sum/frames<<endl;
	 qsort (results, frames, sizeof(float), compare);
	 cout<<"Median: "<<results[frames/2]<<endl;
}