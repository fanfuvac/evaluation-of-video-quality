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
PictureData *getVideoInfo(string path){
	PictureData * data=new PictureData;

	string cmd="ffprobe -v error -of flat=s=_ -select_streams v:0 -show_entries stream=width,height -of default=noprint_wrappers=1:nokey=1 "+path;
	FILE *stream = popen(cmd.c_str(), "r");
	char buffer [50];
	fgets(buffer, 10, stream);
	data->width = atoi (buffer);
	fgets(buffer, 10, stream);
	data->height = atoi (buffer);
	return data;
}
FILE * startFfmepg(string path){
	
    string cmd="ffmpeg -i "+path+" -f image2pipe -pix_fmt rgb24 -vcodec rawvideo - 2>/dev/null";
    
    FILE *stream = popen(cmd.c_str(), "rb");
    if (stream){
      
    }
	return stream;
}
int main(){
	string file1 = "../data/W1.mkv";
	string file2 = "../data/W6.mp4";
	/*BMPData frame=getData("frame.bmp");
	BMPData frame2=getData("frame_2.bmp");*/
	int frame_count=7192;
	float * results=new float[frame_count];
	const int MAX_BUFFER = 2048000;
	
	PictureData * frame=getVideoInfo(file1);
	frame->data = new char[frame->width*frame->height*3];
	
	PictureData * frame2=getVideoInfo(file2);
	frame2->data = new char[frame->width*frame->height*3];
	
	FILE * stream=startFfmepg(file1);
	FILE * stream2=startFfmepg(file2);
	for (int j=1;j<=frame_count;j++){
		
			 // memcpy(&c, buffer, bytesize); //copy the data to a more usable structure for bit manipulation later

            if(fread(frame->data,frame->width*frame->height*3 , 1,stream)== NULL)
            {
               cout<<"error"<<endl;
            }
			if (fread(frame2->data,frame->width*frame->height*3 ,1, stream2) == NULL)
            {
               cout<<"error2"<<endl;
            }
       //pclose(stream);
		/*sprintf (buffer, "../data/Pictures%d.png", j);
		PictureData * frame=read_png_file(buffer);
		sprintf (buffer, "../data4/Pictures%03d.png", j);
		PictureData * frame2=read_png_file(buffer);*/
		float sum=0;
		int pixel;
		for(int i=0;i<frame->width*frame->height*3;i++){
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
		cout<<j<<" "<<psnr<<endl;
		
		//unsigned char a=frame->data[198];
		/*cout<<"SUM: "<<sum<<endl;
		cout<<"MSE: "<<mse<<endl;
		cout<<"MAX^2: "<<max*max<<endl;
		cout<<"PSNR: "<<10*log10(max*max/mse)<<endl;*/
	}
	delete frame->data;
		delete frame2->data;
		delete frame;
		delete frame2;
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