#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <math.h>     
using namespace std;

#define header_size 14+40
struct BMPData{
	int size;
	char* data;
	int width;
	int height;
};

//gets array of data from BMP file, skipping the header
BMPData getData(string filename){
	ifstream frame;
	frame.open (filename.c_str(),ios::in|ios::binary);
	//frame.seekg ( 14+4 );
	BMPData data;
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

int main(){
	BMPData frame=getData("frame.bmp");
	BMPData frame2=getData("frame_2.bmp");
	float sum=0;
	int pixel;
	for(int i=0;i<frame.size;i++){
		if(i%1==0){
			pixel=pow((unsigned char)frame.data[i]-(unsigned char)frame2.data[i],2);
		//cout<<int((unsigned char)frame.data[i])<<endl;
			sum+=pixel;
		}
	}
	float max=pow(2,8)-1;
	
	float mse=sum/frame.width/frame.height/3;
	//unsigned char a=frame.data[198];
	cout<<"SUM: "<<sum<<endl;
	cout<<"MSE: "<<mse<<endl;
	cout<<"MAX^2: "<<max*max<<endl;
	cout<<"PSNR: "<<10*log10(max*max/mse)<<endl;
}