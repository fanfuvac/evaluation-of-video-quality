all:
	g++ PSNR.cpp png_decode.cpp -lpng -g -o PSNR
	