

	
all: 	yuv2tif  tif2yuv 
      
yuv2tif : yuv2tif.cpp
	g++ -O3  yuv2tif.cpp -o yuv2tif -ltiff
	
tif2yuv : tif2yuv.cpp
	g++ -O3  tif2yuv.cpp -o tif2yuv -ltiff
	
	
	
clean : 
	rm -v yuv2tif  tif2yuv	
  
