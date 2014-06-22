// set up for 12 bits 


#include <tiffio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
using namespace std;

void Subsample444to420(unsigned short ** src, unsigned short ** dst, short width, short height, short algorithmn, unsigned long minCV, unsigned long maxCV);

main(int argc, char* argv[])
{
    TIFF* tif = TIFFOpen(argv[1], "r");
    if (tif) {
        tdata_t buf1;
        tdata_t buf2;
        tdata_t buf3;
        tdata_t buf4;
        tstrip_t strip;
        uint32* bc;
        uint32 stripsize;
        uint32 pixel;
        uint32 pixelStart;
        uint32 ib;
        short  numStrips;
        short  stripStart;
	     ofstream yuvOut;
	     short D709 = 0;
	     short D2020 = 0;
	     float tmpF = 0.0;
	     short HD = 0; // ==1 if 1920x1080 cutout
	     short qHD = 0; // ==1 if 960x540 cutout
	     short Y500 =0;
	     short Y100 =0;
	     short DXYZ = 1;
	     // Bit Depth
	     short BD = 12; // 12 bit mode
	     
	     // Defaults to process as 16 bit
	     unsigned long Half = 32768; // e.g half value 12 bits would have been 2048
	     unsigned long minCV = 0; //just full range short for use later
	     unsigned long maxCV = 65535;
	     short FIR = 1;
	     short ALPHA = 0; // no alpha channel is default
	     short numChan2 = 3*2; // only 3 channels rgb with no alpha
	     short numChan = 3;
	     short FULLRANGE = 0; // default is video range
	     unsigned short D = 64; //D=4 for 12 bits =1 10bits = 16 for 14 bits
	     unsigned short minVR, maxVR, minVRC,maxVRC;
	     // video range is
	     // Luma and R,G,B:  CV = Floor(876*D*N+64*D+0.5)
	     // Chroma:  CV = Floor(896*D*N+64*D+0.5)
	     // Process will be assume range of input is correct and clip/ignore out of range on input
	     // convert to YUV working in integer or floating point
	     // subsample at 16 bit
	     // range limit subsampled output to 16 bit video legal range
	     // shift down to 10,12,14 bits and write out into YUV planar file
	     
	     //Process Args:
	     short arg = 2;
	     while(arg < argc) {
	     	
		     if(strcmp(argv[arg],"BOX")==0)FIR = 0;
		     if(strcmp(argv[arg],"FULL")==0)FULLRANGE = 1;
		     if(strcmp(argv[arg],"ALPHA")==0){
		     		ALPHA = 1;
		     		numChan2= 4*2;
		     		numChan = 4;
		     }		     
		     		     	     	
		     if(strcmp(argv[arg],"709")==0)D709 = 1;
		     if(strcmp(argv[arg],"2020")==0)D2020 = 1;
		     
		     if(strcmp(argv[arg],"HD1920")==0)HD = 1;
		     if(strcmp(argv[arg],"HD960")==0)qHD = 1;
		     
		     if(strcmp(argv[arg],"B10")==0) {
		     	  BD = 10;
        	     printf("\nprocessing line data 10 bits\n");		     	  
		     }		     

		     if(strcmp(argv[arg],"B14")==0) {
		     	  BD = 14;
        	     printf("\nprocessing line data 14 bits\n");		     	  
		     }	
		     
		     if(strcmp(argv[arg],"Y100")==0)Y100 = 1;
		     if(strcmp(argv[arg],"Y500")==0)Y500 = 1;	
		     
		     if(D709 || D2020 || Y100 || Y500)DXYZ = 0;	     
		     
		     if(strcmp(argv[arg],"-h")==0) {
		     	 printf("\n ARGS:\n 709 (use Rec709 Color Dif)\n 2020 (use Rec2020 Color Dif)\n HD1920 (cut out center 1920x1080 images)\n HD960 (cut out center 960x540 images)\n (no args get Y'DzDx color difference and 3840x2160 cutout)\n\n\n");
		     	 exit(0);
		     }
	     

	     arg++;
	     }
	     if(D709)printf("Processing for Rec709\n");
	     if(D2020)printf("Processing for Rec2020\n");
	     if(HD)printf("Cutting out HD1920x1080\n");
	     if(qHD)printf("Cutting out HD960x540\n");
	     if(ALPHA)printf("Assuming TIF has ALPHA channel\n");	  	     	  
	     
        // Set up for video range if needed
        if(!FULLRANGE) {
        	  printf("Processing for Video Range\n");
           minVR = 64*D;
           maxVR = 876*D+minVR;
           minVRC = minVR;
           maxVRC = 896*D+minVRC;                
           //achromatic point for chroma will be "Half"(e.g. 512 for 10 bits, 2048 for 12 bits etc..)
        
        }        



	     // max pixel value (e.g 12 bits is 4095)
	     
	     // set up alternate differencing equations for Y'DzwDxw
	     float P,Q,RR,S;
	     if (Y100) {
		     printf("Setting PQRS for 100nit PQ color difference point\n");
		     P = -0.5;
		     Q = 0.491722;
		     RR = 0.5;
		     S = -0.49495;  
	     } else if(Y500) {
		     printf("Setting PQRS for 500nit PQ color difference point\n");
		     P = -0.5;
		     Q = 0.493393;
		     RR = 0.5;
		     S = -0.49602;	     	
	     }
        //printf("Stripsize: %d\n",TIFFStripSize(tif));
        printf("NumStrips: %d\n",TIFFNumberOfStrips(tif));
        numStrips = TIFFNumberOfStrips(tif);
        //printf("TIFFTAG_ROWSPERSTRIP: %d\n",TIFFGetField(tif,TIFFTAG_ROWSPERSTRIP));       
        //printf("SamplePerPixel: %d\n",TIFFGetField(tif,TIFFTAG_SAMPLESPERPIXEL)); 
        //printf("TIFFTAG_IMAGELENGTH = %d\n",TIFFGetField(tif,TIFFTAG_IMAGELENGTH));
        

        
      
        //get count of bytes in each strip of image as vector of rows/lines
        TIFFGetField(tif, TIFFTAG_STRIPBYTECOUNTS, &bc);
        stripsize = bc[0];
        printf("Original: Stripsize pixels: %d, Stripsize bytes: %d\n",stripsize/numChan2, stripsize);
        buf1 = _TIFFmalloc(stripsize);
        buf2 = _TIFFmalloc(stripsize);
        buf3 = _TIFFmalloc(stripsize);
        buf4 = _TIFFmalloc(stripsize);
        
        // Adjust numStrips to 2160
        if(numStrips != 2160 && numStrips != 1080)
        {
        	  printf("Error numstrips not 2160 or 1080\n");
        	  exit(0);
        }
        
        stripStart = 0;
        // Adjust stripstart to cut out if needed
        if(HD) stripStart = (numStrips - 1080)/2;
        if(qHD)stripStart = (numStrips - 540)/2;
        
        
        // Adjust stripsize to 3840x4x2bytes 3840xnumChanx2
        // pixelstart will be byte location of first pixel
        // stripsize will be byte location just after last pixel 
        if(stripsize>(960*numChan2))
        {
        	   pixelStart = (stripsize - 3840*numChan2)/2;
        	   
				// if HD cutout adjust:
				if(HD) pixelStart = (stripsize - 1920*numChan2)/2;
				// if qHD cutout adjust:
				if(qHD)pixelStart = (stripsize - 960*numChan2)/2;        	   
        	   
        	   
        	   // stripsize now to mean strip end pixel byte
        	   // scanline is pixelStart,stripsize,pixelstart e.g. = 8*4196
        	   stripsize = stripsize - pixelStart;
        	   printf("    pixelStart: %d  stripsize: %d  (pixels)\n",pixelStart/numChan2, stripsize/numChan2 - pixelStart/numChan2);
        	   printf("    pixelStart: %d  stripsize: %d  (bytes)\n",pixelStart, stripsize - pixelStart);
        }
        
        // Read 16 bit tiff values into unsigned short into unsigned int for headroom
        // assumption is that RGB = X'Y'Z' from ODT output into 16 bit tiff (as 12 bit clamped values)
        // other tests will clamp ODT at 14 bits and maybe 10 bits
        uint32 R;
        uint32 G;
        uint32 B;
        uint32 A;
        // intermediate values for color differencing
        uint32 Y;
        long Dz;
        long Dx,RC,BC;
        // intermediate 16 bit values for final Y' and subsamples color difference
        unsigned short Y4[4][4];
        int Dz4[4][4];
        int Dx4[4][4];
        // subsampled
        int Dz0[2][2];
        int Dx0[2][2];
        // arrar to store line of output for writing process
        unsigned short  *Line;
        
        unsigned short** YP;
        int arraySizeX = stripsize/numChan2 - pixelStart/numChan2; // eg 3840
        int arraySizeY = numStrips-2*stripStart;
        printf("Frame size: %d x %d\n",arraySizeX, arraySizeY);
        int arraySizeXH = arraySizeX/2; // eg 1920 cols
        int arraySizeYH = arraySizeY/2; // eg 1080 rows
	     YP = (unsigned short**) malloc(arraySizeX*sizeof(unsigned short*)); // cols
        for (int i = 0; i < arraySizeX; i++)
            YP[i] = (unsigned short*) malloc(arraySizeY*sizeof(unsigned short)); //rows

		  // Source 444 chroma arrays            
        unsigned short** Cb444;
        Cb444 = (unsigned short**) malloc(arraySizeX*sizeof(unsigned short*)); // cols
        for (int i = 0; i < arraySizeX; i++)
            Cb444[i] = (unsigned short*) malloc(arraySizeY*sizeof(unsigned short)); //rows            

        unsigned short** Cr444;
        Cr444 = (unsigned short**) malloc(arraySizeX*sizeof(unsigned short*)); // cols
        for (int i = 0; i < arraySizeX; i++)
            Cr444[i] = (unsigned short*) malloc(arraySizeY*sizeof(unsigned short)); //rows      
            
        //Destination 420 chroma arrays
        unsigned short** DzP;
	     DzP = (unsigned short**) malloc(arraySizeXH*sizeof(unsigned short*)); // cols
        for (int i = 0; i < arraySizeXH; i++)
            DzP[i] = (unsigned short*) malloc(arraySizeYH*sizeof(unsigned short)); // rows
            
        unsigned short** DxP;
	     DxP = (unsigned short**) malloc(arraySizeXH*sizeof(unsigned short*)); // cols
        for (int i = 0; i < arraySizeXH; i++)
            DxP[i] = (unsigned short*) malloc(arraySizeYH*sizeof(unsigned short)); // rows
       
       // set output line array to unsigned short              
       Line =  (unsigned short*) malloc(arraySizeX*sizeof(unsigned short));

		  // Open Binary PlanarYUV file for writing:
	     yuvOut.open("YDzDx.yuv", ios::ate | ios::app | ios::out | ios::binary);
	     printf("Opened YDzDx.yuv to file end and appending:\n");
	     
	     // Reduce numStrips to what is needed to center cut:
	     numStrips = numStrips - stripStart;
    
        // read lines to cover arraySizeY number of lines (e.g. 2160)
        for (strip = stripStart; strip < numStrips; strip+=4) {
        	
            // strip is equal to the scanline from 0 to 2159
            
            
            
            // Read 4x4 section into planar yuv420
            TIFFReadRawStrip(tif, strip,   buf1, bc[0]);
            TIFFReadRawStrip(tif, strip+1, buf2, bc[0]);
            TIFFReadRawStrip(tif, strip+2, buf3, bc[0]);
            TIFFReadRawStrip(tif, strip+3, buf4, bc[0]);
                  
            //printf("\rReading Line: %d-%d, Bytes: %d",strip+1,strip+4, bc[0]);      
            	
            // chunk scan line into 4 RGBA or RGB 16bit shorts
            // process blocks of 4 lines by 4 pixels blocks 
            // so skip 16 values per block or 32 bytes
            // 32768 bytes is 4096 * 4(rgba) * 16bits = 4096 * 4 * 2
            for (pixel = pixelStart; pixel < stripsize; pixel+=(numChan2*4)) {
            // step line in blocks of numChan pixels:
            // ib is the pixel counter from 1-4 such 
            // that numChan2*ib is byte location of pixel start
            for (ib = 0;ib<4;ib+=1) {              
            
              // read rgba
              // 12 bit YDzDx formular from octave script:
              //  trans_new = [ 0 1 0; 0 -0.5 0.5; 0.5 -0.5 0 ]
              //  output_new = trans_new * source_vector + [0; 2048; 2048];
              //  2048 = 2**12 / 2
              // Y = G = Y
              // Dz = -0.5*G + 0.5*B  = 0.5*Z - 0.5*Y + 2048
              // Dx = 0.5*R -0.5*G  = 0.5*X - 0.5*Y  +2048
              
              // [consider changing RGBA aka XYZ to 12 bits here first]
              // pixel/2 is unsigned short width index to first R 
              // in block of 4 pixels of 4 values RGBA
              // buf1-4 are the 4 scan lines reading
              // ib=0-3 * 4 is pixel (unsigned short width) offset to pixel value

              R = static_cast<unsigned short *>(buf1)[pixel/2+ib*numChan];
              G = static_cast<unsigned short *>(buf1)[pixel/2+1+ib*numChan];
              B = static_cast<unsigned short *>(buf1)[pixel/2+2+ib*numChan];
              if(ALPHA) A = static_cast<unsigned short *>(buf1)[pixel/2+3+ib*numChan];          
              
              // Insure RGB are clipped to video range if required:
              if(!FULLRANGE) {
              R = (R<minVR) ? minVR : R;
              G = (G<minVR) ? minVR : G;
              B = (B<minVR) ? minVR : B;
              
              R = (R>maxVR) ? maxVR : R;
              G = (G>maxVR) ? maxVR : G;
              B = (B>maxVR) ? maxVR : B;
                            
              }

          
 				 if(DXYZ) {	

	              Y = G;
	              Dz = (int)(-((float)G)/2.0 +((float)B)/2.0 +0.5); //-(int)((G+1)>>1) + (int)((B+1)>>1);
	              Dx =  (int)(-((float)G)/2.0 +((float)R)/2.0 +0.5); //(int)((R+1)>>1) - (int)((G+1)>>1); 
	              //printf("Y: %d, Dz: %d, Dz: %d\n",Y,Dz,Dx);
	              
              } else if(D2020){
              	  tmpF = (0.2627*(float)R + 0.6780*(float)G + 0.0593*(float)B) +0.5;
              	  Y = (uint32)(tmpF);
              	  Dz = (int)((((float)B) - tmpF)/1.8814 + 0.5);
              	  Dx = (int)((((float)R) - tmpF)/1.4746 + 0.5);
              } else if(D709) {
              	  tmpF = (0.2126*(float)R + 0.7152*(float)G + 0.0722*(float)B) +0.5;
              	  Y = (uint32)(tmpF);
              	  Dz = (int)((((float)B) - tmpF)/1.8556 + 0.5);
              	  Dx = (int)((((float)R) - tmpF)/1.5748 + 0.5);              	
              } else if(Y100 || Y500) {
              	  Y = (uint32)(G);
              	  Dz = (int)( P *((float)G) + Q *((float)B) + 0.5);
              	  Dx = (int)( RR *((float)R) + S *((float)G) + 0.5);
              } else {
              	 printf("Can't determine color difference to use?\n\n");
              	 exit(0);
              }

              Dz = Dz + Half - 1;
              Dx = Dx + Half - 1;
              //printf("Y: %d, Dz: %d, Dz: %d\n",Y,Dz,Dx);
              // clamp to full range                         
              if(Y>maxCV)Y=maxCV;
              if(Dz>maxCV)Dz=maxCV;
              if(Dz<minCV)Dz=minCV;
              if(Dx>maxCV)Dx=maxCV;
              if(Dx<minCV)Dx=minCV;
 
              if(!FULLRANGE) {
              Y= (Y<minVR) ? minVR : Y;
              Dz = (Dz<minVRC) ? minVR : Dz;
              Dx = (Dx<minVRC) ? minVR : Dx;
              
              Y = (Y>maxVR) ? maxVR : Y;
              Dz = (Dz>maxVRC) ? maxVRC : Dz;
              Dx = (Dx>maxVRC) ? maxVRC : Dx;
              }
              //printf("Y: %d, Dz: %d, Dz: %d\n",Y,Dz,Dx);
              Y4[0][ib] = Y;
              Dz4[0][ib] = Dz;
              Dx4[0][ib] = Dx;
              
              // Repeat for all 4 rows read
              R = static_cast<unsigned short *>(buf2)[pixel/2+ib*numChan];
              G = static_cast<unsigned short *>(buf2)[pixel/2+1+ib*numChan];
              B = static_cast<unsigned short *>(buf2)[pixel/2+2+ib*numChan];
              if(ALPHA) A = static_cast<unsigned short *>(buf2)[pixel/2+3+ib*numChan];          

              // Insure RGB are clipped to video range if required:
              if(!FULLRANGE) {
              R = (R<minVR) ? minVR : R;
              G = (G<minVR) ? minVR : G;
              B = (B<minVR) ? minVR : B;
              
              R = (R>maxVR) ? maxVR : R;
              G = (G>maxVR) ? maxVR : G;
              B = (B>maxVR) ? maxVR : B;
                            
              }

          
 				 if(DXYZ) {	
	              Y = G;
	              Dz = (int)(-((float)G)/2.0 +((float)B)/2.0 +0.5); //-(int)((G+1)>>1) + (int)((B+1)>>1);
	              Dx =  (int)(-((float)G)/2.0 +((float)R)/2.0 +0.5); //(int)((R+1)>>1) - (int)((G+1)>>1); 
	              //printf("Y: %d, Dz: %d, Dz: %d\n",Y,Dz,Dx);
	              
              } else if(D2020){
              	  tmpF = (0.2627*(float)R + 0.6780*(float)G + 0.0593*(float)B) +0.5;
              	  Y = (uint32)(tmpF);
              	  Dz = (int)((((float)B) - tmpF)/1.8814 + 0.5);
              	  Dx = (int)((((float)R) - tmpF)/1.4746 + 0.5);
              } else if(D709) {
              	  tmpF = (0.2126*(float)R + 0.7152*(float)G + 0.0722*(float)B) +0.5;
              	  Y = (uint32)(tmpF);
              	  Dz = (int)((((float)B) - tmpF)/1.8556 + 0.5);
              	  Dx = (int)((((float)R) - tmpF)/1.5748 + 0.5);              	
              } else if(Y100 || Y500) {
              	  Y = (uint32)(G);
              	  Dz = (int)( P *((float)G) + Q *((float)B) + 0.5);
              	  Dx = (int)( RR *((float)R) + S *((float)G) + 0.5);
              } else {
              	 printf("Can't determine color difference to use?\n\n");
              	 exit(0);
              }

              Dz = Dz + Half - 1;
              Dx = Dx + Half - 1;
              // clamp to full range                         
              if(Y>maxCV)Y=maxCV;
              if(Dz>maxCV)Dz=maxCV;
              if(Dz<minCV)Dz=minCV;
              if(Dx>maxCV)Dx=maxCV;
              if(Dx<minCV)Dx=minCV;
 
 
              if(!FULLRANGE) {
              Y= (Y<minVR) ? minVR : Y;
              Dz = (Dz<minVRC) ? minVR : Dz;
              Dx = (Dx<minVRC) ? minVR : Dx;
              
              Y = (Y>maxVR) ? maxVR : Y;
              Dz = (Dz>maxVRC) ? maxVRC : Dz;
              Dx = (Dx>maxVRC) ? maxVRC : Dx;
              }
           
              Y4[1][ib] = Y;
              Dz4[1][ib] = Dz;
              Dx4[1][ib] = Dx;              

              R = static_cast<unsigned short *>(buf3)[pixel/2+ib*numChan];
              G = static_cast<unsigned short *>(buf3)[pixel/2+1+ib*numChan];
              B = static_cast<unsigned short *>(buf3)[pixel/2+2+ib*numChan];
              if(ALPHA) A = static_cast<unsigned short *>(buf3)[pixel/2+3+ib*numChan];         

              // Insure RGB are clipped to video range if required:
              if(!FULLRANGE) {
              R = (R<minVR) ? minVR : R;
              G = (G<minVR) ? minVR : G;
              B = (B<minVR) ? minVR : B;
              
              R = (R>maxVR) ? maxVR : R;
              G = (G>maxVR) ? maxVR : G;
              B = (B>maxVR) ? maxVR : B;
                            
              }


          
 				 if(DXYZ) {	

	              Y = G;
	              Dz = (int)(-((float)G)/2.0 +((float)B)/2.0 +0.5); //-(int)((G+1)>>1) + (int)((B+1)>>1);
	              Dx =  (int)(-((float)G)/2.0 +((float)R)/2.0 +0.5); //(int)((R+1)>>1) - (int)((G+1)>>1); 
	              //printf("Y: %d, Dz: %d, Dz: %d\n",Y,Dz,Dx);
	              
              } else if(D2020){
              	  tmpF = (0.2627*(float)R + 0.6780*(float)G + 0.0593*(float)B) +0.5;
              	  Y = (uint32)(tmpF);
              	  Dz = (int)((((float)B) - tmpF)/1.8814 + 0.5);
              	  Dx = (int)((((float)R) - tmpF)/1.4746 + 0.5);
              } else if(D709) {
              	  tmpF = (0.2126*(float)R + 0.7152*(float)G + 0.0722*(float)B) +0.5;
              	  Y = (uint32)(tmpF);
              	  Dz = (int)((((float)B) - tmpF)/1.8556 + 0.5);
              	  Dx = (int)((((float)R) - tmpF)/1.5748 + 0.5);              	
              } else if(Y100 || Y500) {
              	  Y = (uint32)(G);
              	  Dz = (int)( P *((float)G) + Q *((float)B) + 0.5);
              	  Dx = (int)( RR *((float)R) + S *((float)G) + 0.5);
              } else {
              	 printf("Can't determine color difference to use?\n\n");
              	 exit(0);
              }

              Dz = Dz + Half - 1;
              Dx = Dx + Half - 1;
              // clamp to full range                         
              if(Y>maxCV)Y=maxCV;
              if(Dz>maxCV)Dz=maxCV;
              if(Dz<minCV)Dz=minCV;
              if(Dx>maxCV)Dx=maxCV;
              if(Dx<minCV)Dx=minCV;

              if(!FULLRANGE) {
              Y= (Y<minVR) ? minVR : Y;
              Dz = (Dz<minVRC) ? minVR : Dz;
              Dx = (Dx<minVRC) ? minVR : Dx;
              
              Y = (Y>maxVR) ? maxVR : Y;
              Dz = (Dz>maxVRC) ? maxVRC : Dz;
              Dx = (Dx>maxVRC) ? maxVRC : Dx;
              }
             
              Y4[2][ib] = Y;
              Dz4[2][ib] = Dz;
              Dx4[2][ib] = Dx;                      
 
              R = static_cast<unsigned short *>(buf4)[pixel/2+ib*numChan];
              G = static_cast<unsigned short *>(buf4)[pixel/2+1+ib*numChan];
              B = static_cast<unsigned short *>(buf4)[pixel/2+2+ib*numChan];
              if(ALPHA) A = static_cast<unsigned short *>(buf4)[pixel/2+3+ib*numChan];   
              //printf("R = %d, G = %d, B = %d, A = %d\n",R,G,B,A);         

              // Insure RGB are clipped to video range if required:
              if(!FULLRANGE) {
              R = (R<minVR) ? minVR : R;
              G = (G<minVR) ? minVR : G;
              B = (B<minVR) ? minVR : B;
              
              R = (R>maxVR) ? maxVR : R;
              G = (G>maxVR) ? maxVR : G;
              B = (B>maxVR) ? maxVR : B;
                            
              }
                        
 				 if(DXYZ) {	

	              Y = G;
	              Dz = (int)(-((float)G)/2.0 +((float)B)/2.0 +0.5); //-(int)((G+1)>>1) + (int)((B+1)>>1);
	              Dx =  (int)(-((float)G)/2.0 +((float)R)/2.0 +0.5); //(int)((R+1)>>1) - (int)((G+1)>>1); 
	              //printf("Y: %d, Dz: %d, Dz: %d\n",Y,Dz,Dx);
	              
              } else if(D2020){
              	  tmpF = (0.2627*(float)R + 0.6780*(float)G + 0.0593*(float)B) +0.5;
              	  Y = (uint32)(tmpF);
              	  Dz = (int)((((float)B) - tmpF)/1.8814 + 0.5);
              	  Dx = (int)((((float)R) - tmpF)/1.4746 + 0.5);
              } else if(D709) {
              	  tmpF = (0.2126*(float)R + 0.7152*(float)G + 0.0722*(float)B) +0.5;
              	  Y = (uint32)(tmpF);
              	  Dz = (int)((((float)B) - tmpF)/1.8556 + 0.5);
              	  Dx = (int)((((float)R) - tmpF)/1.5748 + 0.5);              	
              } else if(Y100 || Y500) {
              	  Y = (uint32)(G);
              	  Dz = (int)( P *((float)G) + Q *((float)B) + 0.5);
              	  Dx = (int)( RR *((float)R) + S *((float)G) + 0.5);
              } else {
              	 printf("Can't determine color difference to use?\n\n");
              	 exit(0);
              }

              Dz = Dz + Half - 1;
              Dx = Dx + Half - 1;
              // clamp to full range                         
              if(Y>maxCV)Y=maxCV;
              if(Dz>maxCV)Dz=maxCV;
              if(Dz<minCV)Dz=minCV;
              if(Dx>maxCV)Dx=maxCV;
              if(Dx<minCV)Dx=minCV;

              if(!FULLRANGE) {
              Y= (Y<minVR) ? minVR : Y;
              Dz = (Dz<minVRC) ? minVR : Dz;
              Dx = (Dx<minVRC) ? minVR : Dx;
              
              Y = (Y>maxVR) ? maxVR : Y;
              Dz = (Dz>maxVRC) ? maxVRC : Dz;
              Dx = (Dx>maxVRC) ? maxVRC : Dx;
              }
              Y4[3][ib] = Y;
              Dz4[3][ib] = Dz;
              Dx4[3][ib] = Dx; 
              
              //printf("Y = %d, Dz = %d, Dx = %d\n",Y,Dz,Dx);
                          
            }   //end ib for loop 
            
            // Fill Source 444 Chroma Arrays Cb444 and Cr444
            
              short stripSave = strip;
              strip = strip - stripStart;
              YP[pixel/numChan2 - pixelStart/numChan2][strip] = Y4[0][0];
              YP[1+pixel/numChan2 - pixelStart/numChan2][strip] = Y4[0][1]; 
              YP[2+pixel/numChan2 - pixelStart/numChan2][strip] = Y4[0][2];
              YP[3+pixel/numChan2 - pixelStart/numChan2][strip] = Y4[0][3];
              YP[pixel/numChan2 - pixelStart/numChan2][strip+1] = Y4[1][0];
              YP[1+pixel/numChan2 - pixelStart/numChan2][strip+1] = Y4[1][1];
              YP[2+pixel/numChan2 - pixelStart/numChan2][strip+1] = Y4[1][2];
              YP[3+pixel/numChan2 - pixelStart/numChan2][strip+1] = Y4[1][3]; 
              YP[pixel/numChan2 - pixelStart/numChan2][strip+2] = Y4[2][0];
              YP[1+pixel/numChan2 - pixelStart/numChan2][strip+2] = Y4[2][1];
              YP[2+pixel/numChan2 - pixelStart/numChan2][strip+2] = Y4[2][2];
              YP[3+pixel/numChan2 - pixelStart/numChan2][strip+2] = Y4[2][3];
              YP[pixel/numChan2 - pixelStart/numChan2][strip+3] = Y4[3][0];
              YP[1+pixel/numChan2 - pixelStart/numChan2][strip+3] = Y4[3][1]; 
              YP[2+pixel/numChan2 - pixelStart/numChan2][strip+3] = Y4[3][2];
              YP[3+pixel/numChan2 - pixelStart/numChan2][strip+3] = Y4[3][3];
                                                  
              Cb444[pixel/numChan2 - pixelStart/numChan2][strip] = Dz4[0][0];
              Cb444[1+pixel/numChan2 - pixelStart/numChan2][strip] = Dz4[0][1]; 
              Cb444[2+pixel/numChan2 - pixelStart/numChan2][strip] = Dz4[0][2];
              Cb444[3+pixel/numChan2 - pixelStart/numChan2][strip] = Dz4[0][3];
              Cb444[pixel/numChan2 - pixelStart/numChan2][strip+1] = Dz4[1][0];
              Cb444[1+pixel/numChan2 - pixelStart/numChan2][strip+1] = Dz4[1][1];
              Cb444[2+pixel/numChan2 - pixelStart/numChan2][strip+1] = Dz4[1][2];
              Cb444[3+pixel/numChan2 - pixelStart/numChan2][strip+1] = Dz4[1][3]; 
              Cb444[pixel/numChan2 - pixelStart/numChan2][strip+2] = Dz4[2][0];
              Cb444[1+pixel/numChan2 - pixelStart/numChan2][strip+2] = Dz4[2][1];
              Cb444[2+pixel/numChan2 - pixelStart/numChan2][strip+2] = Dz4[2][2];
              Cb444[3+pixel/numChan2 - pixelStart/numChan2][strip+2] = Dz4[2][3];
              Cb444[pixel/numChan2 - pixelStart/numChan2][strip+3] = Dz4[3][0];
              Cb444[1+pixel/numChan2 - pixelStart/numChan2][strip+3] = Dz4[3][1]; 
              Cb444[2+pixel/numChan2 - pixelStart/numChan2][strip+3] = Dz4[3][2];
              Cb444[3+pixel/numChan2 - pixelStart/numChan2][strip+3] = Dz4[3][3];

              Cr444[pixel/numChan2 - pixelStart/numChan2][strip] = Dx4[0][0];
              Cr444[1+pixel/numChan2 - pixelStart/numChan2][strip] = Dx4[0][1]; 
              Cr444[2+pixel/numChan2 - pixelStart/numChan2][strip] = Dx4[0][2];
              Cr444[3+pixel/numChan2 - pixelStart/numChan2][strip] = Dx4[0][3];
              Cr444[pixel/numChan2 - pixelStart/numChan2][strip+1] = Dx4[1][0];
              Cr444[1+pixel/numChan2 - pixelStart/numChan2][strip+1] = Dx4[1][1];
              Cr444[2+pixel/numChan2 - pixelStart/numChan2][strip+1] = Dx4[1][2];
              Cr444[3+pixel/numChan2 - pixelStart/numChan2][strip+1] = Dx4[1][3]; 
              Cr444[pixel/numChan2 - pixelStart/numChan2][strip+2] = Dx4[2][0];
              Cr444[1+pixel/numChan2 - pixelStart/numChan2][strip+2] = Dx4[2][1];
              Cr444[2+pixel/numChan2 - pixelStart/numChan2][strip+2] = Dx4[2][2];
              Cr444[3+pixel/numChan2 - pixelStart/numChan2][strip+2] = Dx4[2][3];
              Cr444[pixel/numChan2 - pixelStart/numChan2][strip+3] = Dx4[3][0];
              Cr444[1+pixel/numChan2 - pixelStart/numChan2][strip+3] = Dx4[3][1]; 
              Cr444[2+pixel/numChan2 - pixelStart/numChan2][strip+3] = Dx4[3][2];
              Cr444[3+pixel/numChan2 - pixelStart/numChan2][strip+3] = Dx4[3][3];


              strip = stripSave;

                            
              
            }   // end of 4 scanline blocking
        } // end of loop of scanlines
        
        _TIFFfree(buf1);
        _TIFFfree(buf2);
        _TIFFfree(buf3);
        _TIFFfree(buf4);
        TIFFClose(tif);

    // Subsample =1 means use FIR filter
    
	  if(!FULLRANGE) {
		Subsample444to420(Cb444, DzP, arraySizeX, arraySizeY,FIR,minVRC,maxVRC);
		Subsample444to420(Cr444, DxP, arraySizeX, arraySizeY,FIR,minVRC,maxVRC);
	  } else {
	     Subsample444to420(Cb444, DzP, arraySizeX, arraySizeY,FIR,minCV,maxCV);
         Subsample444to420(Cr444, DxP, arraySizeX, arraySizeY,FIR,minCV,maxCV);
      }      
       // Write YUV Frame
       // write planer output:
        
        // Y
        if(BD==10) {
        	printf("writing line data 10 bits\n");
        } else if(BD==12) {
        	printf("writing line data 12 bits\n");
        } else if(BD==14) {
         printf("Writing line data 14 bits\n");
        }




        // limit video range at 16 bits
        // shift data down to desired bit depth
        for(unsigned short r = 0; r < arraySizeY;r++)
        {
          for(unsigned short c = 0; c < arraySizeX;c++)
          {
              // Fill line array 
              // (10 bits) Line[c] = (((YP[c][r]+2)>>2));
        // limit to video range (for 16 bits then shift later)
		        if(!FULLRANGE) {
			        YP[c][r] = (YP[c][r]<minVR) ? minVR : YP[c][r];	        
			        YP[c][r] = (YP[c][r]>maxVR) ? maxVR : YP[c][r];                      
		        }

              Line[c] = YP[c][r]>>(16-BD);
              //printf(" Line[%d] = %d,  YP[%d][%d] = %d  ",c,Line[c],c,r,YP[c][r]);
          }            
          //    write line arrary yuvOut.write(
          yuvOut.write((char *)Line,2*arraySizeX);
        }
        // Dz of 420
        for(unsigned short r = 0; r < arraySizeY/2;r++)
        {
          for(unsigned short c = 0; c < arraySizeX/2;c++)
          {
              // Fill line array 
              // (10 bits) Line[c] = (((DzP[c][r]+2)>>2));

		        if(!FULLRANGE) {
			        DzP[c][r] = (DzP[c][r]<minVRC) ? minVRC : DzP[c][r];	        
			        DzP[c][r] = (DzP[c][r]>maxVRC) ? maxVRC : DzP[c][r];                      
		        }
              Line[c] = DzP[c][r]>>(16-BD);
          }            
          //    write line arrary yuvOut.write(
          yuvOut.write((char *)Line,2*arraySizeXH);
        }        
        
        // Dx of 420
        for(unsigned short r = 0; r < arraySizeY/2;r++)
        {
          for(unsigned short c = 0; c < arraySizeX/2;c++)
          {
              // Fill line array 
              // (10 bits) Line[c] = (((DxP[c][r]+2)>>2));

		        if(!FULLRANGE) {
			        DxP[c][r] = (DxP[c][r]<minVRC) ? minVRC : DxP[c][r];	        
			        DxP[c][r] = (DxP[c][r]>maxVRC) ? maxVRC : DxP[c][r];                      
		        }
              Line[c] = DxP[c][r]>>(16-BD);
          }            
          //    write line arrary yuvOut.write(
          yuvOut.write((char *)Line,2*arraySizeXH);
        }        
		  yuvOut.close(); // close file	


 
    }
    else {
    	cout << "\nerror\n";
    }

    cout << "Completed" << argv[1] << "\n\n";
}




void Subsample444to420(unsigned short ** src, unsigned short ** dst, short width, short height,short algorithmn,unsigned long minCV, unsigned long maxCV)
{
	


   // algorithmn 0 is box
	if (algorithmn == 0)
	{
		unsigned short C444[4][4];
		unsigned short C420[2][2];

	    for (int strip = 0; strip < height; strip+=4) 
	    {

            for (int pixel = 0; pixel < width; pixel+=4) {
            
            //C444 is row by column
            // src and dst are x by y e.g. column by row (enjoy the madness!)
            C444[0][0]  = src[pixel][strip];
            C444[0][1]  = src[pixel+1][strip];
            C444[0][2]  = src[pixel+2][strip];
            C444[0][3]  = src[pixel+3][strip];                                    

            C444[1][0]  = src[pixel][strip+1];
            C444[1][1]  = src[pixel+1][strip+1];
            C444[1][2]  = src[pixel+2][strip+1];
            C444[1][3]  = src[pixel+3][strip+1]; 

            C444[2][0]  = src[pixel][strip+2];
            C444[2][1]  = src[pixel+1][strip+2];
            C444[2][2]  = src[pixel+2][strip+2];
            C444[2][3]  = src[pixel+3][strip+2]; 

            C444[3][0]  = src[pixel][strip+3];
            C444[3][1]  = src[pixel+1][strip+3];
            C444[3][2]  = src[pixel+2][strip+3];
            C444[3][3]  = src[pixel+3][strip+3]; 
 

            // step line in blocks of 4 pixels:
            // ib is the pixel counter from 1-4 such 
            // that 8*ib is byte location of pixel start
            for (int ib = 0;ib<4;ib+=1) {   


				// Perform box averaging
			   // Create 420 blocks of 2x2 from 4x4
			   // rows = scanlines columns are pixels across
			   // should i do higher precision math here?
			   //  0 1 2 3
			   //  1 
			   //  2
			   //  3
			   C420[0][0] = ((unsigned long)C444[0][0] + (unsigned long)C444[0][1] + (unsigned long)C444[1][0] + (int)C444[1][1])/4;
			   C420[0][1] = ((unsigned long)C444[0][2] + (unsigned long)C444[0][3] + (unsigned long)C444[1][2] + (int)C444[1][3])/4;
			   C420[1][0] = ((unsigned long)C444[2][0] + (unsigned long)C444[2][1] + (unsigned long)C444[3][0] + (int)C444[3][1])/4;
			   C420[1][1] = ((unsigned long)C444[2][2] + (unsigned long)C444[2][3] + (unsigned long)C444[3][2] + (int)C444[3][3])/4;

           // write 2x2 Dz elements
           //unsigned char Dz22[4];
           dst[pixel/2 ][strip/2] = C420[0][0];
           dst[1+pixel/2 ][strip/2] = C420[0][1];
           dst[pixel/2 ][1+strip/2] = C420[1][0];
           dst[1+pixel/2 ][1+strip/2] = C420[1][1];
			  }}}
	
    } else {  // FIR filter subsampling
    
    // perform 444 to 422 conversion horizontally
    // array is aligned as src & dst[width][height]
    printf("FIR Filter subsampling...\n");
    
	int w, i, j, jm6, jm5, jm4, jm3, jm2, jm1;
	int jp1, jp2, jp3, jp4, jp5, jp6;
	
	// Allocate dst422 storage
	short w422 = width>>1;
   unsigned short** dst422 = (unsigned short**) malloc(w422*sizeof(unsigned short*)); // cols
        for (i = 0; i < w422; i++)
            dst422[i] = (unsigned short*) malloc(height*sizeof(unsigned short)); //rows
	
  int im5, im4, im3, im2, im1, ip1, ip2, ip3, ip4, ip5, ip6;
   float scale = 512.0;
	float c21 = 21.0/scale;
	float c52  = 52.0/scale;
	float c159 =  159.0/scale;
	float c256 =  256.0/scale;

	float temp;

	 for (j=0; j<height; j++)
	 {
	   for (i=0; i<width; i+=2)
	   {
        im5 = (i<5) ? 0 : i-5;
        im3 = (i<3) ? 0 : i-3;
        im1 = (i<1) ? 0 : i-1;
        ip1 = (i<width-1) ? i+1 : width-1;
        ip3 = (i<width-3) ? i+3 : width-1;
        ip5 = (i<width-5) ? i+5 : width-1;
	

	     temp = c21*((float)(src[im5][j])+(float)(src[ip5][j]))
	     					-c52*((float)(src[im3][j])+(float)(src[ip3][j]))
                     +c159*((float)(src[im1][j])+(float)(src[ip1][j]))
                     +c256*((float)(src[i][j]))+0.5;
                         
	     if(temp>maxCV)temp = maxCV;
	     if(temp<minCV)temp = minCV;
	     dst422[(i>>1)][j] = (unsigned short)temp;
	     
	     //printf("s:%d  d:%d  ",src[i][j],dst422[(i>>1)][j]);
	   }
    }
 
	
	// perform 422 to 420 conversion vertically
      float c228 = 228.0/scale;
      float c70  =  70.0/scale;
      float c37  =  37.0/scale;
            c21  =  21.0/scale;
      float c11  =  11.0/scale;
      float c5   =   5.0/scale;

	
		for (i=0; i<w422; i++)
		{
			for (j=0; j<height; j+=2)
			{
				jm5 = (j<5) ? 0 : j-5;
				jm4 = (j<4) ? 0 : j-4;
				jm3 = (j<3) ? 0 : j-3;
				jm2 = (j<2) ? 0 : j-2;
				jm1 = (j<1) ? 0 : j-1;
				jp1 = (j<height-1) ? j+1 : height-1;
				jp2 = (j<height-2) ? j+2 : height-1;
				jp3 = (j<height-3) ? j+3 : height-1;
				jp4 = (j<height-4) ? j+4 : height-1;
				jp5 = (j<height-5) ? j+5 : height-1;
				jp6 = (j<height-6) ? j+6 : height-1;
				
				/* FIR filter with 0.5 sample interval phase shift */
				temp = c228*((float)(dst422[i][j])+(float)(dst422[i][jp1]))
										  +c70*((float)(dst422[i][jm1])+(float)(dst422[i][jp2]))
										  -c37*((float)(dst422[i][jm2])+(float)(dst422[i][jp3]))
										  -c21*((float)(dst422[i][jm3])+(float)(dst422[i][jp4]))
										  +c11*((float)(dst422[i][jm4])+(float)(dst422[i][jp5]))
										  + c5*((float)(dst422[i][jm5])+(float)(dst422[i][jp6]))+0.5;
	     if(temp>maxCV)temp = maxCV;
	     if(temp<minCV)temp = minCV;
	     dst[i][(j>>1)] = (unsigned short)temp;

			}

		}
    
    }
	
}






