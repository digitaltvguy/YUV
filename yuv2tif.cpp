// set up for 12 bits 
// ASSUMES 3840x2160 ONLY


#include <tiffio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <math.h>

using namespace std;

double PQ10000_f( double V);
double PQ10000_r( double L);
void Subsample420to444(unsigned short ** src, unsigned short ** dst, short width, short height,short algorithmn,unsigned short minCV, unsigned short maxCV);

void D2020CL(int Luma, unsigned short CrB, unsigned short CrR, int &Red, int &Green, int &Blue, short FULLRANGE, short Gamma, unsigned short Half, unsigned short Full);


// Globals
  uint32 YnegRMx=0, YnegBMx=0;
  int Yav;
  int YavSave;
  int Ybar;
  int Y1,Y2,Y3,Y4;
  int Dz;
  int Dx;
  int Rp,Bp;
  uint32 invalidPixels = 0;
  
  // file pointer for tifstripsize = 3840*2*4;
  TIFF* tif;
  
  // hard code 3840 x 2160
  unsigned short stripsize;
  unsigned short pixelStart = 0;
  unsigned short numStrips = 2160;
  short  stripStart;
  ifstream yuvIn;

main(int argc, char* argv[])
{ 
  
  // set up to handle inverse color difference
  short D709 = 0;
  short D2020 = 0;
  short D2020C = 0;
  float tmpF = 0.0;
  short HD = 0; // ==1 if 1920x1080 cutout
  short qHD = 0; // ==1 if 960x540 cutout
  int frames = -999;
  short IPixF = 0;
  short Y500 =0;
  short Y100 =0;
  short DXYZ = 1;
  short BD = 12; // 10 bit mode
  short Gamma = 10000;
  short SR=4;
  unsigned short Half = 2048.0; // e.g half value 12 bits
  unsigned short Full = 4096.0; //  
  short XX = 0;
  unsigned short minCV = 0; //12 bits
  unsigned short maxCV = 4095;
  short FIR = 1;
  short FULLRANGE = 0; // default is video range
  short ALPHA = 0; // no alpha channel is default
  short numChan2 = 3*2; // only 3 channels rgb with no alpha 
  short numChan = 3; 
  unsigned short D = 4; //D=4 for 12 bits =1 10bits = 16 for 14 bits
  unsigned short minVR, maxVR, minVRC,maxVRC;
  // video range is
  // Luma and R,G,B:  CV = Floor(876*D*N+64*D+0.5)
  // Chroma:  CV = Floor(896*D*N+64*D+0.5)	     
	         
  //Process Args:
  short arg = 2;
  while(arg < argc) {


     
     if(strcmp(argv[arg],"-h")==0) {
     	 printf("\n ARGS:\n 709 (use Rec709 Color Dif)\n 2020 (use Rec2020 Color Dif)\n HD1920 (Format is 1920x1080 images)\n HD960 (Format is 960x540 images)\n (no args get Y'DzDx color difference and 3840x2160 cutout)\n\n\n");
     	 exit(0);
     }
  	
	  if(strcmp(argv[arg],"BOX")==0)FIR = 0;
     if(strcmp(argv[arg],"FULL")==0)FULLRANGE = 1;
     if(strcmp(argv[arg],"ALPHA")==0){
     		ALPHA = 1;
     		numChan2= 4*2;
     		numChan = 4;
     }	
		     		     
     if(strcmp(argv[arg],"709")==0)D709 = 1;
     if(strcmp(argv[arg],"2020")==0)D2020 = 1;
     if(strcmp(argv[arg],"2020C")==0)D2020C = 1;
     
     if(strcmp(argv[arg],"HD1920")==0)HD = 1;
     if(strcmp(argv[arg],"HD960")==0)qHD = 1;  
     
		     
     if(strcmp(argv[arg],"Y100")==0)Y100 = 1;
     if(strcmp(argv[arg],"Y500")==0)Y500 = 1;	
     
     if(strcmp(argv[arg],"B10")==0) {
     	  BD = 10;
     	  SR=6;
     	  Half = 512.0;
     	  Full = 1024.0;
     	  minCV = 0;
     	  maxCV = 1023;
     	  D=1;
  	     printf("\nprocessing line data 10 bits\n");		     	  
     }	     

     if(strcmp(argv[arg],"B14")==0) {
		     	  BD = 14;
		     	  SR=2;
		     	  Half = 8192;
		     	  Full = 16384;
		     	  minCV = 0;
		     	  maxCV = 16383;
		     	  D=16;
  	     printf("\nprocessing line data 14 bits\n");		     	  
     }
     
 		     if(strcmp(argv[arg],"G1k")==0){
				 Gamma = 1;
				 printf("\nUsing Gamma PQ 1k\n");
			 }
		     if(strcmp(argv[arg],"G10k")==0){
				 Gamma = 10000;
				 printf("\nUsing Gamma PQ 10k\n");
			 }
		     if(strcmp(argv[arg],"G1886")==0){
				 Gamma = 1886;
				 printf("\nUsing Gamma Red1886\n");
			 }
 
     if(D709 || D2020 || D2020C || Y100 || Y500)DXYZ = 0;	      
     
     if(strcmp(argv[arg],"-I")==0){
		 IPixF = 1;
		 printf("Invalid Pixel printing enabled\n");
	 }
     if(strcmp(argv[arg],"-X")==0)XX = 1;
     
     if(strcmp(argv[arg],"-f")==0) {
			arg++;
			if(arg < argc)frames=atoi(argv[arg]);		     	
     }
     
   arg++;
  }
	 if(DXYZ)printf("Processing for YDzDx with no comp\n");
	 if(D709)printf("Processing for Rec709\n");
	 if(D2020)printf("Processing for Rec2020\n");
	 if(D2020C)printf("Processing for Rec2020 Constant Luminance\n");
	 if(DXYZ)printf("Processing for YDzDx\n");
	 if(Y100)printf("Processing for Y100 comp of YDzDx\n");
	 if(Y500)printf("Processing for Y500 comp of YDzDx\n");
     if(HD)printf("Processing for HD1920x1080\n");
     if(qHD)printf("Processing for HD960x540\n");	 
  
  // Set up for video range if needed
  if(!FULLRANGE) {
  	  printf("Processing for Video Range\n");
     minVR = 64*D;
     maxVR = 876*D+minVR;
     minVRC = minVR;
     maxVRC = 896*D+minVRC;                
     //achromatic point for chroma will be "Half"(e.g. 512 for 10 bits, 2048 for 12 bits etc..)
  
  }

  stripsize = 3840*numChan2;
  stripStart =0;
  pixelStart = 0;
  if(HD) {
  	stripsize = 1920*numChan2;
  	numStrips = 1080;
  }
  if(qHD) {
  	stripsize= 960*numChan2;
  	numStrips = 540;
  }
  printf ("Stripsize (bytes): %d, %d (pixels), numStrips %d \n",stripsize, stripsize/numChan2, numStrips);
  
  // Saving Invalid Pixel Data as InvalidPixel.txt:
  ofstream invPix;
  if(IPixF) invPix.open("InvalidPixel.txt");

  // set up alternate differencing equations for Y'DzwDxw
  float T,U,V,W;
  if (Y100) {
  	  printf("Setting TUVW for 100nit PQ color difference point\n");
     T = 0.98989899;
     U = 2.0;
     V = 1.016835017;
     W = 2.03367;  
  } else if(Y500) {
  	  printf("Setting TUVW for 500nit PQ color difference point\n");
     T = 0.99203764;
     U = 2.0;
     V = 1.013391241;
     W = 2.026782;	     	
  }
  
  // Readign array (stripsize/numChan2 = 3840 unsigned shorts)
  unsigned short *yuvLine;
  // Array to store line of output for writing process
  // will be allocated to line width with 4 unsigned shorts
  unsigned short *Line;
  
  // Allocate memory to read each image
  unsigned short** YP;
  int arraySizeX = stripsize/numChan2 - pixelStart/numChan2; // eg 3840
  int arraySizeY = numStrips;
  int arraySizeXH = arraySizeX/2; // eg 1920 cols
  int arraySizeYH = arraySizeY/2; // eg 1080 rows
  printf("Frame Size: %d x %d\n",arraySizeX,arraySizeY);
  YP = (unsigned short**) malloc(arraySizeX*sizeof(unsigned short*)); // cols
  for (int i = 0; i < arraySizeX; i++)
      YP[i] = (unsigned short*) malloc(arraySizeY*sizeof(unsigned short)); //rows

// allocate for 444 chroma

  unsigned short** Cb444;
  Cb444 = (unsigned short**) malloc(arraySizeX*sizeof(unsigned short*)); // cols
  for (int i = 0; i < arraySizeX; i++)
      Cb444[i] = (unsigned short*) malloc(arraySizeY*sizeof(unsigned short)); //rows
 
  unsigned short** Cr444;
  Cr444 = (unsigned short**) malloc(arraySizeX*sizeof(unsigned short*)); // cols
  for (int i = 0; i < arraySizeX; i++)
      Cr444[i] = (unsigned short*) malloc(arraySizeY*sizeof(unsigned short)); //rows

// allocate for 420 chroma
  unsigned short** DzP;
  DzP = (unsigned short**) malloc(arraySizeXH*sizeof(unsigned short*)); // cols
  for (int i = 0; i < arraySizeXH; i++)
      DzP[i] = (unsigned short*) malloc(arraySizeYH*sizeof(unsigned short)); // rows
      
  unsigned short** DxP;
  DxP = (unsigned short**) malloc(arraySizeXH*sizeof(unsigned short*)); // cols
  for (int i = 0; i < arraySizeXH; i++)
      DxP[i] = (unsigned short*) malloc(arraySizeYH*sizeof(unsigned short)); // rows
 
  // set output line array to unsigned short              
  Line =  (unsigned short *)malloc(((numChan2/2)*arraySizeX*sizeof(unsigned short)));

  yuvLine = (unsigned short *) malloc(arraySizeX*sizeof(unsigned short));



	// Open yuv file
	// Open Binary PlanarYUV file for reading:
	yuvIn.open(argv[1], ios::in | ios::binary);
	printf("Opened YDzDx.yuv reading...:\n");
	 
	// process yuv
	int tifNum = 0;  
	char tifName[] = "tifXYZ/XpYpZp000000.tif";  // 6 digits	
	int line = 0;
	int pixel = 0;
	while(yuvIn)
	{
		YnegRMx=0;
		YnegBMx=0;
		// read Y' data
		for ( line = 0;line < arraySizeY;line++)
		{
			yuvIn.read((char *)yuvLine, arraySizeX*sizeof(unsigned short));
			for ( pixel = 0; pixel < arraySizeX;pixel++) {
				YP[pixel][line] = yuvLine[pixel];
				//printf(" YP[%d][%d]= %d ",pixel,line,YP[pixel][line]);
				if(!FULLRANGE) {
					YP[pixel][line] = (YP[pixel][line]<minVR) ? minVR : YP[pixel][line];
					YP[pixel][line] = (YP[pixel][line]>maxVR) ? maxVR : YP[pixel][line];
				}
			}
		}
		
		// read Dz data
		for ( line = 0;line < arraySizeY/2;line++)
		{
			yuvIn.read((char *)yuvLine, arraySizeX*sizeof(unsigned short)/2);
			for ( pixel = 0; pixel < arraySizeX/2;pixel++){
				DzP[pixel][line] = yuvLine[pixel];
				if(!FULLRANGE) {
					DzP[pixel][line] = (DzP[pixel][line]<minVRC) ? minVRC : DzP[pixel][line];
					DzP[pixel][line] = (DzP[pixel][line]>maxVRC) ? maxVRC : DzP[pixel][line];
				}				
			}
		}
				
		// read Dx data
		for ( line = 0;line < arraySizeY/2;line++)
		{
			yuvIn.read((char *)yuvLine, arraySizeX*sizeof(unsigned short)/2);
			for ( pixel = 0; pixel < arraySizeX/2;pixel++) {
				DxP[pixel][line] = yuvLine[pixel];
				if(!FULLRANGE) {
					DxP[pixel][line] = (DxP[pixel][line]<minVRC) ? minVRC : DxP[pixel][line];
					DxP[pixel][line] = (DxP[pixel][line]>maxVRC) ? maxVRC : DxP[pixel][line];
				}	
			}
		}		

		printf("Writing tifXYZ/XpYpZp%06d.tif\n",tifNum);
		// Open TIF File
		sprintf(tifName, "tifXYZ/XpYpZp%06d.tif",tifNum);
		invalidPixels = 0;
		
		tif = TIFFOpen(tifName, "w");
		if(ALPHA) {
		   TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 4);
		} else {
		  TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 3);
		}
		
		TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 16);
		TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
		TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, arraySizeX);
		TIFFSetField(tif, TIFFTAG_IMAGELENGTH, arraySizeY);
		TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, 1);
		TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);	
		
		
		if(!FULLRANGE) { // VideoRange
		Subsample420to444(DzP, Cb444, arraySizeX, arraySizeY, FIR, minVRC,  maxVRC);
		Subsample420to444(DxP, Cr444, arraySizeX, arraySizeY, FIR, minVRC,  maxVRC);
		} else { //Fullrange
		Subsample420to444(DzP, Cb444, arraySizeX, arraySizeY, FIR, minCV,  maxCV);
		Subsample420to444(DxP, Cr444, arraySizeX, arraySizeY, FIR, minCV,  maxCV);
		}  		
		
		for (int line = 0;line < arraySizeY;line++)
		{
			for (unsigned int pixel = 0; pixel < numChan*arraySizeX;pixel+=numChan)
			{

           // Y = G = Y
           // Dz = -0.5*G + 0.5*B  = 0.5*Z - 0.5*Y + 2048
           // Dx = 0.5*R -0.5*G  = 0.5*X - 0.5*Y  +2048
           
           // 4095 = 2047.5 * 2
           // recall 2047.5 came from adding 2047 to Dz, Dx 
           //  and adding 8 during initial round off from 16bits to 12

 

           // Try calculating an average Y to use with subsampled Dz, Dx
           //Yav = (YP[pixel>>2][line>>1 + line>>1] + YP[1+ pixel>>2][line>>1 + line>>1] + YP[1+pixel>>2][line>>1 + line>>1] + YP[1 + pixel>>2][1 + line>>1 + line>>1])/4;
           Yav = (int)(YP[pixel/numChan][line]);
           YavSave = Yav; // keep copy of original Y'
                   
           // Calculate Ybar
           if (XX) {
           	 short p0 = (pixel/numChan) % 2;
           	 short l0 = line % 2;
           	 if(p0 == 0 && l0 == 0) {
           	 	Ybar = (int)(YP[pixel/numChan][line]) + (int)(YP[(pixel+numChan)/numChan][line]) + (int)(YP[pixel/numChan][line+1]) + (int)(YP[(pixel+numChan)/numChan][line+1]);
           	 	Ybar = Ybar/4;
           	 }
           	 if(p0 == 1 && l0 == 0) {
           	 	Ybar = (int)(YP[(pixel-numChan)/numChan][line]) + (int)(YP[pixel/numChan][line]) + (int)(YP[(pixel-numChan)/numChan][line+1]) + (int)(YP[pixel/numChan][line+1]);
           	 	Ybar = Ybar/4;
           	 }           	 
           	 if(p0 == 0 && l0 == 1) {
           	 	Ybar = (int)(YP[pixel/numChan][line-1]) + (int)(YP[(pixel+numChan)/numChan][line-1]) + (int)(YP[pixel/numChan][line]) + (int)(YP[(pixel+numChan)/numChan][line]);
           	 	Ybar = Ybar/4;
           	 }  
           	 if(p0 == 1 && l0 == 1) {
           	 	Ybar = (int)(YP[(pixel-numChan)/numChan][line-1]) + (int)(YP[pixel/numChan][line-1]) + (int)(YP[(pixel-numChan)/numChan][line]) + (int)(YP[pixel/numChan][line]);
           	 	Ybar = Ybar/4;
           	 } 
           	 if(Ybar<1)Ybar = 1;
           	 if(Ybar>(Full-1))Ybar=Full-1;          	          	 
           }
       
           // July 6th 2014 changed from (Half-0.5) to Half-1.0 and 
           // in YDzDx from (Full-1.0) to (Half) done before the multiply by 2
           // this mirrors what is in tif2yuv where (Half) is added to chroma after calculation
           // originally had seen something indicating chroma was to +/- 2047.5 but was never sure
           // how to do that. Looking at rec2020 where everything is integer
           // with form like 224*chroma+128 
           if(DXYZ) {
           	  if(XX) {
           	  	  //reconstruct with Ybar then scale
           	  	  float RED,BLUE;
		           RED = (2.0*((float)(Cr444[pixel/numChan][line]) - (float)(Half)) + (float)Ybar) * ((float)Yav)/((float)Ybar); 
		           BLUE = (2.0*((float)(Cb444[pixel/numChan][line]) - (float)(Half)) + (float)Ybar) * ((float)Yav)/((float)Ybar);
		           if(RED>(Full-1.5))RED=Full-1.0;
		           if(BLUE>(Full-1.5))BLUE=Full-1.0;
		           Rp = RED;
		           Bp = BLUE;
		           
           	  } else {
		           Rp = ((int)2*((int)(Cr444[pixel/numChan][line]) - (int)(Half)) + Yav); 
		           Bp = ((int)2*((int)(Cb444[pixel/numChan][line]) - (int)(Half)) + Yav);
		        } 
	        } else if (D2020C) {
				D2020CL((int)(YP[pixel/numChan][line]), Cb444[pixel/numChan][line], Cr444[pixel/numChan][line], Rp, Yav, Bp, FULLRANGE, Gamma, Half, Full);
			}  else if(D2020){
           	  tmpF = ((float)(Cb444[pixel/numChan][line])-(Half))*1.8814 + Yav;
           	  if(tmpF > (Full-1.0))tmpF = (Full-1.0);
           	  Bp = tmpF;
           	  tmpF = ((float)(Cr444[pixel/numChan][line])-(Half))*1.4746 + Yav;
           	  if(tmpF > (Full-1.0))tmpF = (Full-1.0);
           	  Rp = tmpF;
           	  tmpF = ((float)Yav - 0.0593*(float)Bp -0.2627*(float)Rp)/0.6780 +0.5 ; // green
           	  if(tmpF > (Full-1.0)) tmpF = (Full-1.0);
           	  Yav = tmpF; //green
           } else if(D709) {
           	  tmpF = ((float)(Cb444[pixel/numChan][line])-(Half))*1.8556 + Yav;
           	  if(tmpF > (Full-1.0))tmpF = (Full-1.0);
           	  Bp = tmpF;           	  
           	  tmpF = ((float)(Cr444[pixel/numChan][line])-(Half))*1.5748 + Yav;
           	  if(tmpF > (Full-1.0))tmpF = (Full-1.0);
           	  Rp = tmpF;           	  
           	  tmpF = ((float)Yav - 0.07222*(float)Bp -0.2126*(float)Rp)/0.7152 +0.5 ; // green
           	  if(tmpF > (Full-1.0)) tmpF = (Full-1.0);
           	  Yav = tmpF;  //green         	                	
           } else if(Y100 || Y500) {
	           	  tmpF = ((float)(Cb444[pixel/numChan][line])-(Half))*W + ((float)Yav)*V;
	           	  if(tmpF > (Full-1.0))tmpF = (Full-1.0);
	           	  Bp = tmpF;
	           	  tmpF = ((float)(Cr444[pixel/numChan][line])-(Half))*U + ((float)Yav)*T;
	           	  if(tmpF > (Full-1.0))tmpF = (Full-1.0);
	           	  Rp = tmpF; 
           } else {
           	   printf("Can't determine color difference to use?\n\n");
           	   exit(0);
           }	
           
           if(Yav < 0)
           {

           		if(IPixF)invPix << "Pixel=" << pixel/8 << ", " << line << "  Y'=" << YavSave << "  Dz=" << Cb444[pixel/numChan][line] << "  Dx=" << Cr444[pixel/numChan][line] << "  G'=" << Yav << " !! Gneg\n";
           		
           		Yav = 0;
           		invalidPixels++;
           }        
        
           if(Rp < 0 )
           {
           	   //if(IPixF && YavSave>0)
           	   if(IPixF && YavSave>0)invPix << "Pixel=" << pixel/8 << ", " << line << "  Y'=" << YavSave << "  Dz=" << Cb444[pixel/numChan][line] << "  Dx=" << Cr444[pixel/numChan][line] << "  R'=" << Rp << " !! Rneg\n";
           	   
           		if(Yav > YnegRMx && Yav > 0)
           		{
           			YnegRMx = Yav;
           			printf("Line: %d Pixel: %d, Y %d, Dz: %d, Dx %d, X' %d, YnegRMx %d\n", line,pixel/8, Yav, Cb444[pixel/numChan][line], Cr444[pixel/numChan][line], Rp,YnegRMx);
           		}
           		Rp = 0; // zero invalid pixels whether Yav == 0 or not
           		if(YavSave != 0)invalidPixels++;
           }
           
           if(Bp < 0 )
           {
           	   //if(IPixF && YavSave>0)
           	   if(IPixF && YavSave>0)invPix << "Pixel=" << pixel/8 << ", " << line << "  Y'=" << YavSave << "  Dz=" << Cb444[pixel/numChan][line] << "  Dx=" << Cr444[pixel/numChan][line] << "  B'=" << Bp << " !! Bneg\n";
           	
           		if(Yav > YnegBMx && Yav > 0)
           		{
           			YnegBMx = Yav;
           			printf("Line: %d Pixel: %d, Y %d, Dz: %d, Dx %d, Z' %d,YnegBMx %d\n", line,pixel/8, Yav, Cb444[pixel/numChan][line], Cr444[pixel/numChan][line], Bp, YnegBMx);
           		}
           		Bp = 0; // zero invalid pixels whether Yav == 0 or not
           		if(YavSave != 0)invalidPixels++;
           }


              // Insure RGB are clipped to video range if required:
              if(!FULLRANGE) {
              Rp = (Rp<minVR) ? minVR : Rp;
              Yav = (Yav<minVR) ? minVR : Yav;
              Bp = (Bp<minVR) ? minVR : Bp;
              
              Rp = (Rp>maxVR) ? maxVR : Rp;
              Yav = (Yav>maxVR) ? maxVR : Yav;
              Bp = (Bp>maxVR) ? maxVR : Bp;
                            
              } else {
              Rp = (Rp<minCV) ? minCV : Rp;
              Yav = (Yav<minCV) ? minCV : Yav;
              Bp = (Bp<minCV) ? minCV : Bp;
              
              Rp = (Rp>maxCV) ? maxCV : Rp;
              Yav = (Yav>maxCV) ? maxCV : Yav;
              Bp = (Bp>maxCV) ? maxCV : Bp;				  
		      }
		      


           // Calculate Rp from Cr444
           
           
           // Calculate Bp from Cb444
           
           
           
           // R = X = 2*Dx + Y
			  Line[pixel] = ((unsigned short)Rp) << SR;     // R = X
				
				// G = Y
				Line[pixel+1] = ((unsigned short)Yav) << SR;  //G = Y or inverse 2020/709 equation
				
				// B = X = 2*Dx + Y
				Line[pixel+2] = ((unsigned short)Bp) << SR;   // B = Z
				
				// A				
				if(ALPHA) Line[pixel+3] = 65535;  // A

			   //printf("Rp=%d   Gp=%d   Bp=%d | ",Line[pixel],Line[pixel+1],Line[pixel+2]);
								
			}
			

  			//printf("Writing strip %d with width %d bytes %d pixels\n",line,4*arraySizeX*2,arraySizeX);
			TIFFWriteRawStrip(tif, (tstrip_t)line, (tdata_t)Line, numChan2*arraySizeX);

		}		
		
		TIFFClose(tif);
		
		tifNum++;
		
		printf("Max YnegR = %d, Max YnegB = %d\n",YnegRMx, YnegBMx);
		printf("Invalid Pixels:  %d\n",invalidPixels);
		if(frames > 0) {
			frames--;
			if(frames == 0)exit(0);
		}

	}


if(IPixF) invPix.close();

}		


void Subsample420to444(unsigned short ** src, unsigned short ** dst, short width, short height,short algorithmn,unsigned short minCV, unsigned short maxCV)
{
	if(algorithmn == 0) {
		short widthH = width/2;
		short heightH = height/2;
		for (int line = 0;line < heightH;line++)
		{
			for (unsigned int pixel = 0; pixel < widthH;pixel++) {
				dst[2*pixel][2*line] = src[pixel][line];
				dst[2*pixel+1][2*line] = src[pixel][line];
				dst[2*pixel][2*line+1] = src[pixel][line];
				dst[2*pixel+1][2*line+1] = src[pixel][line];				
		   }
	   }
  } else {
   // Implement FIR filter for 420 to 422 then 444  

  int wH,i, j, j2;
  int jm6, jm5, jm4, jm3, jm2, jm1, jp1, jp2, jp3, jp4, jp5, jp6, jp7;

  // width and height come in as of src array


	// Allocate dst422 storage
	short w422 = width>>1;
	short h420 = height>>1;
   unsigned short** dst422 = (unsigned short**) malloc(w422*sizeof(unsigned short*)); // cols
        for (i = 0; i < w422; i++)
            dst422[i] = (unsigned short*) malloc(height*sizeof(unsigned short)); //rows
            
   float scale = 256.0;
   float c3    = 3.0/scale;
   float c16   = 16.0/scale;
   float c67  =  67.0/scale;
   float c227 =  227.0/scale;
   float c32  =   32.0/scale;
   float c7   =   7.0/scale;
   float temp;

  /* intra frame */
  for (i=0; i<w422; i++) // w here is half width (e.g. 422 width)
  {
    for (j=0; j<h420; j++)
    {
        j2 = j<<1;
        jm3 = (j<3) ? 0 : j-3;
        jm2 = (j<2) ? 0 : j-2;
        jm1 = (j<1) ? 0 : j-1;
        jp1 = (j<h420-1) ? j+1 : h420-1;
        jp2 = (j<h420-2) ? j+2 : h420-1;
        jp3 = (j<h420-3) ? j+3 : h420-1;

        /* FIR filter coefficients (*256): 5 -21 70 228 -37 11 */
        /* New FIR filter coefficients (*256): 3 -16 67 227 -32 7 */
        temp =     c3*((float)(src[i][jm3]))
                             -c16*((float)(src[i][jm2]))
                             +c67*((float)(src[i][jm1]))
                            +c227*((float)(src[i][j]))
                             -c32*((float)(src[i][jp1]))
                             +c7*((float)(src[i][jp2]))+0.5;
	     if(temp>maxCV)temp = maxCV;
	     if(temp<minCV)temp = minCV;
	     dst422[i][j2] = (unsigned short)temp;

        temp = c3*((float)(src[i][jp3]))
                             -c16*((float)(src[i][jp2]))
                             +c67*((float)(src[i][jp1]))
                            +c227*((float)(src[i][j]))
                             -c32*((float)(src[i][jm1]))
                             +c7*((float)(src[i][jm2]))+0.5;
	     if(temp>maxCV)temp = maxCV;
	     if(temp<minCV)temp = minCV;
	     dst422[i][j2+1] = (unsigned short)temp;

      }
    }
    
	// 422 to 444 now
   int im2, im1, ip1, ip2, ip3,i2;
   scale = 256.0;
	float c21 = 21.0/scale;
	float c52  = 52.0/scale;
	float c159 =  159.0/scale;
	float c256 =  256.0/scale;

	
    for (j=0; j<height; j++)
    {
      for (i=0; i<w422; i++) // half width (e.g. 422 width)
      {
        i2 = i<<1;
        im2 = (i<2) ? 0 : i-2;
        im1 = (i<1) ? 0 : i-1;
        ip1 = (i<w422-1) ? i+1 : w422-1;
        ip2 = (i<w422-2) ? i+2 : w422-1;
        ip3 = (i<w422-3) ? i+3 : w422-1;

        /* FIR filter coefficients (*256): 21 0 -52 0 159 256 159 0 -52 0 21 */
        /* even samples (0 0 256 0 0) */
        dst[i2][j] = dst422[i][j];


        /* odd samples (21 -52 159 159 -52 21) */
        temp = c21*(((float)(dst422[im2][j]))+((float)(dst422[ip3][j])))
                        -c52*(((float)(dst422[im1][j]))+((float)(dst422[ip2][j]))) 
                       +c159*(((float)(dst422[i][j]))+((float)(dst422[ip1][j])))+0.5;
	     if(temp>maxCV)temp = maxCV;
	     if(temp<minCV)temp = minCV;
	     dst[i2+1][j] = (unsigned short)temp;

      }

    }    
    
  
  }

}


void D2020CL(int Luma, unsigned short CrB, unsigned short CrR, int &Red, int &Green, int &Blue, short FULLRANGE, short Gamma, unsigned short Half, unsigned short Full)
{
	
  const float OUT_WP_MAX = 10000.0;
  const float OUT_WP_10k = 10000.0;
  const float OUT_WP_1k  = 1000.0;
  const uint32 CV_BLACK = 4096; //64.0*64.0;
  const uint32 CV_WHITE = 60160;
  	
	// CrB and CrR to be sent in as unsigned short from 444 array
	// Luma to be sent in as codevalue range
	
	// Recover Blue and Red
	float CrBF = (float)CrB - (Half);
	if(CrBF <= 0.0) {
		CrBF = CrBF * 1.9404;
	}
    else {
		CrBF = CrBF * 1.5816;
	}
	CrBF = CrBF + (float)Luma; // B' = CrB' + Y'
	CrBF = (CrBF<0.0) ? 0.0 : CrBF;
	CrBF = (CrBF>(Full-1.0)) ? (Full-1.0) : CrBF;
	
	float CrRF = (float)CrR - (Half);
	if(CrRF <= 0.0) {
		CrRF = CrRF * 1.7184;
	}
    else {
		CrRF = CrRF * 0.9936;
	}
	
	CrRF = CrRF + (float)Luma;
	CrRF = (CrRF<0.0) ? 0.0 : CrRF;
	CrRF = (CrRF>(Full-1.0)) ? (Full-1.0) : CrRF;	
	
	// now CrBF is Blue' and CrRF is Red'
	// set return integer values
	Red = (int)(CrRF+0.5);
	Blue = (int)(CrBF+0.5);
	
	
	// Recover Green:
	// set Luma, Blue and Red to FULLRANGE 0-1
	  float LumaF = (float)Luma/(Full-1.0);
	  CrRF  = CrRF/(Full-1.0);
	  CrBF  = CrBF/(Full-1.0);

	if(!FULLRANGE) {
	  LumaF = (LumaF - CV_BLACK/65535.0)*(65535.0/(CV_WHITE-CV_BLACK));
	  CrRF = (CrRF - CV_BLACK/65535.0)*(65535.0/(CV_WHITE-CV_BLACK));
	  CrBF = (CrBF - CV_BLACK/65535.0)*(65535.0/(CV_WHITE-CV_BLACK));	  	  
	}
	
	// insure nothing outside 0-1
	// CrRF, CrBF are Red, Blue with Gamma
	LumaF = (LumaF<0.0) ? 0.0 : LumaF;
	LumaF = (LumaF>1.0) ? 1.0 : LumaF;	 
	CrRF = (CrRF<0.0) ? 0.0 : CrRF;
	CrRF = (CrRF>1.0) ? 1.0 : CrRF;	
	CrBF = (CrBF<0.0) ? 0.0 : CrBF;
	CrBF = (CrBF>1.0) ? 1.0 : CrBF;	
	
			
	// remove gamma on Luma, remove gamma on Blue and Red
    // PQ 10k or 1k or 1886
    // Remove gamma:
  if (Gamma == 1) { //remove 1k nit PQ gamma
	  // PQ10000_r(0.1) is OETF value of 1k nits which is max codevalue
	  CrRF = PQ10000_f(PQ10000_r(0.1)*CrRF)*OUT_WP_MAX;
	  LumaF = PQ10000_f(PQ10000_r(0.1)*LumaF)*OUT_WP_MAX;
	  CrBF = PQ10000_f(PQ10000_r(0.1)*CrBF)*OUT_WP_MAX;
	  
	  // make sure not over max
	  CrRF = (CrRF>OUT_WP_1k) ? OUT_WP_1k : CrRF;
	  LumaF = (LumaF>OUT_WP_1k) ? OUT_WP_1k : LumaF;
	  CrBF = (CrBF>OUT_WP_1k) ? OUT_WP_1k : CrBF; 
  	  
  } else if (Gamma == 10000) { // remove 10k nit PQ gamma
	  
	  CrRF = PQ10000_f(CrRF)*OUT_WP_MAX;
	  LumaF = PQ10000_f(LumaF)*OUT_WP_MAX;
	  CrBF = PQ10000_f(CrBF)*OUT_WP_MAX;
	  
	  // make sure not over max
	  CrRF = (CrRF>OUT_WP_10k) ? OUT_WP_10k : CrRF;
	  LumaF = (LumaF>OUT_WP_10k) ? OUT_WP_10k : LumaF;
	  CrBF = (CrBF>OUT_WP_10k) ? OUT_WP_10k : CrBF;
	  
  } else {printf("Error unknown gamma\n");Green = 65535; return;}	
	
	
	
	// recover Green
	float GreenF = (LumaF - 0.0593*CrBF -0.2627*CrRF)/0.6780; 
	if(GreenF < 0 && LumaF > 50.0) {
		//printf("LumaF: %f, CrBF: %f, CrRF: %f, GreenF: %f, CrR = %d, Luma = %d\n",LumaF,CrBF, CrRF, GreenF, CrR, Luma);
		GreenF = 0.0;
	}
	
	// apply gamma to Green
	if (Gamma == 1) {
       GreenF = PQ10000_r(GreenF/OUT_WP_MAX)/PQ10000_r(0.1);
    } else if (Gamma == 10000) {
	   GreenF = PQ10000_r(GreenF/OUT_WP_MAX);
    } else {printf("Error unknown gamma\n");Green = 65535; return;}

   // set Green to video range if needed
   if(!FULLRANGE){
     // set Green to Video Range from 0-1 full range
     // float to integer round
	 GreenF = GreenF *((CV_WHITE-CV_BLACK)/65535.0) + CV_BLACK/65535.0;
	 Green  = (int)(GreenF*(Full-1.0)+0.5);    
   }  else {
	 // set Green to Full Range
	 // float to integer round
	 Green = (int)((Full- 1.0) * GreenF +0.5);
   }

	
	
	return;
}

// 10000 nits
//  1/gamma-ish, calculate V from Luma
// decode L = (max(,0)/(c2-c3*V**(1/m)))**(1/n)
double PQ10000_f( double V) // EOTF
{
  double L;
  // Lw, Lb not used since absolute Luma used for PQ
  // formula outputs normalized Luma from 0-1
  L = pow(std::max(pow(V, 1.0/78.84375) - 0.8359375 ,0.0)/(18.8515625 - 18.6875 * pow(V, 1.0/78.84375)),1.0/0.1593017578);

  return L;
}

// encode (V^gamma ish calculate L from V)
//  encode   V = ((c1+c2*Y**n))/(1+c3*Y**n))**m
double PQ10000_r( double L) // OETF
{
  double V;
  // Lw, Lb not used since absolute Luma used for PQ
  // input assumes normalized luma 0-1
  V = pow((0.8359375+ 18.8515625*pow((L),0.1593017578))/(1+18.6875*pow((L),0.1593017578)),78.84375);
  return V;
}
 

