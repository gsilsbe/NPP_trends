#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <dirent.h> 

float main(float adg_443, 
	   float aph_443, 
	   float bbp_443, 
	   float bbp_s, 
    	   float chl)
{

int w; /* Wavelength step */ 

/* Assumed solar PAR spectrum */
float PAR_spectrum[] = { 0.00227, 0.00218, 0.00239, 0.00189, 0.00297, 0.00348, 
                         0.00345, 0.00344, 0.00373, 0.00377, 0.00362, 0.00364, 
                         0.00360, 0.00367, 0.00354, 0.00368, 0.00354, 0.00357, 
                         0.00363, 0.00332, 0.00358, 0.00357, 0.00359, 0.00340, 
                         0.00350, 0.00332, 0.00342, 0.00347, 0.00342, 0.00290, 
                         0.00314 };
    
/* Wavelength [nm] */
float wv[] = { 400.0, 410.0, 420.0, 430.0, 440.0, 450.0, 460.0, 470.0, 480.0, 490.0, 500.0, 510.0, 520.0, 
             530.0, 540.0, 550.0, 560.0, 570.0, 580.0, 590.0, 600.0, 610.0, 620.0, 630.0, 640.0, 650.0, 
             660.0, 670.0, 680.0, 690.0, 700.0 };

/* Spectral shape of aphi (Bricaud et al. 1998) */
float A_Bricaud[] = { 0.0241, 0.0287, 0.0328, 0.0359, 0.0378, 0.0350, 0.0328, 
                      0.0309, 0.0281, 0.0254, 0.0210, 0.0162, 0.0126, 0.0103, 
                      0.0085, 0.0070, 0.0057, 0.0050, 0.0051, 0.0054, 0.0052, 
                      0.0055, 0.0061, 0.0066, 0.0071, 0.0078, 0.0108, 0.0174, 
                      0.0161, 0.0069, 0.0025 };

/* Spectral shape of aphi (Bricaud et al. 1998) */
float E_Bricaud[] = { 0.6877, 0.6834, 0.6664, 0.6478, 0.6266, 0.5993, 0.5961, 
                      0.5970, 0.5890, 0.6074, 0.6529, 0.7212, 0.7939, 0.8500, 
                      0.9036, 0.9312, 0.9345, 0.9298, 0.8933, 0.8589, 0.8410, 
                      0.8548, 0.8704, 0.8638, 0.8524, 0.8155, 0.8233, 0.8138, 
                      0.8284, 0.9255, 1.0286 };

/* ----------------------------------------------
	Step 1: Declare all variables and dependent functions                       */

/* Inherent Optical Properties*/

float aphi[31];        /* Phytoplankton absorption coefficient [m-1] */
float adg[31];         /* CDOM and detrital absorption coefficient [m-1]   */ 
float bbp[31];         /* particulate backscaterring coefficient [m-1] */
float aph_mean = 0.0;  
float adg_mean = 0.0;  
float bbp_mean = 0.0;       


/*------------------------------------------------------------------------------  
	Step 2: Derive IOPs at 10 nm increments from 400 to 700 nm                 
 
	Comments:
	IOP data from GIOP-DC (Werdell et al. 2013)
	GIOP-DC assumes slope of adg = 0.018 [m-1]
	GIOP-DC assumes spectral shape of phyto absorption coefficient is a function 
	of Chl (Bricuad et al. 1998)                                                */ 
for (w=0; w<31; w++){
  aphi[w] = aph_443 * A_Bricaud[w] * pow(chl, E_Bricaud[w]) / (0.03711 * pow(chl, 0.61479));
  adg[w] = adg_443 * exp(-0.018 * (wv[w] - 443.0));
  bbp[w]  = bbp_443 * powf(443.0 / wv[w], bbp_s);
}

/* -----------------------------------------------------------------------------  
  Step 3: Integrate over par spectrum                                         */ 

for (w=0; w<30; w++){
  aph_mean += 5.0 * (PAR_spectrum[w] * aphi[w] + PAR_spectrum[w+1] * aphi[w+1]);
  adg_mean += 5.0 * (PAR_spectrum[w] * adg[w] + PAR_spectrum[w+1] * adg[w+1]);
  bbp_mean += 5.0 * (PAR_spectrum[w] * bbp[w] + PAR_spectrum[w+1] * bbp[w+1]);
}

float *results = (float *)malloc(sizeof(float)*3);

results[0] = aph_mean;
results[1] = adg_mean;
results[2] = bbp_mean;

return *results;

free(results);
}
