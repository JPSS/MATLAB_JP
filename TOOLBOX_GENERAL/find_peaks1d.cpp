
 /*[filtered image] = my_highlowfilter[raw image, w, sigma, n ] n is number of standart defiation for gaussian filter  */
#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "stdlib.h"
#include "DataProcessing.h"
#include <vector>


void mexFunction (
                  int nlhs,       mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]
                  )
{
	//variables  
  	int width, absolute;
	double h_min;
  	double *xin,*peaks_out;
  	std::vector<double> x, peaks;
	
  	/* Check for proper number of arguments */
  	if (nrhs != 4 ) {
   		 mexErrMsgTxt ("find_peaks1d requires: x, width, h_min, use absolute YES=1/NO=0");
  	} else if (nlhs != 1 ) {
    		mexErrMsgTxt ("find_peaks1d returns peaks");
  	}
 
	//read input
	xin = mxGetPr(prhs[0]);
	width = (int) mxGetScalar(prhs[1]);
	h_min = mxGetScalar(prhs[2]);
	absolute = (int)(mxGetScalar(prhs[3]));


	
	//fill vector x
	for(int i=0; i< DataProcessing::max(mxGetM(prhs[0]), mxGetN(prhs[0]) ); i++){
		x.push_back(xin[i]);
	}

	DataProcessing::find_peaks1d(&x,  width, h_min, absolute, &peaks); 

    if (absolute==1) {
        mexPrintf("Length = %i \twidth = %i \th_min = %f\t absolute height \tfound peaks: %i\n", x.size(), width, h_min, peaks.size());
    }
    else {
        mexPrintf("Length = %i \twidth = %i \th_min = %f\t rel. height \tfound peaks: %i\n", x.size(), width, h_min, peaks.size());
    }
	//make the output	 
	plhs[0] = mxCreateDoubleMatrix(peaks.size(),1, mxREAL); 
  	peaks_out = mxGetPr(plhs[0]);

	//fill output
	for(int i=0; i<peaks.size(); i++){
		peaks_out[i] = peaks[i]+1;
	}
}
