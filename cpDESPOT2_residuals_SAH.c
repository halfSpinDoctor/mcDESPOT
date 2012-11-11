/***************************************************************************
 * MEXFUNCTION [res] = despot2_residual_SAH(fv, phaseCycle, data_ssfp, alpha, tr, tefix)
 *
 * CPU-BASED DESPOT2-FM Objective Function
 *
 * Inputs:
 *   fv = [PD T1 T2 Omega]      --> [Np Matrix x 4] Parameter Vectors, where N=nParam (number of parameter trials)
 *   phaseCycle                 --> [scalar] Phase cycle increment of SSFP 
 *   data_ssfp                  --> [Np x 1] MRI Data, # of flip angles by 1
 *   alpha_spgr, alpha_ssfp     --> [scalar] Flip angles (corrected /w fam) in degrees
 *   tr_ssfp                    --> [scalar] TR times, in ms
 *   tefix                      --> [scalar] 0 = do not use sqrt(exp(-TR/T2)) correction factor   1 = use correction factor
 *
 * Outputs:
 *   res
 *
 * Based on Sean Deoni's DESPOT2-FM Paper (but NOT source code).
 *
 * MATLAB COMPILE COMMAND (R2009b, GLNX): mex CFLAGS="\$CFLAGS -std=c99" -lm cpDESPOT2_residuals_SAH.c
 *
 * Samuel A. Hurley
 * University of Wisconsin
 * v2.0 10-Nov-2012
 *
 * Changelog:
 *    v1.0 - initial code, based on cpMCDESPOT_residuals_SAH.c
 *    v1.1 - fixed missing sqrt() in the tefix portion of the code (Jun-2011) 
 *    v2.0 - Refactor to allow arbitrary SSFP phase cycle to be selected
 ***************************************************************************/

/* Includes MEX for Matlab and Lib MATH headers */
#include <mex.h>
#include <math.h>

// Program Header
#include "cpDESPOT2_residuals_SAH.h"

// Debug Flag

// Define maximum number of SPGR and SSFP data points we can have,
// in order to use constant memory effectively
#define MAX_ALPHA_SSFP 40

// Universal Constants
#define PI         3.14159265358979323846264338327950288419716939937510
#define DEG_TO_RAD 0.017453292519943

// Basic Math Functions
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))

// ==== Global variables (for data shared between all fv evaluations) ====
double d_phaseCycle[1];
double d_data_ssfp[MAX_ALPHA_SSFP];
double d_alpha_ssfp[MAX_ALPHA_SSFP];
double d_tr_ssfp[1];
double d_tefix[1];

int    d_nAlphaSSFP[1];

/* Main Function  ************************************************************/
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  
  /**** 0. Variable Declrations **********************************/
  int i;
  int nParam, nAlphaSSFP;
  
  // Size of data vectors in memory
  size_t paramSize;
  size_t ssfpSize;
  
  // --- Host Data --- //
  // Parameter vector
  double **fv;
  
  // MR Data & Control Paramters (Defined as global vars above)
  double *phaseCycle;
  double *data_ssfp;
  double *alpha_ssfp;
  double *tr_ssfp;
  double *tefix;
  
  // Residual solutions
  double *res;
  
  /* I. Error checking ****************************************************************/
  
  // Check for 6 inputs
  if (nrhs != 6)
    mexErrMsgTxt("Requires 6 inputs: [fv phaseCycle data_ssfp alpha_ssfp tr_ssfp tefix]");
  
  // Check for 1 outputs
  if (nlhs != 1)
    mexErrMsgTxt("Requires one output: [res]");
  
  // Check that input #1 is a [4 x N] matrix
  if (mxGetN(prhs[0]) != 4) {
    mexErrMsgTxt("First input (paramter vector) must be 4xN (4 columns wide)");
  }
  
  // Check that all other inputs are Nx1 (single column)
  for (i=1; i<6; i++) {
   if (mxGetN(prhs[i]) > 1)
     mexErrMsgTxt("All inputs must be Nx1 (one column wide)");
  }
  
  // Grab nParam -- number of different paramter guesses
  nParam = mxGetM(prhs[0]);
  
  // Grab the number of SSFP flip angles
  nAlphaSSFP = mxGetM(prhs[3]);            // prhs[3] is alpha_ssfp
  
  // Check that the number of supplied flip angles matches
  if (mxGetM(prhs[2]) != nAlphaSSFP)
    
    mexErrMsgTxt("Number of supplied SSFP flip angles does not match number of SSFP data points");
  
  // Check that the TR and phasCycle are sca
  if (mxGetN(prhs[1]) !=1 || mxGetM(prhs[1]) !=1)
    mexErrMsgTxt("phaseCycle must be a scalar (in degrees)");
  
  if (mxGetN(prhs[4]) !=1 || mxGetM(prhs[4]) !=1)
    mexErrMsgTxt("TR must be a scalar (in sec)");
  
  // Check that tefix is a scalar
  if (mxGetN(prhs[5]) !=1 || mxGetM(prhs[5]) !=1)
    mexErrMsgTxt("tefix must be a scalar (0 or 1)");
  
  /* II. Load In Data *******************************************************/
  
  // Size of Param Vector
  paramSize = sizeof(double)*nParam;
  
  // Alloc Param Vector
  fv      = (double**) mxMalloc(paramSize);
  for (i=0; i<nParam; i++) {
    fv[i] = (double*)  mxMalloc(4*sizeof(double));
  }
  
  // Alloc Output Vector as mxDoubleMatrix
  plhs[0]          = mxCreateDoubleMatrix(0, 0, mxREAL);
  mxSetM(plhs[0],    nParam);
  mxSetN(plhs[0],    1);
  mxSetData(plhs[0], mxMalloc(paramSize));
  
  // Grab pointers to output data
  res              = mxGetPr(plhs[0]);

  // Size of SPGR & SSFP Data
  ssfpSize = sizeof(double)*nAlphaSSFP;
  
  // Load in parameter data matrix
  load_mrhs(fv, 0, prhs);
  
  // Load in MRI data
  phaseCycle    = mxGetPr(prhs[1]);
  data_ssfp     = mxGetPr(prhs[2]);
  alpha_ssfp    = mxGetPr(prhs[3]);
  tr_ssfp       = mxGetPr(prhs[4]);
  tefix         = mxGetPr(prhs[5]);
  
  if (!(tefix[0] == 0.0 || tefix[0] == 1.0)) {
    mexErrMsgTxt("Error: tefix must be 0 or 1");
  }
  
  // Check that the size of spgr/ssfp data are not larger than the hard-coded maxima 
  if (ssfpSize > sizeof(double)*MAX_ALPHA_SSFP) {
    mexErrMsgTxt("Error: This algorithm has a hard-coded limit of 40 SSFP data points, and you used too many!\nDo not panic, the authorities are on their way...");
  }
  
  // Copy SSFP data into constant memory 
  d_nAlphaSSFP[0] = nAlphaSSFP;
  d_phaseCycle    = (double) phaseCycle[0];
  d_tr_ssfp[0]    = (double) tr_ssfp[0];
  d_tefix[0]      = (double) tefix[0];
  
  for (i=0; i<nAlphaSSFP; i++) {
    d_data_ssfp[i]  = data_ssfp[i];
    d_alpha_ssfp[i] = alpha_ssfp[i];
  }
  
  /* III. Compute Residuals *****************************************************/
  
  for (i=0; i<nParam; i++) {
    // SSFP
    calcDESPOT2(fv[i], &res[i]);
  }
  
  /* IV. Cleanup & Free Memory ************************************************/
  
  // Parameter vectors
  for (i=0; i<nParam; i++) {
    mxFree(fv[i]);
  }
  mxFree(fv);
  
} // </Main mexFunction>


/* CPU Functions  ***********************************************************/

void calcDESPOT2(double* d_fv, double* d_res) {
  
  #ifdef DEBUG_FLAG
    mexPrintf("== Entering calcDESPOT2 ==\n");
  #endif
  
  // Declare Some Variables
  int i, j;
  double sina[MAX_ALPHA_SSFP];
  double cosa[MAX_ALPHA_SSFP];
  double ssfp_signal[MAX_ALPHA_SSFP];
  double Mx, My;
  double sinb, cosb, beta, denom;
  
  // Pull out model parameters
  double PD    = d_fv[0];
  double T1    = d_fv[1];
  double T2    = d_fv[2];
  double Omega = d_fv[3];
  
  // Flip angle terms
  double E1    = exp(-d_tr_ssfp[0] / T1);
  double E2    = exp(-d_tr_ssfp[0] / T2);
  
  // Sines and Cosines
  for (i=0; i<d_nAlphaSSFP[0]; i++) {
    sina[i]    = sin(d_alpha_ssfp[i]);
    cosa[i]    = cos(d_alpha_ssfp[i]);
  }
  
  // SSFP Signals
  beta  = Omega*2*PI*d_tr_ssfp[0] + (d_phaseCycle[0]*DEG_TO_RAD);
  
  sinb = sin(beta);
  cosb = cos(beta);
  
  for (i=0; i<d_nAlphaSSFP[0]; i++) {
    denom  = (1-E1*cosa[i]) * (1-E2*cosb) - E2*(E1-cosa[i])*(E2-cosb);
    
    Mx = PD * ((1-E1) * E2 * sina[i] * sinb)      / denom;
    My = PD * ((1-E1) * E2 * sina[i] * (cosb-E2)) / denom;
     
    if (d_tefix[0] == 1) {
      // With extra sqrt(E2) term to account for center-echo readout
      ssfp_signal[i] = sqrt(Mx*Mx + My*My) * sqrt(exp(d_tr_ssfp[0]/T2));
       
    } else {
      // Standard Freeman-Hill Formula
      ssfp_signal[i] = sqrt(Mx*Mx + My*My);
    }
  }
  
  // Compute Residual
  d_res[0] = 0.0;
  for (i=0; i<d_nAlphaSSFP[0]; i++) {
     d_res[0] += pow(ssfp_signal[i]   - d_data_ssfp[i],   2.0);
  }
  
  return;
}


/* Host Helper Functions  ****************************************************/

/* Load in data from prhs[i] matrix and convert to an matrix of doubles *******/
/* P. Mossahebi v1.0 12-Jul-2010 */
void load_mrhs(double** output_double, int argin, const mxArray *prhs[])
{
  // Get a pointer to the double-floating point data
  double *input_double = mxGetPr(prhs[argin]);
  
  // Get the number of elements
  int nSize = mxGetN(prhs[argin]);
  int mSize = mxGetM(prhs[argin]);
  
  int i, j;
  for (i = 0; i < nSize; i++) {
    for(j = 0; j < mSize; j++) {
      output_double[j][i] = input_double[(j*nSize)+i];
    }
  }
}
