/***************************************************************************
 * MEXFUNCTION [resSPGR resSSFP_0 resSSFP_180 medSSFP_0 medSSFP_180] = mcDESPOT_residual_SD(T1m, T1f, T2m, T2f, Fm, Tau_m, Omega, PD_spgr, PD_ssfp, data_spgr, data_ssfp_0, data_ssfp_180, alpha_spgr, alpha_ssfp, tr_spgr, tr_ssfp)
 *
 * T1m, T1f, T2m, T2f, Fm, Tau_m, Omega, PD_spgr, PD_ssfp --> [Nx1]  Parameter Vectors, where N=nParam (number of parameter trials)
 * data_spgr, data_ssfp_0 data_ssfp_180 --> [Npx1] MRI Data, # of flip angles by 1
 * alpha_spgr, alpha_ssfp        --> [scalar] Flip angles (corrected /w B1err) in degrees
 * tr_spgr, tr_ssfp              --> [scalar] TR times, in ms
 *
 * Based on Sean Deoni's mcDESPOT C-Code
 *
 * COMPILE COMMAND: mex CFLAGS="\$CFLAGS -std=c99" -lm cpMCDESPOT_residuals_SD.c
 *
 * Samuel A. Hurley
 * Pouria Mossahebi
 * University of Wisconsin
 * v3.1 Jul-2010
 *
 * Changelog:
 *    v1.0 - initial code, based on gpMCDESPOT_residual
 *    v1.1 - added some debugging lines (Mar-2010)
 *    v2.0 - convert generalMatrixInverse to jxMatrixInverse (Jordan Exchange) method
 *    v3.0 - fitting explicit pd term for spgr and ssfp, based on v1.0 (Jun-2010)
 *    v3.1 - Implemented Double capabilities (Jul-2010)
 ***************************************************************************/

/* Includes MEX for Matlab and CUDA headers */
#include <mex.h>
#include <math.h>

// Application Header
#include "cpMCDESPOT_residuals_SD.h"

// Define maximum number of SPGR and SSFP data points we can have,
// in order to use constant memory effectively
#define MAX_ALPHA_SPGR 40
#define MAX_ALPHA_SSFP 40

// Everyone's favourite number...
#define PI 3.1415926535897932384626433
// Re-scale signal before taking difference for residuals
#define SIGNALSCALE 1.00     


// Global variables
double d_data_spgr[MAX_ALPHA_SPGR];
double d_data_ssfp_0[MAX_ALPHA_SSFP];
double d_data_ssfp_180[MAX_ALPHA_SSFP];
double d_alpha_spgr[MAX_ALPHA_SPGR];
double d_alpha_ssfp[MAX_ALPHA_SSFP];
double d_tr_spgr[1];
double d_tr_ssfp[1];
int   d_nAlphaSPGR[1];
int   d_nAlphaSSFP[1];

/* Main Function  ************************************************************/
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  
  /**** 0. Variable Declrations **********************************/
  int i;
  int nParam, nAlphaSPGR, nAlphaSSFP;
  
  // Size of data vectors in memory
  size_t paramSize;
  size_t spgrSize;
  size_t ssfpSize;
  
  // --- Host Data --- //
  // Parameter vector
  double *T1m, *T1f, *T2m, *T2f, *Fm, *Tau_m, *Omega, *PD_spgr, *PD_ssfp;
  
  // Residual solutions
  double *resSPGR, *resSSFP_0, *resSSFP_180;
  double *medSSFP_0, *medSSFP_180;
  
  // MR Data
  double *data_spgr,  *data_ssfp_0, *data_ssfp_180;
  double *alpha_spgr, *alpha_ssfp;
  double *tr_spgr,    *tr_ssfp;
  
  /* I. Error checking ****************************************************************/
  // Check for a single output
  if (nlhs != 5)
    mexErrMsgTxt("Requires five outputs [resSPGR resSSFP_0 resSSFP_180 medSSFP_0 medSSFP_180]");
  
  // Check for 16 inputs
  if (nrhs != 16)
    mexErrMsgTxt("Requires 16 inputs: [T1m, T1f, T2m, T2f, Fm, Tau_m, Omega, PD_spgr, PD_ssfp, data_spgr, data_ssfp_0, data_ssfp_180, alpha_spgr, alpha_ssfp, tr_spgr, tr_ssfp]");
  
  // Check that all inputs are Nx1 (single column)
  for (i=0; i<16; i++) {
   if (mxGetN(prhs[i]) > 1)
     mexErrMsgTxt("All inputs must be Nx1 (one column wide)");
  }
  
  // Grab nParam
  nParam = mxGetM(prhs[0]);
  
  // Check that first nine inputs all have same length
  for (i=0; i<9; i++) {
   if (mxGetM(prhs[i]) != nParam)
     mexErrMsgTxt("Parameter inputs (#1-9) must have same length");
  }
  
  // Grab the number of SPGR & SSFP flip angles
  nAlphaSPGR = mxGetM(prhs[10]);
  nAlphaSSFP = mxGetM(prhs[11]);
  
  // Check that the number of supplied flip angles matches
  if (mxGetM(prhs[9]) != nAlphaSPGR)
    mexErrMsgTxt("Number of supplied SPGR flip angles does not match number of SPGR data points");
  
  if (mxGetM(prhs[10]) != nAlphaSSFP)
    mexErrMsgTxt("Number of supplied SSFP   0-phase flip angles does not match number of SSFP data points");
  
  if (mxGetM(prhs[11]) != nAlphaSSFP)
    mexErrMsgTxt("Number of supplied SSFP 180-phase flip angles does not match number of SSFP data points");
  
  // Check that the TR's are scalars
  if (mxGetN(prhs[12]) !=1 || mxGetM(prhs[14]) !=1)
    mexErrMsgTxt("SPGR flip angle and TR must be a scalar (in milliseconds)");
  
  if (mxGetN(prhs[13]) !=1 || mxGetM(prhs[15]) !=1)
    mexErrMsgTxt("SSFP flip angle and TR must be a scalar (in milliseconds)");
  
  /* II. Load In Data *******************************************************/
  
  // Size of param vectors
  paramSize = sizeof(double)*nParam;
  
  // Alloc param vectors on host
  T1m     = (double*) mxMalloc(paramSize);
  T1f     = (double*) mxMalloc(paramSize);
  T2m     = (double*) mxMalloc(paramSize);
  T2f     = (double*) mxMalloc(paramSize);
  Fm      = (double*) mxMalloc(paramSize);
  Tau_m   = (double*) mxMalloc(paramSize);
  Omega   = (double*) mxMalloc(paramSize);
  PD_spgr = (double*) mxMalloc(paramSize);
  PD_ssfp = (double*) mxMalloc(paramSize);
  
  resSPGR     = (double*) mxMalloc(paramSize);
  resSSFP_0   = (double*) mxMalloc(paramSize);
  resSSFP_180 = (double*) mxMalloc(paramSize);
  medSSFP_0   = (double*) mxMalloc(paramSize);
  medSSFP_180 = (double*) mxMalloc(paramSize);

  // Size of SPGR & SSFP Data
  spgrSize = sizeof(double)*nAlphaSPGR;
  ssfpSize = sizeof(double)*nAlphaSSFP;
  
  // Alloc MRI data on host
  data_spgr     = (double*) mxMalloc(spgrSize);
  data_ssfp_0   = (double*) mxMalloc(ssfpSize);
  data_ssfp_180 = (double*) mxMalloc(ssfpSize);
  
  alpha_spgr = (double*) mxMalloc(spgrSize);
  alpha_ssfp = (double*) mxMalloc(ssfpSize);
  
  tr_spgr = (double*) mxMalloc(sizeof(double));  // Scalar value
  tr_ssfp = (double*) mxMalloc(sizeof(double));  // Scalar value
  
  
  // Load in parameter data
  load_rhs(T1m,     0, prhs);
  load_rhs(T1f,     1, prhs);
  load_rhs(T2m,     2, prhs);
  load_rhs(T2f,     3, prhs);
  load_rhs(Fm,      4, prhs);
  load_rhs(Tau_m,   5, prhs);
  load_rhs(Omega,   6, prhs);
  load_rhs(PD_spgr, 7, prhs);
  load_rhs(PD_ssfp, 8, prhs);
  
  // Load in MRI data
  load_rhs(data_spgr,      9, prhs);
  load_rhs(data_ssfp_0,   10, prhs);
  load_rhs(data_ssfp_180, 11, prhs);
  load_rhs(alpha_spgr,    12, prhs);
  load_rhs(alpha_ssfp,    13, prhs);
  load_rhs(tr_spgr,       14, prhs);
  load_rhs(tr_ssfp,       15, prhs);
  
  
  // Check that the size of spgr/ssfp data are not larger than the hard-coded maxima
  if (spgrSize > sizeof(double)*MAX_ALPHA_SPGR) {
    mexErrMsgTxt("Error: This algorithm has a hard-coded limit of SPGR data points, and you used too many!\nDo not panic, the local authorities are on their way...");
  }
  
  if (ssfpSize > sizeof(double)*MAX_ALPHA_SSFP) {
    mexErrMsgTxt("Error: This algorithm has a hard-coded limit of SSFP data points, and you used too many!\nDo not panic, the local authorities are on their way...");
  }
  
  
  // Copy SSFP & SPGR data into constant memory 
  d_nAlphaSPGR[0] = nAlphaSPGR;
  d_nAlphaSSFP[0] = nAlphaSSFP;
  
  d_tr_spgr[0] = (double) tr_spgr[0];
  d_tr_ssfp[0] = (double) tr_ssfp[0];
  
  for (i=0; i<nAlphaSPGR; i++) {
    d_data_spgr[i]  = data_spgr[i];
    d_alpha_spgr[i] = alpha_spgr[i];
  }
  
  for (i=0; i<nAlphaSSFP; i++) {
    d_data_ssfp_0[i]   = data_ssfp_0[i];
    d_data_ssfp_180[i] = data_ssfp_180[i];
    
    d_alpha_ssfp[i]   = alpha_ssfp[i];
  }

  
  /* III. Data Processing *****************************************************/
  
  // SPGR
  for (i=0; i<nParam; i++) {
    calcSPGR(&T1m[i], &T1f[i], &T2m[i], &T2f[i], &Fm[i], &Tau_m[i], &PD_spgr[i], &resSPGR[i]);
  
    // SSFP with Off-Resonance Correction
    calcSSFP(&T1m[i], &T1f[i], &T2m[i], &T2f[i], &Fm[i], &Tau_m[i], &Omega[i], &PD_ssfp[i], &resSSFP_0[i],   &medSSFP_0[i]  , 0);
    calcSSFP(&T1m[i], &T1f[i], &T2m[i], &T2f[i], &Fm[i], &Tau_m[i], &Omega[i], &PD_ssfp[i], &resSSFP_180[i], &medSSFP_180[i], PI);
  }
  
  /* IV. Place results in output vector *************/
  plhs[0] = mxCreateDoubleMatrix(nParam,1,mxREAL);
  plhs[1] = mxCreateDoubleMatrix(nParam,1,mxREAL);
  plhs[2] = mxCreateDoubleMatrix(nParam,1,mxREAL);
  plhs[3] = mxCreateDoubleMatrix(nParam,1,mxREAL);
  plhs[4] = mxCreateDoubleMatrix(nParam,1,mxREAL);
  
  convert_double2double(resSPGR    ,  mxGetPr(plhs[0]),  nParam);
  convert_double2double(resSSFP_0  ,  mxGetPr(plhs[1]),  nParam);
  convert_double2double(resSSFP_180,  mxGetPr(plhs[2]),  nParam);
  convert_double2double(medSSFP_0   , mxGetPr(plhs[3]),  nParam);
  convert_double2double(medSSFP_180 , mxGetPr(plhs[4]),  nParam);
  
  
  /* V. Cleanup & Free Memory ************************************************/
  // Parameter vectors - host
  mxFree(T1m);
  mxFree(T1f);
  mxFree(T2m);
  mxFree(T2f);
  mxFree(Fm);
  mxFree(Tau_m);
  mxFree(Omega);
  mxFree(PD_spgr);
  mxFree(PD_ssfp);
  
  // Output vectors
  mxFree(resSPGR);
  mxFree(resSSFP_0  );
  mxFree(resSSFP_180);
  mxFree(medSSFP_0);
  mxFree(medSSFP_180);
  
  // MRI Data - host
  mxFree(data_spgr);
  mxFree(data_ssfp_0);
  mxFree(data_ssfp_180);
  mxFree(alpha_spgr);
  mxFree(alpha_ssfp);
  mxFree(tr_spgr);
  mxFree(tr_ssfp);
  
} // </Main mexFunction>

/* CPU Functions  ***********************************************************/

// Kernel to initiate residual calcuations
void calcSPGR(double* d_T1m, double* d_T1f, double* d_T2m, double* d_T2f, double* d_Fm, double* d_Tau_m, double* d_PD_spgr, double* d_resSPGR) {
  
  /* Declare Common Variables used within this kernel **********/
  int i, j;
  int a;
  
  // SPGR Signal Calculation
  
  /************<setUpTwoComponentSPGRRelaxationAndRecoveryMatrices>**************/
  
  // Map array paramter values onto Sean's variable naming
	double myelinT1 = d_T1m[0];
	double freeT1   = d_T1f[0];
	// double myelinT2 = d_T2m[idx];  Not required for SPGR
	// double freeT2   = d_T2f[idx];  Not required for SPGR
	double vMyelin  = d_Fm[0];
	double myelinResidenceTime = d_Tau_m[0];
    double spgrProtonDensity = d_PD_spgr[0];
  
	// calculate the dependent parameters
	double vFree = 1.00 - vMyelin;
	double exchangeMyelinToFree = (1.00/ myelinResidenceTime);
	double exchangeFreeToMyelin = ((vMyelin*exchangeMyelinToFree)/ vFree);
  
  //spgr relaxation and recovery matricies
  double A[2][2];
  double invA[2][2] = {0, 0, 0, 0};  // Initilize to zero
  
    A[0][0] = (-1.00/ myelinT1) -  exchangeMyelinToFree;
	A[1][0] = exchangeFreeToMyelin;
	A[0][1] = exchangeMyelinToFree;
	A[1][1] = (-1.00/ freeT1) - exchangeFreeToMyelin;
  
  // Create spgrRelaxation = e^A*TR
  double spgrRelaxation[2][2];
  
  // <---------specialMatrixExponential(spgrRelaxation, A, spgrTR, 2)----------------->

  for (i=0; i<2; i++){
    for (j=0; j<2; j++) spgrRelaxation[i][j] = 0.00;
  }

  double I[2][2];
  double scalarMatrix[2][2];
  double squareMatrix[2][2];
  double tripleMatrix[2][2];
  double quadMatrix[2][2];
  double quintMatrix[2][2];

  // first, apply the scalar to the incoming matrix
  for (i=0; i<2; i++) {
    for (j=0; j<2; j++) scalarMatrix[i][j] = d_tr_spgr[0]*A[i][j];
  }

  // define the identity matrix
  I[0][0] = 1.0;
  I[1][1] = 1.0;
  I[0][1] = 0.0;
  I[1][0] = 0.0;

  // now calculate the square, triple, quad and quint matrices
  // start with the square		
  for (i=0; i<2; i++) {
    for (j=0; j<2; j++) {

      squareMatrix[i][j] = 0.00;
      for (a=0; a<2; a++) squareMatrix[i][j] += ((scalarMatrix[a][j]*scalarMatrix[i][a])/ 2.00);
    }
  }
  // next the triple
  for (i=0; i<2; i++) {
    for (j=0; j<2; j++) {

      tripleMatrix[i][j] = 0.00;
      for (a=0; a<2; a++) tripleMatrix[i][j] += ((squareMatrix[a][j]*scalarMatrix[i][a])/ 3.00);
    }
  }
  // next the quad
  for (i=0; i<2; i++) {
    for (j=0; j<2; j++) {

      quadMatrix[i][j] = 0.00;
      for (a=0; a<2; a++) quadMatrix[i][j] += ((tripleMatrix[a][j]*scalarMatrix[i][a])/ 4.00);
    }
  }
  // and, finally, the quint
  for (i=0; i<2; i++) {
    for (j=0; j<2; j++) {

      quintMatrix[i][j] = 0.00;
      for (a=0; a<2; a++) quintMatrix[i][j] += ((quadMatrix[a][j]*scalarMatrix[i][a])/ 5.00);
    }
  }

  // and, finally, sum them together appropriately
  for (i=0; i<2; i++) {
    for (j=0; j<2; j++) spgrRelaxation[i][j] = I[i][j] + scalarMatrix[i][j] + squareMatrix[i][j] + tripleMatrix[i][j] + quadMatrix[i][j] + quintMatrix[i][j];
  }
  // <---------------------- End specialMatrixExponential ------------------------>
  
  
	//<----------------------generalMatrixInverse(invA, A, 2)----------------------->
	double scalar = (1.00/ (A[0][0]*A[1][1] - A[1][0]*A[0][1]));
	
	invA[0][0] = scalar * A[1][1];
	invA[1][1] = scalar * A[0][0];
	invA[1][0] = scalar * -1.00*A[1][0];
	invA[0][1] = scalar * -1.00*A[0][1];
  //<------------------------End generalMatrixInverse------------------------->
  
  
  // now, tackle the recovery matricies
  // spgrRecovery = inv(A)*C = M0
  double spgrRecovery[2];
	spgrRecovery[0] = (invA[0][0]*spgrProtonDensity*vMyelin/ myelinT1) + (invA[1][0]*spgrProtonDensity*vFree/ freeT1);
	spgrRecovery[1] = (invA[0][1]*spgrProtonDensity*vMyelin/ myelinT1) + (invA[1][1]*spgrProtonDensity*vFree/ freeT1);
  
  /***********</SetupTwoComponantSPGRRelaxationAndRecoveryMatricies>*********************/
  
  
  
  /********************<calculateAllTwoComponentSPGRSignalPoints>**********************/
  // Sin and Cos of flip angles
  double sina, cosa;
  int flipAngle;
  
  // Holds final signal [# of SPGR flip angles]
  double theoreticalSPGRSignal[MAX_ALPHA_SPGR];   // Must define size at compile-time.
  
  // spgrRelaxationMinusIdentity = (e^(A*TR)-I)
  double spgrRelaxationMinusIdentity[2][2];
  
  // spgrIdentityMinusRelaxation = (I - e^(A*TR)*cos(alpha))
  double spgrIdentityMinusRelaxation[2][2];
  
  // spgrInvIdentityMinusRelaxation = inv(I - e^(A*TR)*cos(alpha))
  double spgrInvIdentityMinusRelaxation[2][2];

  // spgrRelaxationMinusIdentityTimesRecovery = (I-e^(A*TR)*M0*sin(alpha))
  double spgrRelaxationMinusIdentityTimesRecovery[2];

  for (i=0; i<2; i++) {
		for (j=0; j<2; j++) spgrRelaxationMinusIdentity[i][j] = spgrRelaxation[i][j] - I[i][j];
	}
  
  int N = d_nAlphaSPGR[0];
  
  // Loop over each flip angle & calculate signal
	for (flipAngle=0; flipAngle<N; flipAngle++) {
    // Compute sin & cos using fast device code
    sina = sin(d_alpha_spgr[flipAngle]*PI / 180.00);
    cosa = cos(d_alpha_spgr[flipAngle]*PI / 180.00);
		
		// calculate the identity minus relaxation * cosa matrix and the relaxation minus identity matrix
		for (i=0; i<2; i++) {
			for (j=0; j<2; j++) spgrIdentityMinusRelaxation[i][j] = I[i][j] - (spgrRelaxation[i][j]*cosa);
		}
		
		// calculate the inverse of the identityMinusRelaxation matrix
    // <------------generalMatrixInverse(spgrInvIdentityMinusRelaxation, spgrIdentityMinusRelaxation, 2)------------------>
    scalar = (1.00/ (spgrIdentityMinusRelaxation[0][0]*spgrIdentityMinusRelaxation[1][1] - spgrIdentityMinusRelaxation[1][0]*spgrIdentityMinusRelaxation[0][1]));
	
    spgrInvIdentityMinusRelaxation[0][0] = scalar * spgrIdentityMinusRelaxation[1][1];
    spgrInvIdentityMinusRelaxation[1][1] = scalar * spgrIdentityMinusRelaxation[0][0];
    spgrInvIdentityMinusRelaxation[1][0] = scalar * -1.00*spgrIdentityMinusRelaxation[1][0];
    spgrInvIdentityMinusRelaxation[0][1] = scalar * -1.00*spgrIdentityMinusRelaxation[0][1];
    //<------------------------End generalMatrixInverse------------------------->
		
		// calculate the product of the relaxationMinusIdentity x the recovery matrix
		spgrRelaxationMinusIdentityTimesRecovery[0] = (spgrRelaxationMinusIdentity[0][0]*spgrRecovery[0] + spgrRelaxationMinusIdentity[1][0]*spgrRecovery[1])*sina;
		spgrRelaxationMinusIdentityTimesRecovery[1] = (spgrRelaxationMinusIdentity[0][1]*spgrRecovery[0] + spgrRelaxationMinusIdentity[1][1]*spgrRecovery[1])*sina;
		
		// now, multiply invIdentityMinusRelaxation and relaxationMinusIdentityTimesRecovery together
		theoreticalSPGRSignal[flipAngle] = spgrInvIdentityMinusRelaxation[0][0]*spgrRelaxationMinusIdentityTimesRecovery[0] + spgrInvIdentityMinusRelaxation[1][0]*spgrRelaxationMinusIdentityTimesRecovery[1];
	  theoreticalSPGRSignal[flipAngle] += spgrInvIdentityMinusRelaxation[0][1]*spgrRelaxationMinusIdentityTimesRecovery[0] + spgrInvIdentityMinusRelaxation[1][1]*spgrRelaxationMinusIdentityTimesRecovery[1];
	}
  /******************</calculateAllTwoComponentSPGRSignalPoints********************/
  
	// calculate the mean signal intensity
  double spgrResiduals = 0.00;
		
  // calculate sum-of-square residuals.
	for (i=0; i<N; i++) {
    spgrResiduals += (SIGNALSCALE*theoreticalSPGRSignal[i] - SIGNALSCALE*d_data_spgr[i]) * (SIGNALSCALE*theoreticalSPGRSignal[i] - SIGNALSCALE*d_data_spgr[i]);
  }
  
  // Return Result!  Hoorah!
  d_resSPGR[0] = spgrResiduals;
  
} // </calcSPGR>


void calcSSFP(double* d_T1m, double* d_T1f, double* d_T2m, double* d_T2f, double* d_Fm, double* d_Tau_m, double* d_Omega, double* d_PD_ssfp, double* d_resSSFP, double* d_medSSFP, double rfPulsePhase) {
  
  // SSFP Signal Calculation 

  /***************<setUpTwoComponentSSFPRelaxationAndRecoveryMatrices>*************/
  double Aupper[4][4];
  double Alower[2][2];
	double w;
  int i, j, a;
  double sina, cosa;
  
  
  // copy over independant parameters
  // Map array paramter values onto Sean's variable naming
	double myelinT1 = d_T1m[0];
	double freeT1   = d_T1f[0];
	double myelinT2 = d_T2m[0];
	double freeT2   = d_T2f[0];
	double vMyelin  = d_Fm[0];
	double myelinResidenceTime = d_Tau_m[0];
    double freeWaterOffResonanceAngle = d_Omega[0];
    double ssfpProtonDensity = d_PD_ssfp[0];
  
  
	// calculate the dependent parameters
	double vFree = 1.00 - vMyelin;
	double exchangeMyelinToFree = (1.00 / myelinResidenceTime);
	double exchangeFreeToMyelin = ((vMyelin*exchangeMyelinToFree) / vFree);
  
  //spgr relaxation and recovery matricies
  double A[2][2];
  double invA[2][2];
  
  A[0][0] = (-1.00 / myelinT1) -  exchangeMyelinToFree;
	A[1][0] = exchangeFreeToMyelin;
	A[0][1] = exchangeMyelinToFree;
	A[1][1] = (-1.00 / freeT1) - exchangeFreeToMyelin;
  
  // compute Alower (using SSFP TR)
  // <---------specialMatrixExponential(Alower, A, ssfpTR, 2)----------------->

  for (i=0; i<2; i++){
    for (j=0; j<2; j++) Alower[i][j] = 0.00;
  }

  double I[2][2] = {1, 0, 0, 1};
  double scalarMatrix[2][2];
  double squareMatrix[2][2];
  double tripleMatrix[2][2];
  double quadMatrix[2][2];
  double quintMatrix[2][2];

  // first, apply the scalar to the incoming matrix
  for (i=0; i<2; i++) {
    for (j=0; j<2; j++) scalarMatrix[i][j] = d_tr_ssfp[0]*A[i][j];
  }


  // now calculate the square, triple, quad and quint matrices
  // start with the square		
  for (i=0; i<2; i++) {
    for (j=0; j<2; j++) {

      squareMatrix[i][j] = 0.00;
      for (a=0; a<2; a++) squareMatrix[i][j] += ((scalarMatrix[a][j]*scalarMatrix[i][a])/ 2.00);
    }
  }
  // next the triple
  for (i=0; i<2; i++) {
    for (j=0; j<2; j++) {

      tripleMatrix[i][j] = 0.00;
      for (a=0; a<2; a++) tripleMatrix[i][j] += ((squareMatrix[a][j]*scalarMatrix[i][a])/ 3.00);
    }
  }
  // next the quad
  for (i=0; i<2; i++) {
    for (j=0; j<2; j++) {

      quadMatrix[i][j] = 0.00;
      for (a=0; a<2; a++) quadMatrix[i][j] += ((tripleMatrix[a][j]*scalarMatrix[i][a])/ 4.00);
    }
  }
  // and, finally, the quint
  for (i=0; i<2; i++) {
    for (j=0; j<2; j++) {

      quintMatrix[i][j] = 0.00;
      for (a=0; a<2; a++) quintMatrix[i][j] += ((quadMatrix[a][j]*scalarMatrix[i][a])/ 5.00);
    }
  }

  // and, finally, sum them together appropriately
  for (i=0; i<2; i++) {
    for (j=0; j<2; j++) Alower[i][j] = I[i][j] + scalarMatrix[i][j] + squareMatrix[i][j] + tripleMatrix[i][j] + quadMatrix[i][j] + quintMatrix[i][j];
  }
  // <---------------------- End specialMatrixExponential ------------------------>
  
  // Compute inverse of A
	//<----------------------generalMatrixInverse(invA, A, 2)----------------------->
	double scalar = (1.00/ (A[0][0]*A[1][1] - A[1][0]*A[0][1]));
	
	invA[0][0] = scalar * A[1][1];
	invA[1][1] = scalar * A[0][0];
	invA[1][0] = scalar * -1.00*A[1][0];
	invA[0][1] = scalar * -1.00*A[0][1];
  //<------------------------End generalMatrixInverse------------------------->
  
  
  // ssfpRecovery = A^-1*C  where C is [0 0 0 0 Ff/T1f Fs/T1s]
  double ssfpRecovery[6] = {0, 0, 0, 0, 0, 0};
	ssfpRecovery[4] = (invA[0][0]*ssfpProtonDensity*vMyelin/ myelinT1) + (invA[1][0]*ssfpProtonDensity*vFree/ freeT1);
	ssfpRecovery[5] = (invA[0][1]*ssfpProtonDensity*vMyelin/ myelinT1) + (invA[1][1]*ssfpProtonDensity*vFree/ freeT1);
  
  // now, calculate the upper 4x4 portion of the relaxation matrix
  double A2[4][4];
  
  // free water off-resonance, assume 0 for simplified case
  w = (freeWaterOffResonanceAngle + rfPulsePhase)/d_tr_ssfp[0];
  
  A2[0][0] = (-1.00/ myelinT2) - exchangeMyelinToFree;
	A2[1][0] = exchangeFreeToMyelin;
	A2[2][0] = w;
  A2[3][0] = 0.00;
	
	A2[0][1] = exchangeMyelinToFree;
	A2[1][1] = (-1.00/ freeT2) - exchangeFreeToMyelin;
	A2[2][1] = 0.00;
  A2[3][1] = w;
  
	A2[0][2] = -w;
  A2[1][2] = 0.00;
	A2[2][2] = (-1.00/ myelinT2) - exchangeMyelinToFree;
	A2[3][2] = exchangeFreeToMyelin;
	
  A2[0][3] = 0.00;
	A2[1][3] = -w;
	A2[2][3] = exchangeMyelinToFree;
	A2[3][3] = (-1.00/ freeT2) - exchangeFreeToMyelin;
  
  
  // Compute Aupper
  scalar = d_tr_ssfp[0];
  double size = 4;
  //<--------------generalMatrixExponential(Aupper, A2, ssfpTR, 4);---------------->
  {
    double scalarMatrix[4][4], X[4][4], cX[4][4], E[4][4], D[4][4], tmp[4][4];
    double c;
    int i, j, k, q, p;
    int a, b;

    c = 0.5;
    p = 1;
    q = 6;

    // define the identity matrix
    double I[4][4] = {1, 0, 0, 0,
                     0, 1, 0, 0,
                     0, 0, 1, 0,
                     0, 0, 0, 1};

    // first, apply the scalar to the incoming matrix and the scale down by a large factor of 2 (256)
    for (i=0; i<size; i++) {
      for (j=0; j<size; j++) {
        scalarMatrix[i][j] = (scalar*A2[i][j]) / 256;
        X[i][j] = scalarMatrix[i][j];
        E[i][j] = I[i][j] + c*scalarMatrix[i][j];
        D[i][j] = I[i][j] - c*scalarMatrix[i][j];
      }
    }


    for (k=2; k<=q; k++) {
      c = c*(q-k+1.00) / (k*(2.00*q - k+1));

      // multiply A and X
      for (i=0; i<size; i++) {
        for (j=0; j<size; j++) {

          tmp[i][j] = 0.00;
          for (a=0; a<size; a++) tmp[i][j] += (X[a][j]*scalarMatrix[i][a]);
        }
      }
      for (i=0; i<size; i++) {
        for (j=0; j<size; j++) {
          X[i][j] = tmp[i][j];
          cX[i][j] = c*X[i][j];
          E[i][j] = E[i][j] + cX[i][j];

          if (p == 1) {
            D[i][j] = D[i][j] + cX[i][j];
          }
          else {
            D[i][j] = D[i][j] - cX[i][j];
          }

        }
      }

      p *= -1;

    }

    // now we need to calculate the inverse of D
    //<--------------------generalMatrixInverse(tmp, D, size);-------------------->
    {
      double maxValue, pointValue, tmpPoint;
      int pivot[4];
      int p, k, ptmp, pivotPoint;
      int i,j,a;

      for (i=0; i<size; i++) pivot[i] = i;

      for (i=0; i<size; i++) {
        for (j=0; j<size; j++) tmp[i][j] = 0.00;
      }

      for (k=0; k<size-1; k++) {
        p = k;
        maxValue = abs(D[k][k]);
        // scan the lower part of this column
        for (i=k; i<size; i++) {
          pointValue = abs(D[k][i]);
          if (pointValue > maxValue) {
            maxValue = pointValue;
            p = i;
          }
        }


        // swap rows if max value not in row k
        if (p != k) {
          for (a=0; a<size; a++) {
            tmpPoint = D[a][p];
            D[a][p] = D[a][k];
            D[a][k] = tmpPoint;
          }
          ptmp = pivot[k];
          pivot[k] = pivot[p];
          pivot[p] = ptmp;
        }


        // divide the lower trianglar part of column by max
        if (D[k][k] != 0) {
          for (i=k+1; i<size; i++) D[k][i] = D[k][i] / D[k][k];
        }
        // subtract multiple of column from remaining columes
        for (j=k+1; j<size; j++) {
          for (i=k+1; i<size; i++) {
            D[j][i] = D[j][i] - D[k][i]*D[j][k];
          }
        }

      }


      // now solve the inverse for each column
      for (k=0; k<size; k++) {
        pivotPoint = pivot[k];

        tmp[pivotPoint][k] = 1.00;

        for (j=k; j<size; j++) {
          if (tmp[pivotPoint][j] != 0) {
            for (i=j+1; i<size; i++) tmp[pivotPoint][i] = tmp[pivotPoint][i] - tmp[pivotPoint][j]*D[j][i];
          }
        }

        for (j=size-1; j>=0; j--) {
          tmp[pivotPoint][j] = tmp[pivotPoint][j] / D[j][j];
          if (tmp[pivotPoint][j] != 0) {
            for (i=0; i<j; i++) tmp[pivotPoint][i] = tmp[pivotPoint][i] - tmp[pivotPoint][j]*D[j][i];
          }
        }
      }
    }
    //<--------------------End generalMatrixInverse------------------------------->

    for (i=0; i<size; i++) {
      for (j=0; j<size; j++) {

        D[i][j] = 0.00;
        for (a=0; a<size; a++) D[i][j] += E[a][j]*tmp[i][a];
      }
    }


    // now, undo the initial scaling by repeated squaring - number of iterations = s where 256 = 2^s
    for (b=0; b<8; b++) {
      for (i=0; i<size; i++) {
        for (j=0; j<size; j++) {

          tmp[i][j] = 0.00;
          for (a=0; a<size; a++) tmp[i][j] += (D[a][j]*D[i][a]);
        }
      }
      for (i=0; i<size; i++) {
        for (j=0; j<size; j++) D[i][j] = tmp[i][j];
      }

    }

    // and, finally, sum them together appropriately
    for (i=0; i<size; i++) {
      for (j=0; j<size; j++) Aupper[i][j] = D[i][j];
    }
  }
  //<-----------------end generalMatrixExponential----------------------------->
  
  
  
  // Combine Aupper and Alower into ssfpRelaxation
  double ssfpRelaxation[6][6];
	for (i=0; i<4; i++) {
		for (j=0; j<4; j++) ssfpRelaxation[i][j] = Aupper[i][j];
	}
	for (i=0; i<2; i++) {
		for (j=0; j<2; j++) ssfpRelaxation[i+4][j+4] = Alower[i][j];
	}
  
  // Fill in the rest with zero values
  for (i=4; i<6 ;i++) {
   for (j=0; j<4 ;j++) {
     ssfpRelaxation[i][j] = 0.00;
     ssfpRelaxation[j][i] = 0.00;
   }
  }
  /***************</setUpTwoComponentSSFPRelaxationAndRecoveryMatrices>*************/
  
  
  /*****************<calculateAllTwoComponentSSFPSignalPoints>***********************/
	double rfRotation[6][6];
  // Initilize to zerossfp banding
  for (i=0; i<6 ;i++) {
    for (j=0; j<6 ;j++) {
      rfRotation[i][j] = 0.00;
    }
  }
  
	double ssfpMagnetization[6];
  double ssfpRelaxationMinusIdentity[6][6];
  double ssfpRelaxationTimesRotation[6][6];
  double ssfpIdentityMinusRelaxationTimesRotation[6][6];
  double ssfpRelaxationMinusIdentityTimesRecovery[6];
  double ssfpInvIdentityMinusRelaxationTimesRotation[6][6] = {1, 0, 0, 0, 0, 0,
                                                             0, 1, 0, 0, 0, 0,
                                                             0, 0, 1, 0, 0, 0,
                                                             0, 0, 0, 1, 0, 0,
                                                             0, 0, 0, 0, 1, 0,
                                                             0, 0, 0, 0, 0, 1};
  
  double theoreticalSSFPSignal[MAX_ALPHA_SSFP];   // Must define size at compile-time.
  
  // 6x6 Identity
  double I6[6][6] = {1, 0, 0, 0, 0, 0,
                    0, 1, 0, 0, 0, 0,
                    0, 0, 1, 0, 0, 0,
                    0, 0, 0, 1, 0, 0,
                    0, 0, 0, 0, 1, 0,
                    0, 0, 0, 0, 0, 1};
                   
  for (i=0; i<6; i++) {
		for (j=0; j<6; j++) ssfpRelaxationMinusIdentity[i][j] = ssfpRelaxation[i][j] - I6[i][j];
	}
                   

  // Look over each flip angle
  int N = d_nAlphaSSFP[0];
  
  int flipAngle;
  for (flipAngle=0; flipAngle<N; flipAngle++) {
        // Compute sin & cos using fast device code
    sina = sin(d_alpha_ssfp[flipAngle]*PI / 180.00);
    cosa = cos(d_alpha_ssfp[flipAngle]*PI / 180.00);
		
		// define the rf rotation matrix
		rfRotation[0][0] = 1.00;
		rfRotation[1][1] = 1.00;
		rfRotation[2][2] = cosa;
		rfRotation[3][3] = cosa;
		rfRotation[4][4] = cosa;
		rfRotation[5][5] = cosa;
		rfRotation[4][2] = sina;
		rfRotation[5][3] = sina;
		rfRotation[2][4] = -sina;
		rfRotation[3][5] = -sina;
		
		// calculate the relaxationTimesRotation matrix
		for (i=0; i<6; i++) {
			for (j=0; j<6; j++) ssfpRelaxationTimesRotation[j][i] = ssfpRelaxation[0][i]*rfRotation[j][0] + ssfpRelaxation[1][i]*rfRotation[j][1] + ssfpRelaxation[2][i]*rfRotation[j][2] + ssfpRelaxation[3][i]*rfRotation[j][3] + ssfpRelaxation[4][i]*rfRotation[j][4] + ssfpRelaxation[5][i]*rfRotation[j][5];
		}
		
		// calculate the identityMinusRelaxationTimesRotation and relaxationMinusIdentity matrices
		for (i=0; i<6; i++) {
			for (j=0; j<6; j++) ssfpIdentityMinusRelaxationTimesRotation[i][j] = I6[i][j]-ssfpRelaxationTimesRotation[i][j];
		}    
    
		//-----------------generalMatrixInverse(ssfpInvIdentityMinusRelaxationTimesRotation, ssfpIdentityMinusRelaxationTimesRotation, 6);-----------
    double maxValue, pointValue, tmpPoint;
    int pivot[6];
    int p, k, ptmp, pivotPoint;
    int size = 6;

    for (i=0; i<size; i++) pivot[i] = i;

    for (i=0; i<size; i++) {
      for (j=0; j<size; j++) ssfpInvIdentityMinusRelaxationTimesRotation[i][j] = 0.00;
    }


    for (k=0; k<size-1; k++) {
      p = k;
      maxValue = abs(ssfpIdentityMinusRelaxationTimesRotation[k][k]);
      // scan the lower part of this column
      for (i=k; i<size; i++) {
        pointValue = abs(ssfpIdentityMinusRelaxationTimesRotation[k][i]);
        if (pointValue > maxValue) {
				  maxValue = pointValue;
				  p = i;
        }
      }


      // swap rows if max value not in row k
      if (p != k) {
        for (a=0; a<size; a++) {
          tmpPoint = ssfpIdentityMinusRelaxationTimesRotation[a][p];
          ssfpIdentityMinusRelaxationTimesRotation[a][p] = ssfpIdentityMinusRelaxationTimesRotation[a][k];
          ssfpIdentityMinusRelaxationTimesRotation[a][k] = tmpPoint;
        }
        ptmp = pivot[k];
        pivot[k] = pivot[p];
        pivot[p] = ptmp;
      }


      // divide the lower trianglar part of column by max
      if (ssfpIdentityMinusRelaxationTimesRotation[k][k] != 0) {
        for (i=k+1; i<size; i++) ssfpIdentityMinusRelaxationTimesRotation[k][i] = ssfpIdentityMinusRelaxationTimesRotation[k][i] / ssfpIdentityMinusRelaxationTimesRotation[k][k];
      }
      // subtract multiple of column from remaining columes
      for (j=k+1; j<size; j++) {
        for (i=k+1; i<size; i++) {
          ssfpIdentityMinusRelaxationTimesRotation[j][i] = ssfpIdentityMinusRelaxationTimesRotation[j][i] - ssfpIdentityMinusRelaxationTimesRotation[k][i]*ssfpIdentityMinusRelaxationTimesRotation[j][k];
        }
      }
    }


    // now solve the inverse for each column
    for (k=0; k<size; k++) {
      pivotPoint = pivot[k];

      ssfpInvIdentityMinusRelaxationTimesRotation[pivotPoint][k] = 1.00;

      for (j=k; j<size; j++) {
        if (ssfpInvIdentityMinusRelaxationTimesRotation[pivotPoint][j] != 0) {
          for (i=j+1; i<size; i++) ssfpInvIdentityMinusRelaxationTimesRotation[pivotPoint][i] = ssfpInvIdentityMinusRelaxationTimesRotation[pivotPoint][i] - ssfpInvIdentityMinusRelaxationTimesRotation[pivotPoint][j]*ssfpIdentityMinusRelaxationTimesRotation[j][i];
        }
      }

      for (j=size-1; j>=0; j--) {
        ssfpInvIdentityMinusRelaxationTimesRotation[pivotPoint][j] = ssfpInvIdentityMinusRelaxationTimesRotation[pivotPoint][j] / ssfpIdentityMinusRelaxationTimesRotation[j][j];
        if (ssfpInvIdentityMinusRelaxationTimesRotation[pivotPoint][j] != 0) {
          for (i=0; i<j; i++) ssfpInvIdentityMinusRelaxationTimesRotation[pivotPoint][i] = ssfpInvIdentityMinusRelaxationTimesRotation[pivotPoint][i] - ssfpInvIdentityMinusRelaxationTimesRotation[pivotPoint][j]*ssfpIdentityMinusRelaxationTimesRotation[j][i];
        }

      }

    }
    //----------------------------------------End generalMatrixInverse --------------------------------------------------------------------------
		
		// calculate the relaxationMinusIdentityTimesRecovery matrix
		for (i=0; i<6; i++) ssfpRelaxationMinusIdentityTimesRecovery[i] = ssfpRelaxationMinusIdentity[0][i]*ssfpRecovery[0] + ssfpRelaxationMinusIdentity[1][i]*ssfpRecovery[1] + ssfpRelaxationMinusIdentity[2][i]*ssfpRecovery[2] + ssfpRelaxationMinusIdentity[3][i]*ssfpRecovery[3] + ssfpRelaxationMinusIdentity[4][i]*ssfpRecovery[4] + ssfpRelaxationMinusIdentity[5][i]*ssfpRecovery[5];
		
		// multiply invIdentityMinusRelaxationTimesRotation with relaxationMinusIdentityTimesRecovery to get the final magnetization
		for (i=0; i<6; i++) ssfpMagnetization[i] = ssfpInvIdentityMinusRelaxationTimesRotation[0][i]*ssfpRelaxationMinusIdentityTimesRecovery[0] + ssfpInvIdentityMinusRelaxationTimesRotation[1][i]*ssfpRelaxationMinusIdentityTimesRecovery[1]  + ssfpInvIdentityMinusRelaxationTimesRotation[2][i]*ssfpRelaxationMinusIdentityTimesRecovery[2]  + ssfpInvIdentityMinusRelaxationTimesRotation[3][i]*ssfpRelaxationMinusIdentityTimesRecovery[3]  + ssfpInvIdentityMinusRelaxationTimesRotation[4][i]*ssfpRelaxationMinusIdentityTimesRecovery[4]  + ssfpInvIdentityMinusRelaxationTimesRotation[5][i]*ssfpRelaxationMinusIdentityTimesRecovery[5];
		
    // Use fast square root round-to-nearest-even function
		theoreticalSSFPSignal[flipAngle] = sqrt(ssfpMagnetization[0]*ssfpMagnetization[0] + ssfpMagnetization[2]*ssfpMagnetization[2]) + sqrt( ssfpMagnetization[1]*ssfpMagnetization[1] + ssfpMagnetization[3]*ssfpMagnetization[3]);
	}
  /*****************</calculateAllTwoComponentSSFPSignalPoints>***********************/       
	
  double meanSignal    = 0.00;
  double ssfpResiduals = 0.00;
  
	for (i=0; i<N; i++) meanSignal += theoreticalSSFPSignal[i];
	meanSignal /= N;
  
  // determine the median value of the theoretical signal
	int middle;
  double tmp;
	middle = floor(N/2.00);
	
	double sortedArray[MAX_ALPHA_SSFP];
  
	for (i=0; i<N; i++) sortedArray[i] = theoreticalSSFPSignal[i];
	
	for (j=0; j<N; j++) {
		for (i=0; i<N-1; i++) {
			
			if (sortedArray[i] > sortedArray[i+1]) {
				tmp = sortedArray[i];
				sortedArray[i] = sortedArray[i+1];
				sortedArray[i+1] = tmp;
			}
		}
	}
		
  // compute sum of squared residuals
	for (i=0; i<N; i++) {
    if (rfPulsePhase == 0) {
      ssfpResiduals += (SIGNALSCALE*theoreticalSSFPSignal[i] - SIGNALSCALE*d_data_ssfp_0[i])   * (SIGNALSCALE*theoreticalSSFPSignal[i] - SIGNALSCALE*d_data_ssfp_0[i]);
    } else {
      ssfpResiduals += (SIGNALSCALE*theoreticalSSFPSignal[i] - SIGNALSCALE*d_data_ssfp_180[i]) * (SIGNALSCALE*theoreticalSSFPSignal[i] - SIGNALSCALE*d_data_ssfp_180[i]);
    }
  }
  
  // Return result
  d_resSSFP[0] = ssfpResiduals;
  d_medSSFP[0] = sortedArray[middle];
  
}


/* Host Helper Functions  ****************************************************/
void convert_double2double( double *input_double, double *output_double, int Ntot)
{
   int i;
   for (i = 0; i < Ntot; i++) {output_double[i] = (double) input_double[i];}
}


/*
 * Load in data from prhs[i] and convert to an array of double
 */
void load_rhs(double *output_double, int argin, const mxArray *prhs[]) {
  
  // Get a pointer to the double-doubleing point data
  double *input_double = mxGetPr(prhs[argin]);
  
  // Get the number of elements
  int nElem = mxGetM(prhs[argin]);
  
  // Convert elements to single-doubleing point, place in output array
  int i;
  for (i=0; i<nElem; i++) {output_double[i] = (double) input_double[i];}
}
