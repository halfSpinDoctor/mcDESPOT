/***************************************************************************
 * NON CUDA-Based mcDESPOT Algorithm (cpMCDESPOT)
 *
 * Returns shape of SPGR or SSFP signal curve
 *
 * [resSPGR resSSFP_0 resSSFP_180] = mcDESPOT_signal(T1m, T1f, T2m, T2f, Fm, Tau_m, Omega, alpha_spgr, alpha_ssfp, tr_spgr, tr_ssfp)
 *
 * T1m, T1f, T2m, T2f, Fm, Tau_m, Omega --> [Nx1]  Parameter Vectors, where N=nParam (number of parameter trials)
 *
 * data_spgr, data_ssfp_0 data_ssfp_180 --> [Npx1] MRI Data, # of flip angles by 1
 *
 * alpha_spgr, alpha_ssfp        --> [scalar] Flip angles (corrected /w B1err) in degrees
 *
 * tr_spgr, tr_ssfp              --> [scalar] TR times, in ms
 *
 * COMPILE COMMAND: mex CFLAGS="\$CFLAGS -std=c99" -lm cpMCDESPOT_signal.c
 *
 *
 * Samuel A. Hurley
 * Pouria Mossahebi
 * University of Wisconsin
 * v1.0  25-Mar-2010
 *
 * Changelog:
 *    v1.0 - initial code, based on cpMCDESPOT_residual (Mar-2010)
 ***************************************************************************/

/* Includes MEX for Matlab and CUDA headers */
#include <mex.h>
#include <math.h>

// Application Header
#include "cpMCDESPOT_signal.h"

// Define maximum number of SPGR and SSFP data points we can have,
// in order to use constant memory effectively
#define MAX_ALPHA_SPGR 40
#define MAX_ALPHA_SSFP 40

// Everyone's favourite number...
#define PI 3.1415926535897932384626433
// Re-scale signal before taking difference for residuals
#define SIGNALSCALE 1000.00

// Global variables
float d_data_spgr[MAX_ALPHA_SPGR];
float d_data_ssfp_0[MAX_ALPHA_SSFP];
float d_data_ssfp_180[MAX_ALPHA_SSFP];

float d_alpha_spgr[MAX_ALPHA_SPGR];
float d_alpha_ssfp[MAX_ALPHA_SSFP];

float alpha_spgr[MAX_ALPHA_SPGR];
float alpha_ssfp[MAX_ALPHA_SSFP];

float d_tr_spgr[1];
float d_tr_ssfp[1];
int   d_nAlphaSPGR[1];
int   d_nAlphaSSFP[1];

/* Main Function  ************************************************************/
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  
  /**** 0. Variable Declrations **********************************/
  int i;
  int nParam, nAlphaSPGR, nAlphaSSFP;
  
  // Size of data vectors in memory
  size_t paramSize;
  size_t flipOutputSize;
  size_t spgrSize;
  size_t ssfpSize;
  
  // --- Host Data --- //
  // Parameter vector
  float *T1m, *T1f, *T2m, *T2f, *Fm, *Tau_m, *Omega;
  
  // Residual solutions
  float *resSPGR, *resSSFP_0, *resSSFP_180;
  
  // MR Data
  float *tr_spgr,    *tr_ssfp;
  
  // --- Device Data (deonte with d_) --- //
  // Parameter vector
  // float *d_T1m, *d_T1f, *d_T2m, *d_T2f, *d_Fm, *d_Tau_m, *d_Omega;
  
  // Residual Solutions
  // float *d_resSPGR, *d_resSSFP_0, *d_resSSFP_180;
  // float *d_medSSFP_0, *d_medSSFP_180;
  
  // MR Data -- Declared as global constant memory (see top of this file)
  // float *d_data_spgr,  *d_data_ssfp_0;
  //                      *d_data_ssfp_180;
  // float *d_alpha_spgr, *d_alpha_ssfp;
  // float *d_tr_spgr,    *d_tr_ssfp;
  
  // CUDA-Device Error
  // cudaError_t cError = cudaSuccess;
  
  
  /* I. Error checking ****************************************************************/
  // Check for a single output
  if (nlhs != 3)
    mexErrMsgTxt("Requires three outputs [sigSPGR sigSSFP_0 sigSSFP_180 ]");
  
  // Check for 11 inputs
  if (nrhs != 11)
    mexErrMsgTxt("Requires 11 inputs: [T1m, T1f, T2m, T2f, Fm, Tau_m, Omega, alpha_spgr, alpha_ssfp, tr_spgr, tr_ssfp]");
  
  
  // Check that all inputs are Nx1 (single column)
  for (i=0; i<11; i++) {
   if (mxGetN(prhs[i]) > 1)
     mexErrMsgTxt("All inputs must be Nx1 (one column wide)");
  }
  
  // Grab nParam
  nParam = mxGetM(prhs[0]);
  
  // Check that first seven inputs all have same length
  for (i=0; i<7; i++) {
   if (mxGetM(prhs[i]) != nParam)
     mexErrMsgTxt("Parameter inputs (#1-7) must have same length");
  }
  
  // Grab the number of SPGR & SSFP flip angles
  nAlphaSPGR = mxGetM(prhs[7]);
  nAlphaSSFP = mxGetM(prhs[8]);
  
  // Check that the TR's are scalars
  if (mxGetN(prhs[9]) !=1 || mxGetM(prhs[9]) !=1)
    mexErrMsgTxt("SPGR TR must be a scalar (in milliseconds)");
  
  if (mxGetN(prhs[10]) !=1 || mxGetM(prhs[10]) !=1)
    mexErrMsgTxt("SSFP TR must be a scalar (in milliseconds)");
  
  /* II. Load In Data *******************************************************/
  
  // Size of param vectors
  paramSize = sizeof(float)*nParam;
  
  // Alloc param vectors on host
  T1m   = (float*) mxMalloc(paramSize);
  T1f   = (float*) mxMalloc(paramSize);
  T2m   = (float*) mxMalloc(paramSize);
  T2f   = (float*) mxMalloc(paramSize);
  Fm    = (float*) mxMalloc(paramSize);
  Tau_m = (float*) mxMalloc(paramSize);
  Omega = (float*) mxMalloc(paramSize);
  
  /*
  cudaSetDeviceFlags(cudaDeviceMapHost);
  cError = cudaHostAlloc((void **)&resSPGR, paramSize, cudaHostAllocMapped);
  cError = cudaHostAlloc((void **)&resSSFP_0, paramSize, cudaHostAllocMapped);
  cError = cudaHostAlloc((void **)&medSSFP_0, paramSize, cudaHostAllocMapped);
  cError = cudaHostAlloc((void **)&resSSFP_180, paramSize, cudaHostAllocMapped);
  cError = cudaHostAlloc((void **)&medSSFP_180, paramSize, cudaHostAllocMapped);
  */
  
  flipOutputSize = sizeof(float)*nAlphaSPGR;
  
  if (nAlphaSPGR != nAlphaSSFP)
    mexErrMsgTxt("Must use same number of SPGR and SSFP flip angles for signal calculation");

  resSPGR     = (float*) mxMalloc(flipOutputSize);
  resSSFP_0   = (float*) mxMalloc(flipOutputSize);
  resSSFP_180 = (float*) mxMalloc(flipOutputSize);

  
  // Alloc param vectors on GPU device
  /*
  cError = cudaMalloc((void**)&d_T1m,   paramSize);
  cError = cudaMalloc((void**)&d_T1f,   paramSize);
  cError = cudaMalloc((void**)&d_T2m,   paramSize);
  cError = cudaMalloc((void**)&d_T2f,   paramSize);
  cError = cudaMalloc((void**)&d_Fm,    paramSize);
  cError = cudaMalloc((void**)&d_Tau_m, paramSize);
  cError = cudaMalloc((void**)&d_Omega, paramSize);
  */
  
  // OLD -- for page-locked memory
  // Results: Grab a pointer to the device memory
  /*
  cError = cudaHostGetDevicePointer((void **)&d_resSPGR, (void *)resSPGR, 0);
  cError = cudaHostGetDevicePointer((void **)&d_resSSFP_0, (void *)resSSFP_0, 0);
  cError = cudaHostGetDevicePointer((void **)&d_medSSFP_0, (void *)medSSFP_0, 0);
  cError = cudaHostGetDevicePointer((void **)&d_resSSFP_180, (void *)resSSFP_180, 0);
  cError = cudaHostGetDevicePointer((void **)&d_medSSFP_180, (void *)medSSFP_180, 0);
  */
  /*
  cError = cudaMalloc((void**)&d_resSPGR,     paramSize);
  cError = cudaMalloc((void**)&d_resSSFP_0,   paramSize);
  cError = cudaMalloc((void**)&d_resSSFP_180, paramSize);
  cError = cudaMalloc((void**)&d_medSSFP_0  , paramSize);
  cError = cudaMalloc((void**)&d_medSSFP_180, paramSize);
  */
  
  // Size of SPGR & SSFP Data
  spgrSize = sizeof(float)*nAlphaSPGR;
  ssfpSize = sizeof(float)*nAlphaSSFP;
  
  tr_spgr = (float*) mxMalloc(sizeof(float));  // Scalar value
  tr_ssfp = (float*) mxMalloc(sizeof(float));  // Scalar value
  
  // Constant memory already reserved on host, do not need cudaMalloc
  
  // Check for errors performing cudaAlloc
  /*
  if(cError != cudaSuccess) {
    mexPrintf("CUDA Error Occured : %s\n", cudaGetErrorString(cudaGetLastError()));
    mexErrMsgTxt("CUDA Fatal: Error occured attempting to allocate memory on device");
	}*/
  
  // Load in parameter data
  load_rhs(T1m,   0, prhs);
  load_rhs(T1f,   1, prhs);
  load_rhs(T2m,   2, prhs);
  load_rhs(T2f,   3, prhs);
  load_rhs(Fm,    4, prhs);
  load_rhs(Tau_m, 5, prhs);
  load_rhs(Omega, 6, prhs);
  
  // Load in MRI data
  load_rhs(alpha_spgr,   7, prhs);
  load_rhs(alpha_ssfp,   8, prhs);
  load_rhs(tr_spgr,      9, prhs);
  load_rhs(tr_ssfp,     10, prhs);
  
  // -- GPU Device: Copy Parameter data --
  // Reset Error Flag
  /*
  cError = cudaSuccess;
  
  cError = cudaMemcpy(d_T1m,   T1m,   paramSize, cudaMemcpyHostToDevice);
  cError = cudaMemcpy(d_T1f,   T1f,   paramSize, cudaMemcpyHostToDevice);
  cError = cudaMemcpy(d_T2m,   T2m,   paramSize, cudaMemcpyHostToDevice);
  cError = cudaMemcpy(d_T2f,   T2f,   paramSize, cudaMemcpyHostToDevice);
  cError = cudaMemcpy(d_Fm,    Fm,    paramSize, cudaMemcpyHostToDevice);
  cError = cudaMemcpy(d_Tau_m, Tau_m, paramSize, cudaMemcpyHostToDevice);
  cError = cudaMemcpy(d_Omega, Omega, paramSize, cudaMemcpyHostToDevice);
  
  // Check for errors performing cudaMemcpy
  if(cError != cudaSuccess) {
    mexPrintf("CUDA Error Occured : %s\n", cudaGetErrorString(cudaGetLastError()));
    mexErrMsgTxt("CUDA Fatal: Error occured attempting to copy host->device dynamic memory (cudaMemcpy)");
	}
  
  // -- GPU Device: Copy MRI data to constant memory --
  // Reset error flag
  cError = cudaSuccess;
   **/
  
  // Check that the size of spgr/ssfp data are not larger than the hard-coded maxima
  if (spgrSize > sizeof(float)*MAX_ALPHA_SPGR) {
    mexErrMsgTxt("Error: This algorithm has a hard-coded limit of SPGR data points, and you used too many!\nDo not panic, the local authorities are on their way...");
  }
  
  if (ssfpSize > sizeof(float)*MAX_ALPHA_SSFP) {
    mexErrMsgTxt("Error: This algorithm has a hard-coded limit of SSFP data points, and you used too many!\nDo not panic, the local authorities are on their way...");
  }
  
  /*
  
  cError = cudaMemcpyToSymbol(d_data_spgr, data_spgr, spgrSize, 0, cudaMemcpyHostToDevice);
  cError = cudaMemcpyToSymbol(d_data_ssfp_0  , data_ssfp_0,   ssfpSize, 0, cudaMemcpyHostToDevice);
  cError = cudaMemcpyToSymbol(d_data_ssfp_180, data_ssfp_180, ssfpSize, 0, cudaMemcpyHostToDevice);
  
  cError = cudaMemcpyToSymbol(d_alpha_spgr, alpha_spgr, spgrSize, 0, cudaMemcpyHostToDevice);
  cError = cudaMemcpyToSymbol(d_alpha_ssfp, alpha_ssfp, ssfpSize, 0, cudaMemcpyHostToDevice);
  
  cError = cudaMemcpyToSymbol(d_tr_spgr, tr_spgr, sizeof(float), 0, cudaMemcpyHostToDevice);
  cError = cudaMemcpyToSymbol(d_tr_ssfp, tr_ssfp, sizeof(float), 0, cudaMemcpyHostToDevice);
  
  cError = cudaMemcpyToSymbol(d_nAlphaSPGR, &nAlphaSPGR, sizeof(int), 0, cudaMemcpyHostToDevice);
  cError = cudaMemcpyToSymbol(d_nAlphaSSFP, &nAlphaSSFP, sizeof(int), 0, cudaMemcpyHostToDevice);
  */
  
  // Copy SSFP & SPGR data into constant memory -- probably a more efficient way to do this but i just don't care
  d_nAlphaSPGR[0] = nAlphaSPGR;
  d_nAlphaSSFP[0] = nAlphaSSFP;
  
  d_tr_spgr[0] = (float) tr_spgr[0];
  d_tr_ssfp[0] = (float) tr_ssfp[0];

  
  /* III. Data Processing *****************************************************/
  // Call CPU-based functions
  
  // Reset error codes
  // cError = cudaSuccess;
  
  i = 0;
  
  // SPGR
  calcSPGR(&T1m[i], &T1f[i], &T2m[i], &T2f[i], &Fm[i], &Tau_m[i], &resSPGR[i]);
  
  // SSFP with Off-Resonance Correction
  calcSSFP(&T1m[i], &T1f[i], &T2m[i], &T2f[i], &Fm[i], &Tau_m[i], &Omega[i], &resSSFP_0[i],   0);
  calcSSFP(&T1m[i], &T1f[i], &T2m[i], &T2f[i], &Fm[i], &Tau_m[i], &Omega[i], &resSSFP_180[i], PI);
  
  /* IV. Place results in output vector *************/
  plhs[0] = mxCreateDoubleMatrix(nAlphaSPGR,1,mxREAL);
  plhs[1] = mxCreateDoubleMatrix(nAlphaSPGR,1,mxREAL);
  plhs[2] = mxCreateDoubleMatrix(nAlphaSPGR,1,mxREAL);
  
  convert_float2double(resSPGR    ,  mxGetPr(plhs[0]),  nAlphaSPGR);
  convert_float2double(resSSFP_0  ,  mxGetPr(plhs[1]),  nAlphaSPGR);
  convert_float2double(resSSFP_180,  mxGetPr(plhs[2]),  nAlphaSPGR);
  
  
  /* V. Cleanup & Free Memory ************************************************/
  // Parameter vectors - host
  mxFree(T1m);
  mxFree(T1f);
  mxFree(T2m);
  mxFree(T2f);
  mxFree(Fm);
  mxFree(Tau_m);
  mxFree(Omega);
  
  // Output vectors
  mxFree(resSPGR);
  mxFree(resSSFP_0  );
  mxFree(resSSFP_180);
  
  // MRI Data - host
  mxFree(tr_spgr);
  mxFree(tr_ssfp);
  
} // </Main mexFunction>

/* CPU Functions  ***********************************************************/

// Kernel to initiate residual calcuations
void calcSPGR(float* d_T1m, float* d_T1f, float* d_T2m, float* d_T2f, float* d_Fm, float* d_Tau_m, float* d_resSPGR) {
  
  /* Declare Common Variables used within this kernel **********/
  int i, j;
  int a;
  
  // SPGR Signal Calculation
  
  /************<setUpTwoComponentSPGRRelaxationAndRecoveryMatrices>**************/
  
  // Map array paramter values onto Sean's variable naming
	float myelinT1 = d_T1m[0];
	float freeT1   = d_T1f[0];
	// float myelinT2 = d_T2m[idx];  Not required for SPGR
	// float freeT2   = d_T2f[idx];  Not required for SPGR
	float vMyelin  = d_Fm[0];
	float myelinResidenceTime = d_Tau_m[0];
  
	// calculate the dependent parameters
	float vFree = 1.00 - vMyelin;
	float exchangeMyelinToFree = (1.00/ myelinResidenceTime);
	float exchangeFreeToMyelin = ((vMyelin*exchangeMyelinToFree)/ vFree);
  
  //spgr relaxation and recovery matricies
  float A[2][2];
  float invA[2][2] = {0, 0, 0, 0};  // Initilize to zero
  
  A[0][0] = (-1.00/ myelinT1) -  exchangeMyelinToFree;
	A[1][0] = exchangeFreeToMyelin;
	A[0][1] = exchangeMyelinToFree;
	A[1][1] = (-1.00/ freeT1) - exchangeFreeToMyelin;
  
  // Create spgrRelaxation = e^A*TR
  float spgrRelaxation[2][2];
  
  // <---------specialMatrixExponential(spgrRelaxation, A, spgrTR, 2)----------------->

  for (i=0; i<2; i++){
    for (j=0; j<2; j++) spgrRelaxation[i][j] = 0.00;
  }

  float I[2][2];
  float scalarMatrix[2][2];
  float squareMatrix[2][2];
  float tripleMatrix[2][2];
  float quadMatrix[2][2];
  float quintMatrix[2][2];

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
	float scalar = (1.00/ (A[0][0]*A[1][1] - A[1][0]*A[0][1]));
	
	invA[0][0] = scalar * A[1][1];
	invA[1][1] = scalar * A[0][0];
	invA[1][0] = scalar * -1.00*A[1][0];
	invA[0][1] = scalar * -1.00*A[0][1];
  //<------------------------End generalMatrixInverse------------------------->
  
  
  // now, tackle the recovery matricies
  // spgrRecovery = inv(A)*C = M0
  float spgrRecovery[2];
	spgrRecovery[0] = (invA[0][0]*vMyelin/ myelinT1) + (invA[1][0]*vFree/ freeT1);
	spgrRecovery[1] = (invA[0][1]*vMyelin/ myelinT1) + (invA[1][1]*vFree/ freeT1);
  
  /***********</SetupTwoComponantSPGRRelaxationAndRecoveryMatricies>*********************/
  
  
  
  /********************<calculateAllTwoComponentSPGRSignalPoints>**********************/
  // Sin and Cos of flip angles
  float sina, cosa;
  int flipAngle;
  
  // Holds final signal [# of SPGR flip angles]
  float theoreticalSPGRSignal[MAX_ALPHA_SPGR];   // Must define size at compile-time.
  
  // spgrRelaxationMinusIdentity = (e^(A*TR)-I)
  float spgrRelaxationMinusIdentity[2][2];
  
  // spgrIdentityMinusRelaxation = (I - e^(A*TR)*cos(alpha))
  float spgrIdentityMinusRelaxation[2][2];
  
  // spgrInvIdentityMinusRelaxation = inv(I - e^(A*TR)*cos(alpha))
  float spgrInvIdentityMinusRelaxation[2][2];

  // spgrRelaxationMinusIdentityTimesRecovery = (I-e^(A*TR)*M0*sin(alpha))
  float spgrRelaxationMinusIdentityTimesRecovery[2];

  for (i=0; i<2; i++) {
		for (j=0; j<2; j++) spgrRelaxationMinusIdentity[i][j] = spgrRelaxation[i][j] - I[i][j];
	}
  
  
  // Loop over each flip angle & calculate signal
	for (flipAngle=0; flipAngle<d_nAlphaSPGR[0]; flipAngle++) {
    // Compute sin & cos using fast device code
    sina = sin(alpha_spgr[flipAngle]*PI / 180.00);
    cosa = cos(alpha_spgr[flipAngle]*PI / 180.00);
		
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
    
    // Return Result!  Hoorah!
    d_resSPGR[flipAngle] = theoreticalSPGRSignal[flipAngle];
	}
  /******************</calculateAllTwoComponentSPGRSignalPoints********************/
  
	// calculate the mean signal intensity
	// float meanSignal = 0.00;
  // float spgrResiduals = 0.00;
  
  // find mean signal
	// for (i=0; i<N; i++) meanSignal += theoreticalSPGRSignal[i];
	// meanSignal /= N;
		
	// divide the true signal values by the mean
  // calculate sum-of-square residuals.  asumes that data_spgr has already been normalized
	// for (i=0; i<N; i++) {
  //   theoreticalSPGRSignal[i] /= meanSignal;
  //   spgrResiduals += (SIGNALSCALE*theoreticalSPGRSignal[i] - SIGNALSCALE*d_data_spgr[i]) * (SIGNALSCALE*theoreticalSPGRSignal[i] - SIGNALSCALE*d_data_spgr[i]);
  // }
  
  
} // </calcSPGR>


void calcSSFP(float* d_T1m, float* d_T1f, float* d_T2m, float* d_T2f, float* d_Fm, float* d_Tau_m, float* d_Omega, float* d_resSSFP, float rfPulsePhase) {
  
  // SSFP Signal Calculation 

  /***************<setUpTwoComponentSSFPRelaxationAndRecoveryMatrices>*************/
  float Aupper[4][4];
  float Alower[2][2];
	float w;
  int i, j, a;
  float sina, cosa;
  
  
  // copy over independant parameters
  // Map array paramter values onto Sean's variable naming
	float myelinT1 = d_T1m[0];
	float freeT1   = d_T1f[0];
	float myelinT2 = d_T2m[0];
	float freeT2   = d_T2f[0];
	float vMyelin  = d_Fm[0];
	float myelinResidenceTime = d_Tau_m[0];
  float freeWaterOffResonanceAngle = d_Omega[0];
  
  
	// calculate the dependent parameters
	float vFree = 1.00 - vMyelin;
	float exchangeMyelinToFree = (1.00 / myelinResidenceTime);
	float exchangeFreeToMyelin = ((vMyelin*exchangeMyelinToFree) / vFree);
  
  //spgr relaxation and recovery matricies
  float A[2][2];
  float invA[2][2];
  
  A[0][0] = (-1.00 / myelinT1) -  exchangeMyelinToFree;
	A[1][0] = exchangeFreeToMyelin;
	A[0][1] = exchangeMyelinToFree;
	A[1][1] = (-1.00 / freeT1) - exchangeFreeToMyelin;
  
  // compute Alower (using SSFP TR)
  // <---------specialMatrixExponential(Alower, A, ssfpTR, 2)----------------->

  for (i=0; i<2; i++){
    for (j=0; j<2; j++) Alower[i][j] = 0.00;
  }

  float I[2][2] = {1, 0, 0, 1};
  float scalarMatrix[2][2];
  float squareMatrix[2][2];
  float tripleMatrix[2][2];
  float quadMatrix[2][2];
  float quintMatrix[2][2];

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
	float scalar = (1.00/ (A[0][0]*A[1][1] - A[1][0]*A[0][1]));
	
	invA[0][0] = scalar * A[1][1];
	invA[1][1] = scalar * A[0][0];
	invA[1][0] = scalar * -1.00*A[1][0];
	invA[0][1] = scalar * -1.00*A[0][1];
  //<------------------------End generalMatrixInverse------------------------->
  
  
  // ssfpRecovery = A^-1*C  where C is [0 0 0 0 Ff/T1f Fs/T1s]
  float ssfpRecovery[6] = {0, 0, 0, 0, 0, 0};
	ssfpRecovery[4] = (invA[0][0]*vMyelin/ myelinT1) + (invA[1][0]*vFree/ freeT1);
	ssfpRecovery[5] = (invA[0][1]*vMyelin/ myelinT1) + (invA[1][1]*vFree/ freeT1);
  
  // now, calculate the upper 4x4 portion of the relaxation matrix
  float A2[4][4];
  
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
  int size = 4;
  //<--------------generalMatrixExponential(Aupper, A2, ssfpTR, 4);---------------->
  {
    float scalarMatrix[4][4], X[4][4], cX[4][4], E[4][4], D[4][4], tmp[4][4];
    float c;
    int i, ii, j, jj, k, q, p, m, n;
    int a, b;

    c = 0.5;
    p = 1;
    q = 6;

    // define the identity matrix
    float I[4][4] = {1, 0, 0, 0,
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
    //<--------------------jxMatrixInverse(tmp, D, size);-------------------->
    // Pouria's new Jordan Exchange Inverse Implementation
    // Input is D, output is tmp
    
    {
      // Compute inverse
      // Change Aupper to input matrix, invAupper to output matrix
      for (i=0; i<size; i++) {
        for (ii=0; ii<size; ii++) {

          // Precompute x
          float x = D[i][ii]/D[i][i];

          for (jj=0; jj<size; jj++) {
            if (ii == i && jj == i)      // Pivot
              tmp[i][i] = 1.00 / D[i][i];
            else if (ii == i && jj != i) // Row
              tmp[jj][i] = -D[jj][i]/D[i][i];
            else if (ii != i && jj == i) // Column
              tmp[i][ii] = x;
            else                       // Off-diag
              tmp[jj][ii] = D[jj][ii] - x * D[jj][i];
          }
        }

        // Copy result back to Aupper
        // Note that this will mess /w Aupper, so we can't use it in subsequent calculations
        for (m=0; m<size; m++) {
          for (n=0; n<size; n++) {
            D[m][n] = tmp[m][n];
          }
        }
      }
    }
                     
    //<--------------------End jxMatrixInverse------------------------------->

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
  float ssfpRelaxation[6][6];
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
	float rfRotation[6][6];
  // Initilize to zerossfp banding
  for (i=0; i<6 ;i++) {
    for (j=0; j<6 ;j++) {
      rfRotation[i][j] = 0.00;
    }
  }
  
	float ssfpMagnetization[6];
  float ssfpRelaxationMinusIdentity[6][6];
  float ssfpRelaxationTimesRotation[6][6];
  float ssfpIdentityMinusRelaxationTimesRotation[6][6];
  float ssfpRelaxationMinusIdentityTimesRecovery[6];
  float ssfpInvIdentityMinusRelaxationTimesRotation[6][6] = {1, 0, 0, 0, 0, 0,
                                                             0, 1, 0, 0, 0, 0,
                                                             0, 0, 1, 0, 0, 0,
                                                             0, 0, 0, 1, 0, 0,
                                                             0, 0, 0, 0, 1, 0,
                                                             0, 0, 0, 0, 0, 1};
  
  float theoreticalSSFPSignal[MAX_ALPHA_SSFP];   // Must define size at compile-time.
  
  // 6x6 Identity
  float I6[6][6] = {1, 0, 0, 0, 0, 0,
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
    sina = sin(alpha_ssfp[flipAngle]*PI / 180.00);
    cosa = cos(alpha_ssfp[flipAngle]*PI / 180.00);
		
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
    float maxValue, pointValue, tmpPoint;
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
    
    // Return Signal
    d_resSSFP[flipAngle] = theoreticalSSFPSignal[flipAngle];
    
	}
  /*****************</calculateAllTwoComponentSSFPSignalPoints>***********************/       
	
	// calculate the mean signal intensity
	// float meanSignal = 0.00;
  // float ssfpResiduals = 0.00;
  
	// for (i=0; i<N; i++) meanSignal += theoreticalSSFPSignal[i];
	// meanSignal /= N;
  
  // determine the median value of the theoretical signal
	// int middle;
  // float tmp;
	// middle = floor(N/2.00);
	
	// float sortedArray[MAX_ALPHA_SSFP];
  
	// for (i=0; i<N; i++) sortedArray[i] = theoreticalSSFPSignal[i];
	
	// for (j=0; j<N; j++) {
		// for (i=0; i<N-1; i++) {
			
			// if (sortedArray[i] > sortedArray[i+1]) {
				// tmp = sortedArray[i];
				// sortedArray[i] = sortedArray[i+1];
				// sortedArray[i+1] = tmp;
			// }
		// }
	// }
	
		
	// divide the true signal values by the mean
  // compute sum of squared residuals
	// for (i=0; i<N; i++) {
    // theoreticalSSFPSignal[i] /= meanSignal;
    // if (rfPulsePhase == 0) {
      // ssfpResiduals += (SIGNALSCALE*theoreticalSSFPSignal[i] - SIGNALSCALE*d_data_ssfp_0[i])   * (SIGNALSCALE*theoreticalSSFPSignal[i] - SIGNALSCALE*d_data_ssfp_0[i]);
    // } else {
      // ssfpResiduals += (SIGNALSCALE*theoreticalSSFPSignal[i] - SIGNALSCALE*d_data_ssfp_180[i]) * (SIGNALSCALE*theoreticalSSFPSignal[i] - SIGNALSCALE*d_data_ssfp_180[i]);
    // }
  // }
  
}


/* Host Helper Functions  ****************************************************/
void convert_float2double( float *input_float, double *output_double, int Ntot)
{
   int i;
   for (i = 0; i < Ntot; i++) {output_double[i] = (double) input_float[i];}
}


/*
 * Load in data from prhs[i] and convert to an array of double
 */
void load_rhs(float *output_float, int argin, const mxArray *prhs[]) {
  
  // Get a pointer to the double-floating point data
  double *input_double = mxGetPr(prhs[argin]);
  
  // Get the number of elements
  int nElem = mxGetM(prhs[argin]);
  
  // Convert elements to single-floating point, place in output array
  int i;
  for (i=0; i<nElem; i++) {output_float[i] = (float) input_double[i];}
}
