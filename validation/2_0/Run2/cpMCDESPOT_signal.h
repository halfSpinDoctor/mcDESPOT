/***************************************************************************
 * NON CUDA-Based mcDESPOT Signal Curves (cpMCDESPOT)
 *
 ***************************************************************************/

// Includes MEX for Matlab
#include <mex.h>

/* Function Prototypes ********************************************************/
// Main Function
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

// CPU Functions
void calcSPGR(float* d_T1m, float* d_T1f, float* d_T2m, float* d_T2f, float* d_Fm, float* d_Tau_m, float* d_resSPGR);
void calcSSFP(float* d_T1m, float* d_T1f, float* d_T2m, float* d_T2f, float* d_Fm, float* d_Tau_m, float* d_Omega, float* d_resSSFP, float rfPulsePhase);

// Conversion Functions
void convert_float2double( float *input_float, double *output_double, int Ntot);

// Load in data from prhs
void load_rhs(float *output_float, int argin, const mxArray *prhs[]);