/***************************************************************************
 * NON CUDA-Based mcDESPOT Algorithm (cpMCDESPOT)
 *
 * v4.0  26-Oct-2009
 ***************************************************************************/

// Includes MEX for Matlab
#include <mex.h>

/* Function Prototypes ********************************************************/
// Main Function
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

// CPU Functions
void calcSPGR(double* d_T1m, double* d_T1f, double* d_T2m, double* d_T2f, double* d_Fm, double* d_Tau_m, double* d_PD_spgr, double* d_resSPGR);
void calcSSFP(double* d_T1m, double* d_T1f, double* d_T2m, double* d_T2f, double* d_Fm, double* d_Tau_m, double* d_Omega, double* d_PD_ssfp, double* d_resSSFP, double* d_medSSFP, double rfPulsePhase);

// Conversion Functions
void convert_double2double( double *input_double, double *output_double, int Ntot);

// Load in data from prhs
void load_rhs(double *output_double, int argin, const mxArray *prhs[]);
