/***************************************************************************
 * CPU-BASED DESPOT2-FM OBJECTIVE FUNCTION
 *
 * v1.0  15-Feb-2011
 ***************************************************************************/

// Includes MEX for Matlab
#include <mex.h>

/* Function Prototypes ********************************************************/
// Main Function
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

// CPU Functions
void calcDESPOT2(double* fv, double* d_res);

// Helper Functions
void load_mrhs(double** output_double, int argin, const mxArray *prhs[]);
