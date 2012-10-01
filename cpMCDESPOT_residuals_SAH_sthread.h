/***************************************************************************
 * CPU-BASED mcDESPOT OBJECTIVE FUNCTION (cpMCDESPOT)
 *
 * v4.0  9-Nov-2010
 ***************************************************************************/

// Includes MEX for Matlab
#include <mex.h>

/* Function Prototypes ********************************************************/
// Main Function
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

// CPU Functions
void calcSPGR(double* fv, double* d_resSPGR);
void calcSSFP(double* fv, double* d_resSSFP, double rfPulsePhase);

// Helper Functions
void load_mrhs(double** output_double, int argin, const mxArray *prhs[]);
void matrixExponential(double TR, double Tm, double Tf, double Kmf, double Kfm, double MWF, double A[][2]);
void initZeros6by6(double matrix[][6]);
void matrixMultiply6by6(double X[][6], double Y[][6], double Z[][6]);

// Different ways of handling A\X=B problem
void matrixJXleastSquares6by6(double A[][6], double B[6], double X[6]);
void matrixGeneralInverse6by6(double A[][6], double B[6], double X[6]);
void matrixAnalyticalInverse6by6(double A[][6], double B[6], double X[6]);

// Debug Functions
void debugPrint2by2(double A[][2]);
void debugPrint2by1(double A[]);
void debugPrint6by6(double A[][6]);
void debugPrint6by1(double A[]);
