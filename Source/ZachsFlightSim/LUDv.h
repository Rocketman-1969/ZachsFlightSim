#include <cmath>
#define ZERO 1.0e-12
using namespace std;

/*------------------------- Function Prototypes -----------------------------*/

double LUDecomp( double** a, int n, double* pvt);
void LUSolve( double** a, int n, double* b,double* pvt);
