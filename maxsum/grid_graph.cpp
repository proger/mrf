#include <stdlib.h>
//#include "C:/Matlab/extern/include/mex.h"
#include "mex.h"

/*================================================================================================
  E = grid_graph(M,N)
================================================================================================*/
void mexFunction( int nargout, mxArray *argout[], int nargin, const mxArray *argin[] )
{
  if ( nargin!=1 ) mexErrMsgTxt("1 argument expected.");
  if ( !mxIsDouble(argin[0]) || mxGetM(argin[0])*mxGetN(argin[0])!=2 ) mexErrMsgTxt("The argument must be a double 3-vector.");
  double *dN= (double*)mxGetData(argin[0]);
  const unsigned M= (unsigned)dN[0], N= (unsigned)dN[1];

  unsigned nE= M*(N-1) + N*(M-1);
  const int dims[] = { 3, nE };
  argout[0] = mxCreateNumericArray( 2, dims, mxUINT32_CLASS, mxREAL );
  unsigned *E= (unsigned*)mxGetData(argout[0]), i, j;

  for ( j=0; j<N; j++ )
    for ( i=0; i<M-1; i++ ) {
      E[0]= i+0+M*j;
      E[1]= i+1+M*j;
      E[2]= 0;
      E+= 3;
    }
  for ( i=0; i<M; i++ )
    for ( j=0; j<N-1; j++ ) {
      E[0]= i+M*(j+0);
      E[1]= i+M*(j+1);
      E[2]= 1;
      E+= 3;
    }
}
