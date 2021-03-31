/*
 * random_JB.cpp
 *
 * Function generating uniform and gaussian random numbers.
 * Initially created by John Buckart, but modified by myself
 * in order to be more compliant with current standards
 *
 *  Created on: 18 Mar 2016
 *      Author: obenomar
 */
#include <iostream>
#include <iomanip>
# include <Eigen/Dense>

using namespace std;

double *r8vec_normal_01 ( int n, int *seed );
double *uniform_01 ( int n, int *seed );

double *r8vec_normal_01 ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_NORMAL_01 returns a unit pseudonormal R8VEC.
//
//  Discussion:
//
//    The standard normal probability distribution function (PDF) has
//    mean 0 and standard deviation 1.
//
//    This routine can generate a vector of values on one call.  It
//    has the feature that it should provide the same results
//    in the same order no matter how we break up the task.
//
//    Before calling this routine, the user may call RANDOM_SEED
//    in order to set the seed of the random number generator.
//
//    The Box-Muller method is used, which is efficient, but
//    generates an even number of values each time.  On any call
//    to this routine, an even number of new values are generated.
//    Depending on the situation, one value may be left over.
//    In that case, it is saved for the next call.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of values desired.  If N is negative,
//    then the code will flush its internal memory; in particular,
//    if there is a saved value to be used on the next call, it is
//    instead discarded.  This is useful if the user has reset the
//    random number seed, for instance.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8VEC_NORMAL_01_NEW[N], a sample of the standard normal PDF.
//
//  Local parameters:
//
//    Local, int MADE, records the number of values that have
//    been computed.  On input with negative N, this value overwrites
//    the return value of N, so the user can get an accounting of
//    how much work has been done.
//
//    Local, real R(N+1), is used to store some uniform random values.
//    Its dimension is N+1, but really it is only needed to be the
//    smallest even number greater than or equal to N.
//
//    Local, int SAVED, is 0 or 1 depending on whether there is a
//    single saved value left over from the previous call.
//
//    Local, int X_LO, X_HI, records the range of entries of
//    X that we need to compute.  This starts off as 1:N, but is adjusted
//    if we have a saved value that can be immediately stored in X(1),
//    and so on.
//
//    Local, real Y, the value saved from the previous call, if
//    SAVED is 1.
//
{
# define PI 3.141592653589793

  int i;
  int m;
  static int made = 0;
  double *r;
  static int saved = 0;
  double *x;
  int x_hi;
  int x_lo;
  static double y = 0.0;

  x = new double[n];
//
//  I'd like to allow the user to reset the internal data.
//  But this won't work properly if we have a saved value Y.
//  I'm making a crock option that allows the user to signal
//  explicitly that any internal memory should be flushed,
//  by passing in a negative value for N.
//
  if ( n < 0 )
  {
    made = 0;
    saved = 0;
    y = 0.0;
    return NULL;
  }
  else if ( n == 0 )
  {
    return NULL;
  }
//
//  Record the range of X we need to fill in.
//
  x_lo = 1;
  x_hi = n;
//
//  Use up the old value, if we have it.
//
  if ( saved == 1 )
  {
    x[0] = y;
    saved = 0;
    x_lo = 2;
  }
//
//  Maybe we don't need any more values.
//
  if ( x_hi - x_lo + 1 == 0 )
  {
  }
//
//  If we need just one new value, do that here to avoid null arrays.
//
  else if ( x_hi - x_lo + 1 == 1 )
  {
    r = uniform_01 ( 2, seed );

    x[x_hi-1] = sqrt ( - 2.0 * log ( r[0] ) ) * cos ( 2.0 * PI * r[1] );
    y =         sqrt ( - 2.0 * log ( r[0] ) ) * sin ( 2.0 * PI * r[1] );

    saved = 1;

    made = made + 2;

    delete [] r;
  }
//
//  If we require an even number of values, that's easy.
//
  else if ( ( x_hi - x_lo + 1 ) % 2 == 0 )
  {
    m = ( x_hi - x_lo + 1 ) / 2;

    r = uniform_01 ( 2*m, seed );

    for ( i = 0; i <= 2*m-2; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * PI * r[i+1] );
      x[x_lo+i  ] = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * PI * r[i+1] );
    }
    made = made + x_hi - x_lo + 1;

    delete [] r;
  }
//
//  If we require an odd number of values, we generate an even number,
//  and handle the last pair specially, storing one in X(N), and
//  saving the other for later.
//
  else
  {
    x_hi = x_hi - 1;

    m = ( x_hi - x_lo + 1 ) / 2 + 1;

    r = uniform_01 ( 2*m, seed );

    for ( i = 0; i <= 2*m-4; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * PI * r[i+1] );
      x[x_lo+i  ] = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * PI * r[i+1] );
    }

    i = 2*m - 2;

    x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * PI * r[i+1] );
    y           = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * PI * r[i+1] );

    saved = 1;

    made = made + x_hi - x_lo + 2;

    delete [] r;
  }

  return x;
# undef PI
}
//****************************************************************************80

//****************************************************************************80


double *uniform_01 ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_01_NEW returns a new unit pseudorandom R8VEC.
//  Modified:
//
//    19 August 2004
//  Comments:
//     
//    The original function written by John Burkardt was a horrible
//    random generator. I have used more 'standard' way of generating
//    uniform number here. There is better random generator referenced
//    in the litterature using stdC++11 standard (std::random_device).
//    but this seems not performing well.
//
//  Author:
//     Othman Benomar
{
  double *r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_01_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new double[n];

  for (int i = 0; i < n; i++ )
  {

    r[i] = ((double)rand()/(double)RAND_MAX);
  }

  return r;
}
