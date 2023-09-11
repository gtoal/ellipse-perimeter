// cc -o ellipse-perimeter -Wall ellipse-perimeter.c -lm

// The code below calculates the 'exact' perimeter of an ellipse, at least to the
// precision of the machine word, by evaluating the sum of an infinite series, and
// stopping once a sum term is smaller than the lowest bit in the floating point
// representation.  Surprisingly very few terms are needed (7 for "long double")
// and the calculation is quite efficient.  The floating point type can be changed
// by the programmer to "float" or "double" if desired.


//  I started writing this after reading somewhere that there was no exact solution to
//  the elliptical perimeter problem.  In hindsight I realise now they meant no exact
//  solution that does not involve an infinite series, or that if you are using an
//  infinite series to calculate the perimeter, it is invalid because no-one can
//  ever fully evaluate an infinite series.  However this is academic because our
//  computers cannot calculate with infinite precision and a converging infinite series
//  *can* be evaluated down to the last bit of precision supported by our hardware.
//  
//  So ... my approach to calculating the perimeter of an ellipse was this:
//
//    We know the perimeter of a circle (2*Pi*r) and we know that any ellipse can
//    be formed by applying a shear on one axis to a circle and then rotating
//    the axes.  The rotation is trivial and we can make our calcuation simpler
//    by assuming an axis-aligned ellipse.
//  
//    When a circle of radius 'r' is stretched uniformly in one dimension by a
//  factor of k, it transforms into an ellipse with semi-major axis a = kr and
//  semi-minor axis b = r.
//  
//  While searching for an equation to calculate the effect of the shear, I found
//  the web page:   https://www.mathsisfun.com/geometry/ellipse-perimeter.html
//  (which gave a formula for the perimeter of an ellipse directly and doesn't
//  mention the shear transformation) so I ended up implementing that instead!
//  (using the infinite series #2 from that page)
//  

//  To compute the perimeter using the information from www.mathsisfun.com, first we
//  calculate h:
//
//       (a-b)^2
//   h = -------
//       (a+b)^2
//   
//   Then we use the sum of the infinite series:
//   
//         p = Pi*(a+b) * Sigma(n={0,Inf}, Binom(0.5,n)^2*h^n)
//   
//   which expands to:
//   
//                             1               1               1              25              49                441
//         p = Pi*(a+b) * (   --- * h^0   +   --- * h^1   +   --- * h^2  +   --- * h^4   +   ----- * h^5   +   ----- * h^6   +   ...   )
//                             1               4               64            256             16384             65536 
//
//   We will add terms until the terms are so small that they no longer affect the cumulative sum.
//   This is sure to happen as the series is known to converge monotonically.
//
//   To calculate this in a C program we need to generate the numerators and denominators of this series.
//   Because of the size of the numbers generated in these sequences we must use floating point data types
//   rather than integers.
//   
//   The denominators are relatively easy: the numbers are all of the form 4^n and taking log_4 of these numbers
//   reveals the n's to be: 0, 1, 3, 4, 7, 8, 10, 11, 15,
//
//   This series (which we can find at: https://oeis.org/A005187 ) can be formed by a recurrence relation,
//    A005187(n) = n + A005187(floor(n/2));
//   
//   The numerators are a little more tricky...
//
//   We find the sequence at https://oeis.org/A056981 where it is defined in terms of another series as:
//   
//    A056981(n) = A002596(n)^2.
//   
//   The other series is found at https://oeis.org/A002596 where it is defined as:
//   
//    A002596(n+2) = C(n+1)/2^k(n+1)
//   
//    when n >= 0.  C(n) = A000108(n), k(n) = A048881(n).   ('C()' are the Catalan numbers)
//    
//    (  i.e. expressed alternatively, we have: A002596(n+2) = A000108(n+1)/2^A048881(n+1)  )
//   
//    These are defined in turn by:
//   
//    A000108(n) = (2*n)!/(n!*(n+1)!);
//      However we will use the gamma function here instead of factorial, because our data type is floating point:
//      A000108(n) = f(2*n)/(f(n)*f(n+1));    f(n) is available in C as 'tgamma(n+1)'.
//   
//    and
//   
//    A048881(n) = A000120(n+1) - 1
//   
//    Once again we need another series:
//   
//    A000120: 1's-counting sequence: number of 1's in binary expansion of n (or the binary weight of n).
//    Since our data is floating point we can't use bit operations to determine this so our function 'wt'
//    uses division by two and remainder modulo two, to perform an analagous computation.
//   
//   Putting these all together does indeed reproduce the series 1, 1, 1, 1, 25, 49, 441, 1089, 184041, etc.
//
// Having worked all that out, I also found an approximate formula for the perimeter of the ellipse and
// used the exact code to check it.  The results of the approximation are very accurate and the approximation
// can probably be used for all real-world applications.

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#define LONG_DOUBLE 1
#define DOUBLE 2
#define FLOAT 3

#ifndef PRECISION  // can be overridden on command-line with -DPRECISION="FLOAT"
#define PRECISION LONG_DOUBLE     // Pick a precision!
#endif

#define NAME1(x) #x
#define NAME(x) NAME1(x)

#if PRECISION == LONG_DOUBLE
#undef PRECISION
#define PRECISION long double
#define FORMAT "%.16Lf"
#define POW(x,y) powl(x,y)
#define FABS(x)  fabsl(x)
#define SQRT(x)  sqrtl(x)
#define TGAMMA(x) tgammal(x)
#elif PRECISION == DOUBLE
#undef PRECISION
#define PRECISION double
#define FORMAT "%.16lf"
#define POW(x,y) pow(x,y)
#define FABS(x)  fabs(x)
#define SQRT(x)  sqrt(x)
#define TGAMMA(x)  tgamma(x)
#elif PRECISION == FLOAT
#undef PRECISION
#define PRECISION float
#define FORMAT "%.16f"
#define POW(x,y) powf(x,y)
#define FABS(x)  fabsf(x)
#define SQRT(x)  sqrtf(x)
#define TGAMMA(x)  tgammaf(x)
#endif

#define pow(x,y) ((y)<=0.0 ? 1.0 : POW(x,y)) // correct a problem with pow. X^0 should always be 1.
#define fabs(x) FABS(x)
#define sqrt(x) SQRT(x)
#define tgamma(x) TGAMMA(x)


static PRECISION Den(PRECISION n) {
  // calculate the denominator of the coefficient of each term
  // see https://oeis.org/A005187
  if (n == 0.0) return 0.0;
  return Den(floor(n/2)) + n;
}

static PRECISION f(PRECISION x) {
  PRECISION gamma = tgamma(fabs(x) + 1.0);

  if (isnan(gamma) || isinf(gamma)) {
    fprintf(stderr, "f(" FORMAT ") = NAN -> 1.0\n", x);
    return 1.0;
  }

  // fprintf(stderr, "f(" FORMAT ") = " FORMAT "\n", x, gamma);
  return gamma;
}

static PRECISION C(PRECISION n) {
  PRECISION r = floor(f(2.0*n))/(floor(f(n))*floor(f(n+1.0)));
  //fprintf(stderr, "C(" FORMAT ") -> " FORMAT "\n", n, r);
  return round(r);
}

static PRECISION wt(PRECISION n) {
  // wt(n) = 1's-counting sequence: number of 1's in binary expansion of n
  PRECISION tot = 0.0;
  if (n < 0.0) n = -n;
  while (n > 0.25) {
    if (round(n - floor(n/2.0)*2.0) == 1.0) tot++;
    n = floor(n/2.0);
  }
  //fprintf(stderr, "wt(" FORMAT ") -> " FORMAT "\n", n, tot);
  return tot;
}

static PRECISION k(PRECISION n) {
  if (n <= 0.0) return 1.0;
  PRECISION r1 = wt(n+1.0);
  PRECISION r = r1 - 1.0;
  //fprintf(stderr, "k(" FORMAT ") -> wt(" FORMAT ")-1.0 -> " FORMAT " - 1.0 -> " FORMAT "\n", n, n+1.0, r1, r);
  return r;
}

static PRECISION a(PRECISION n) {
  PRECISION r, t1, t2, tk;
  n = n-2.0;
  t1=C(n+1.0);
  tk=k(n+1.0);
  t2=pow(2.0, tk);  // n >= 0; where C(n) = A000108(n), k(n) = A048881(n)
  r = round(t1/t2);
  //fprintf(stderr, "a(" FORMAT ") -> C(" FORMAT ")/pow(2, k(" FORMAT ")) -> " FORMAT "/pow(2, " FORMAT ") -> " FORMAT "/" FORMAT " -> " FORMAT "^2 -> " FORMAT "\n",
  //                       n,              n+1,                 n+1,             t1,               tk,            t1,        t2,            r,            r*r);
  return r*r;
}

static PRECISION Num(int n) {
  // calculate the numerator of the coefficient of each term. see https://oeis.org/A056981
  PRECISION r;
  r = a(n);
  //fprintf(stderr, "COMPARE: " FORMAT " vs " FORMAT "\n", lookup[n], a(n));
  return r;
}

// Axis-aligned.  If not, rotate it first until it *is* axis-aligned.
// If the input ellipse is the result of a shear, use Nate's formula
// to convert to a rotation first.  Need to experiment to determine
// if a>b or b>a leads to a more accurate or quicker-to-compute result.

PRECISION exact_ellipse_perimeter(PRECISION a, PRECISION b) {
  int term;

  if (a == 0.0) return b*4.0; // Special-case for ellipse flattened to a line!
  if (b == 0.0) return a*4.0;

  PRECISION h = ((a - b) * (a - b)) / ((a + b) * (a + b));
  PRECISION p = M_PI * (a + b);
  PRECISION result, numerator, denominator, coefficient, infinite_sum = 0.0, prev_result = 0.0;

  // for (term = 0; term < sizeof(lookup)/sizeof(lookup[0]); term++) {
  for (term = 0; ; term++) {
    numerator = Num(term);
    denominator = pow(4.0, Den(term));
    coefficient = numerator/denominator;
    infinite_sum += coefficient * pow(h, term);
    result = p * infinite_sum;
    if (result == prev_result) {
      fprintf(stderr, "\nWe hit the limit of machine accuracy (using type \"" NAME(PRECISION) "\") after %d terms in the expansion of the infinite series\n", term-1);
      break;
    }
    prev_result = result;
  }
  return result;
}

PRECISION approximate_ellipse_perimeter(PRECISION a, PRECISION b) {
  // This is pretty accurate and can probably be used in every real-world application. 
  PRECISION h = pow((a-b)/(a+b), 2.0);
  return M_PI*(a + b)*(1 + (3*h)/(10 + sqrt(4 - 3*h)));
}

int main(int argc, char **argv) {
  // Use the 'exact' calculation to verify the approximate calculation above.
  PRECISION a = 1.2, b = 1.0; // later: get from argv instead.
  
  PRECISION approximate_p = approximate_ellipse_perimeter(a,b);
  PRECISION exact_p       = exact_ellipse_perimeter(a,b);

  PRECISION err = fabs((exact_p - approximate_p)/exact_p);

  fprintf(stderr, "\nThe perimeter of an ellipse whose major axis is rx=" FORMAT " and ry=" FORMAT " is " FORMAT "\n", a, b, exact_p);
  fprintf(stderr, "\nA good approximation to the perimeter is " FORMAT " with error=" FORMAT "%%\n", approximate_p, err*100.0);

  fprintf(stderr, "\nThe cruder approximation of 2 * PI * sqrt( (rx + ry) / 2 ) is " FORMAT " and should never be used!\n\n", 2 * M_PI * sqrt( (a + b) / 2 ));
  exit(EXIT_SUCCESS);
  return 0;
}
