# ellipse-perimeter
Calculate the perimeter of an ellipse down to the accuracy of the type of floating point used.

     cc -o ellipse-perimeter -Wall ellipse-perimeter.c -lm

This C code calculates the 'exact' perimeter of an ellipse, at least to the
precision of the machine word, by evaluating the sum of an infinite series, and
stopping once a sum term is smaller than the lowest bit in the floating point
representation.  Surprisingly very few terms are needed (7 for "long double")
and the calculation is quite efficient.  The floating point type can be changed
by the programmer to "float" or "double" if desired.

I started writing this after reading somewhere that there was no exact solution to
the elliptical perimeter problem.  In hindsight I realise now they meant no exact
solution that does not involve an infinite series, or that if you are using an
infinite series to calculate the perimeter, it is invalid because no-one can
ever fully evaluate an infinite series.  However this is academic because our
computers cannot calculate with infinite precision and a converging infinite series
*can* be evaluated down to the last bit of precision supported by our hardware.
 
So ... my approach to calculating the perimeter of an ellipse was this:

We know the perimeter of a circle (2*&pi;*r) and we know that any ellipse can
be formed by applying a shear on one axis to a circle and then rotating
the axes.  The rotation is trivial and we can make our calcuation simpler
by assuming an axis-aligned ellipse.
  
  When a circle of radius 'r' is stretched uniformly in one dimension by a
factor of k, it transforms into an ellipse with semi-major axis a = kr and
semi-minor axis b = r.
  
While searching for an equation to calculate the effect of the shear, I found
the web page:   https://www.mathsisfun.com/geometry/ellipse-perimeter.html
(which gave a formula for the perimeter of an ellipse directly and doesn't
mention the shear transformation) so I ended up implementing that instead!
(using the infinite series #2 from that page)
  
To compute the perimeter using the information from www.mathsisfun.com, first we
calculate h:

            (a-b)^2
        h = -------
            (a+b)^2
   
Then we use the sum of the infinite series:
   
        p = Pi*(a+b) * Sigma(n={0,Inf}, Binom(0.5,n)^2*h^n)
   
which expands to:
   
                             1               1               1              25              49                441
         p = Pi*(a+b) * (   --- * h^0   +   --- * h^1   +   --- * h^2  +   --- * h^4   +   ----- * h^5   +   ----- * h^6   +   ...   )
                             1               4               64            256             16384             65536 

We will add terms until the terms are so small that they no longer affect the cumulative sum.
This is sure to happen as the series is known to converge monotonically.

To calculate this in a C program we need to generate the numerators and denominators of this series.
Because of the size of the numbers generated in these sequences we must use floating point data types
rather than integers.
   
The denominators are relatively easy: the numbers are all of the form 4^n and taking log_4 of these numbers
reveals the n's to be: 0, 1, 3, 4, 7, 8, 10, 11, 15,

This series (which we can find at: https://oeis.org/A005187 ) can be formed by a recurrence relation,

     A005187(n) = n + A005187(floor(n/2));
   
The numerators are a little more tricky...

We find the sequence at https://oeis.org/A056981 where it is defined in terms of another series as:
   
     A056981(n) = A002596(n)^2.
   
The other series is found at https://oeis.org/A002596 where it is defined as:
   
     A002596(n+2) = C(n+1)/2^k(n+1)
   
when n >= 0.

     C(n) = A000108(n), k(n) = A048881(n).
     
('C()' are the Catalan numbers)
(  i.e. expressed alternatively, we have:

     A002596(n+2) = A000108(n+1)/2^A048881(n+1)
)
   
These are defined in turn by:
   
     A000108(n) = (2*n)!/(n!*(n+1)!);

However we will use the gamma function here instead of factorial, because our data type is floating point:

     A000108(n) = f(2*n)/(f(n)*f(n+1));    f(n) is available in C as 'tgamma(n+1)'.
   
and
   
     A048881(n) = A000120(n+1) - 1
   
Once again we need another series:
   
A000120: 1's-counting sequence: number of 1's in binary expansion of n (or the binary weight of n).
Since our data is floating point we can't use bit operations to determine this so our function 'wt'
uses division by two and remainder modulo two, to perform an analagous computation.
   
Putting these all together does indeed reproduce the series 1, 1, 1, 1, 25, 49, 441, 1089, 184041, etc.

Having worked all that out, I also found an approximate formula for the perimeter of the ellipse and
used the exact code to check it.  The results of the approximation are very accurate and the approximation
can probably be used for all real-world applications.
