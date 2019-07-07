module inverseErf;

import std.math: abs, exp, M_2_SQRTPI;
import std.mathspecial: erf;

struct inErfCacheStruct {
  double guess;
  double fx;
  double fpx;
  double fppx;
}

class inErfClass {

  private:

  static void derivatives(double x, ref double fpx, ref double fppx) {
    // computes d/dx erf(x), and d^2/dx^2 erf(x)
    // because fpx and fppx can be generated from one another,
    // I will use a single local variable to compute them both
    double value = -(x * x);
    value = exp(value);
    value *= M_2_SQRTPI; // M_2_SQRTPI = (2/sqrt(pi))
    // f'(x) is (2/sqrt(pi)) * e^(-x^2)
    // we can set it now
    fpx = value;
    // f''x is (-4/sqrt(pi)) * e^(-x^2) * x
    // which is a factor of -2x away from what we already have
    fppx = value * (-2) * x;
  }

  static double halley(double x, double fx, double fpx, double fppx) {
    // implements Halley's Method
    // https://en.wikipedia.org/wiki/Halley%27s_method
    // given a guess x, and the value of some funciton f(x)
    // and the first and second derivatives of the same function f
    // we can produce a better guess x2
    double numerator = 2 * fx * fpx;
    double denominator = (2 * fpx * fpx) - (fx * fppx);
    double x2 = x - (numerator / denominator);
    return x2;
  }

  inErfCacheStruct zero;
  inErfCacheStruct[5] positiveCache;
  inErfCacheStruct[5] negativeCache;
  double[5] positiveBoundaries;
  double[5] negativeBoundaries;

  static void computeCache(ref inErfCacheStruct cache, double x) {
    cache.guess = x;
    cache.fx = erf(x);
    derivatives(x,cache.fpx,cache.fppx);
  }

  immutable inErfCacheStruct * getInitialGuess(double x) {
    if (x > 0) {
      if (x < positiveBoundaries[0]) {
        return cast(inErfCacheStruct*) &(zero);
      }
      for (long i = 1; i < 5; i++) {
        if (x < positiveBoundaries[i]) {
          return cast(inErfCacheStruct*) &(positiveCache[i-1]);
        }
      }
      return cast(inErfCacheStruct*) &(positiveCache[4]);
    }
    else {
      if (x > negativeBoundaries[0]) {
        return cast(inErfCacheStruct*) &(zero);
      }
      for (long i = 1; i < 5; i++) {
        if (x > negativeBoundaries[i]) {
          return cast(inErfCacheStruct*) &(negativeCache[i-1]);
        }
      }
      return cast(inErfCacheStruct*) &(negativeCache[4]);
    }
  }

  public:

  this() {
    // first do the zero stuff
    computeCache(zero,0);
    // then the other cache objects
    double t;
    for (long i = 0; i < 5; i++) {
      t = i + 1;
      computeCache( positiveCache[i] , t );
      computeCache( negativeCache[i] , -t );
    }
    // then the boundaries
    for (long i = 0; i < 5; i++) {
      t = i;
      t += 0.5;
      positiveBoundaries[i] = erf( t );
      negativeBoundaries[i] = erf( -t );
    }
  }

}

immutable inErfClass inErfObject;

static this() {
  inErfObject = cast(immutable inErfClass) new inErfClass();
}

double inErf(double x, double precision = 0.001) {
  // sanity check
  if (x < -1 || x > 1) {
    throw new Exception("InverseErf input must be in the open interval -1 to 1.");
  }
  if (precision <= 0) {
    throw new Exception("InverseErf precision must be greater than zero.");
  }
  // use cache to quickly find a good first guess
  inErfCacheStruct * initialGuess = inErfObject.getInitialGuess(x);
  // copy our guess to the local scope
  double guess = initialGuess.guess;
  double fx = initialGuess.fx - x;
  double fpx = initialGuess.fpx;
  double fppx = initialGuess.fppx;
  // now the loop
  double new_guess;
  while (true) {
    new_guess = inErfObject.halley(guess,fx,fpx,fppx);
    if (abs(new_guess - guess) < precision) {
      return new_guess;
    }
    else {
      guess = new_guess;
    }
    // we have updated our guess
    fx = erf(guess) - x;
    inErfObject.derivatives(guess,fpx,fppx);
    // we have updated our derivatives
  }
}
