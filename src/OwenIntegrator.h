#ifndef OWENINTEGRATOR_H
#define OWENINTEGRATOR_H

#include "DEIntegrationConstants.h"
#include <float.h>

/*! Numerical integration in one dimension using the double expontial method of M. Mori. */
class OwenIntegrator
{
public:
    static double Integrate
    (
        double H,                       //!< [in] h parameter
        double b,                       //!< [in] right limit of integration (a)
        double targetAbsoluteError,     //!< [in] desired bound on error
        int& numFunctionEvaluations,    //!< [out] number of function evaluations used
        double& errorEstimate           //!< [out] estimated error in integration
    )
    {
        // Apply the linear change of variables x = ct + d
        // $$\int_a^b f(x) dx = c \int_{-1}^1 f( ct + d ) dt$$
        // c = (b-a)/2, d = (a+b)/2
        // b <- a; a <- 0

        double d = 0.5*b;

        return IntegrateCore
        (
            H,
            d,
            targetAbsoluteError,
            numFunctionEvaluations,
            errorEstimate,
            doubleExponentialAbcissas,
            doubleExponentialWeights
        );
    }

    /*! This version overloaded to not require arguments passed in for
        function evaluation counts or error estimates.
        @return The value of the integral.
    */
    static double Integrate
    (
        double H,                       //!< [in] h parameter
        double b,                       //!< [in] right limit of integration (a)
        double targetAbsoluteError      //!< [in] desired bound on error
    )
    {
        int numFunctionEvaluations;
        double errorEstimate;
        return Integrate
        (
            H,
            b,
            targetAbsoluteError,
            numFunctionEvaluations,
            errorEstimate
        );
    }

private:
    static double f(double h, double x)
    {
      return exp (-(h*h)*(1+x*x)/2) / (1+x*x);
    }

    // Integrate f(cx + d) with the given integration constants
    static double IntegrateCore
    (
        double H,
        double d,   // intercept of change of variables
        double targetAbsoluteError,
        int& numFunctionEvaluations,
        double& errorEstimate,
        const double* abcissas,
        const double* weights
    )
    {
        targetAbsoluteError /= d;

        // Offsets to where each level's integration constants start.
        // The last element is not a beginning but an end.
        int offsets[] = {1, 4, 7, 13, 25, 49, 97, 193};
        int numLevels = sizeof(offsets)/sizeof(int) - 1;

        //double newContribution = 0.0;
        //double integral = 0.0;
        errorEstimate = DBL_MAX;

        double integral = f(H, d*abcissas[0] + d) * weights[0];
        int i;
        for (i = offsets[0]; i != offsets[1]; ++i)
            integral += weights[i]*(f(H, d*abcissas[i] + d) + f(H, -d*abcissas[i] + d));

        {
          double h = 1.0;
          double previousDelta, currentDelta = DBL_MAX;
          for (int level = 1; level != numLevels; ++level)
          {
              h *= 0.5;
              double newContribution = 0.0;
              for (i = offsets[level]; i != offsets[level+1]; ++i)
                  newContribution += weights[i]*(f(H, d*abcissas[i] + d) + f(H, -d*abcissas[i] + d));
              newContribution *= h;

              // difference in consecutive integral estimates
              previousDelta = currentDelta;
              currentDelta = fabs(0.5*integral - newContribution);
              integral = 0.5*integral + newContribution;

              // Once convergence kicks in, error is approximately squared at each step.
              // Determine whether we're in the convergent region by looking at the trend in the error.
              if (level == 1)
                  continue; // previousDelta meaningless, so cannot check convergence.

              // Exact comparison with zero is harmless here.  Could possibly be replaced with
              // a small positive upper limit on the size of currentDelta, but determining
              // that upper limit would be difficult.  At worse, the loop is executed more
              // times than necessary.  But no infinite loop can result since there is
              // an upper bound on the loop variable.
              if (currentDelta == 0.0)
                  break;
              double r = log( currentDelta )/log( previousDelta );  // previousDelta != 0 or would have been kicked out previously

              if (r > 1.9 && r < 2.1)
              {
                  // If convergence theory applied perfectly, r would be 2 in the convergence region.
                  // r close to 2 is good enough. We expect the difference between this integral estimate
                  // and the next one to be roughly delta^2.
                  errorEstimate = currentDelta*currentDelta;
              }
              else
              {
                  // Not in the convergence region.  Assume only that error is decreasing.
                  errorEstimate = currentDelta;
              }

              if (errorEstimate < 0.1*targetAbsoluteError)
                  break;
          }
        }
        numFunctionEvaluations = 2*i - 1;
        errorEstimate *= d;
        return d*integral;
    }

};

#endif // include guard
