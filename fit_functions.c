//
// The fit functions
// (to make a 3rd order polynomial)
//

#include "fit_functions.h"

double monomial_f0(double *x)
{
    return 1.0;
}

double monomial_f1(double *x)
{
    return x[0];
}

double monomial_f2(double *x)
{
    return x[0] * x[0];
}

double monomial_f3(double *x)
{
    return x[0] * x[0] * x[0];
}

double monomial_f4(double *x)
{
    double t = x[0] * x[0];
    return t * t;;
}


double (*polynomial_o4[5])(double *x) = {monomial_f0, monomial_f1, monomial_f2, monomial_f3, monomial_f4};
