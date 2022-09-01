#include <stdlib.h>
#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

double rad = 1.0; // radius
size_t ndim = 3; // dimension
double exact = pow(M_PI,ndim/2.0) / gsl_sf_gamma(ndim/2.0+1.0) * pow(rad,ndim);

double integrand (double *k, size_t dim, void *params)
{
  (void)(dim);  (void)(params);
  double r = 0.0;
  for(size_t i=0; i<dim; i++) r += k[i]*k[i];
  if(r <= rad*rad) return 1.0;
  else return 0.0;  
}

int main (void)
{
  // some dummy variables to capture the result
  double res, sig;

  // define upper and lower limit
  double xl[ndim], xu[ndim];
  std::fill_n (xl, ndim, -rad);
  std::fill_n (xu, ndim, +rad);

  // set up integrand function
  gsl_monte_function F = { &integrand, ndim, 0 };

  // set up random number generator
  const gsl_rng_type *T;
  gsl_rng *r;
  gsl_rng_env_setup ();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  // monte-carlo points
  size_t calls = 500000;

  std::cout << "exact result:" << std::endl;
  std::cout << "result: " << exact << "\n\n";

  // plain integration
  std::cout << "plain integration:" << std::endl;
  gsl_monte_plain_state *s_p = gsl_monte_plain_alloc (ndim);
  gsl_monte_plain_integrate (&F, xl, xu, ndim, calls, r, s_p, &res, &sig);
  gsl_monte_plain_free (s_p);
  std::cout << "result: " << res 
    << "\t\t" << "error: " << exact - res 
    << "\t\t" << "sigma: " << sig << "\n\n";
  
  // miser integration
  std::cout << "miser integration:" << std::endl;
  gsl_monte_miser_state *s_m = gsl_monte_miser_alloc (ndim);
  gsl_monte_miser_integrate (&F, xl, xu, ndim, calls, r, s_m, &res, &sig);
  gsl_monte_miser_free (s_m);
  std::cout << "result: " << res 
    << "\t\t" << "error: " << exact - res 
    << "\t\t" << "sigma: " << sig << "\n\n";
  
  // vegas warm-up integration
  std::cout << "vegas warm-up integration:" << std::endl;
  gsl_monte_vegas_state *s_v = gsl_monte_vegas_alloc (ndim);
  gsl_monte_vegas_integrate (&F, xl, xu, ndim, 10000, r, s_v, &res, &sig);
  std::cout << "result: " << res 
    << "\t\t" << "error: " << exact - res 
    << "\t\t" << "sigma: " << sig << "\n\n";

  // vegas final integration
  std::cout << "vegas final integration:" << std::endl;
  size_t it = 0;
  do{
    gsl_monte_vegas_integrate (&F, xl, xu, ndim, 10000, r, s_v, &res, &sig);
    it++;
    std::cout << "iteration = " << it 
      << ",\t|chi^2-1| = " << fabs(gsl_monte_vegas_chisq(s_v)-1.0) 
      << ",\tsd = " << sig << std::endl;
  }while (fabs (gsl_monte_vegas_chisq (s_v) - 1.0) > 0.1);
  gsl_monte_vegas_free (s_v);
  std::cout << "result: " << res 
    << "\t\t" << "error: " << exact - res 
    << "\t\t" << "sigma: " << sig << "\n\n";

  // finalize
  gsl_rng_free (r);
  return 0;
}
