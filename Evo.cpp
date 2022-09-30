#if HAVE_CONFIG_H
# include <config.h>
#endif

#include <cfloat>
#include <cmath>

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_min.h>

#include "Evo.h"

#include <math.h>       /* pow */
#include <iostream>
#include <stdlib.h>
#include <functional>



/*cosmological parameters for evolution factors (applied to analytic solutions) */
struct param_type2 {
	double omega00;
	double omega_cb;
	double hubble;
  double par1;
	double par2;
	double par3;
};

int jac (double a, const double G[], double *dfdy, double dfdt[], void *params)
{
	return GSL_SUCCESS;
}

// fixed omega
// H/H0
static double HAde(double a, double omega0, double w0, double wa){
  double A = -3.*(1.+w0+wa);
  double omegaf = pow(a,A)*exp(3.*(-1.+a)*wa);
	double omegaL= (1.-omega0)*omegaf;
	return  sqrt(omega0/pow(a,3)+omegaL);
}
// aH dH/da / H0^2
static double HA1de(double a, double omega0, double w0, double wa){
  double A = -3.*(1.+w0+wa);
  double omegaf = pow(a,A)*exp(3.*(-1.+a)*wa);
  double omegaL= (1.-omega0)*omegaf;
	return -3*omega0/(2.*pow(a,3)) + (A+3.*a*wa)*omegaL/2.;
}
//HA2de = -dH/dt/H^2 = -a dH/da / H
static double HA2de(double a, double omega0, double w0, double wa){
  return -HA1de(a,omega0,w0,wa)/pow(HAde(a,omega0,w0,wa),2);
}
//3/(2H^2) * Omega_m_0 *a^3 = Omega_m(a)
static double HA2de2(double a, double omega0, double w0, double wa){
  return 3.*omega0/(2.*pow(HAde(a,omega0,w0,wa),2)*pow(a,3));
}


// linear growth factor for wCDM CPL with interaction
int evo_eqs_ide(double a, const double G[], double F[], void *params)
{
	param_type2 p = *(param_type2 *)(params);
	double omega0 = p.omega00;
	double omega_cb = p.omega_cb;
  double hubble=p.hubble;
	double w0 = p.par1;
	double wa = p.par2;
	double xi = p.par3;


	double A, omegaf, omegaL, fric,hade1,hade2,f_cb;

	A = -3.*(1.+w0+wa);
	omegaf = pow(a,A)*exp(3.*(-1.+a)*wa);
	omegaL= (1.-omega0)*omegaf;
	fric = (1.+w0+(1.-a)*wa)*omegaL*xi*hubble/HAde(a,omega0,w0,wa)*0.0974655;
	hade1 = HA2de(a,omega0,w0,wa);
	hade2 = HA2de2(a,omega0,w0,wa);

	f_cb = omega_cb/omega0;

	F[0] = -G[1]/a;
	F[1] =1./a*(-(2.+fric-hade1)*G[1]-f_cb*hade2*G[0]);

	return GSL_SUCCESS;
}


void get_growth(double A, double omega0, double omegacb, double h,double w0, double wa, double xi, double * res, double accuracy)
{

	struct param_type2 params_my = {omega0,omegacb,h,w0,wa,xi};

	//Initial a
	double a = 0.0001;

	double G[2] = { a,-a};

	gsl_odeiv2_system sys = {evo_eqs_ide, jac, 2, &params_my};

	//  this gives the system, the method, the initial interval step, the absolute error and the relative error.
	gsl_odeiv2_driver * d =
	gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
								 1e-3, accuracy, 0.);
	int status1 = gsl_odeiv2_driver_apply (d, &a, A, G);


	res[0] = G[0];
	res[1] = -G[1]/G[0];

	gsl_odeiv2_driver_free(d);


}
