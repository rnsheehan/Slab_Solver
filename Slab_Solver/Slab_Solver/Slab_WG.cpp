#ifndef ATTACH_H
#include "Attach.h"
#endif

// Definitions for the slab base class

slab::slab()
{
	// Default Constructor
	M = 0;

	k_sqr = nc_sqr = ncl_sqr = ns_sqr = nr_sqr = k_sqr_nc_sqr = k_sqr_ns_sqr = k_sqr_ncl_sqr = k_sqr_nr_sqr = 0.0;
	aa = bb = d = l = nc = ns = ncl = nr = 0.0;
	V = na = k = lower = upper = upper = w = efieldint = hfieldint = 0.0;
}

slab::~slab()
{
	// Deconstructor
	betaE.clear();

	betaH.clear(); 
}

int slab::nbeta(bool mode)
{
	// return the number of computed propagation constants for a given polarisation
	// it can be the case that the predicted value M is greater than the actual number of modes in a waveguide

	return static_cast<int>(mode ? betaE.size() : betaH.size() );
}

double slab::h(int i, bool t)
{
	//Wavenumber in Core

	try {

		if (nbeta(t) > 0) {
			if (i<0 || i > nbeta(t)) {
				throw std::range_error("Error: double slab_tl::h(int i, bool t)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {
				if (t) {
					return sqrt(k_sqr_nc_sqr - template_funcs::DSQR(betaE[i]));
				}
				else {
					return sqrt(k_sqr_nc_sqr - template_funcs::DSQR(betaH[i]));
				}
			}
		}
		else {
			throw std::invalid_argument("Error: double slab::h(int i, bool t)\nNo modes have been computed\n");
			return 0; 
		}
	}
	catch (std::range_error &e) {
		std::cerr << e.what();
		return 0;
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab::p(int i, bool t)
{
	//Wavenumber in Substrate

	try {
		if (nbeta(t) > 0) {
			if (i<0 || i > nbeta(t)) {
				throw std::range_error("Error: double slab_tl::p(int i, bool t)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {
				if (t) {
					return sqrt(template_funcs::DSQR(betaE[i]) - k_sqr_ns_sqr);
				}
				else {
					return sqrt(template_funcs::DSQR(betaH[i]) - k_sqr_ns_sqr);
				}
			}
		}
		else {
			throw std::invalid_argument("Error: double slab::p(int i, bool t)\nNo modes have been computed\n");
			return 0;
		}
	}
	catch (std::range_error &e) {
		std::cerr << e.what();
		return 0;
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab::q(int i, bool t)
{
	//Wavenumber in Cladding

	try {

		if (nbeta(t) > 0) {

			if (i<0 || i > nbeta(t)) {
				throw std::range_error("Error: double slab_tl::q(int i, bool t)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {
				if (t) {
					return sqrt(template_funcs::DSQR(betaE[i]) - k_sqr_ncl_sqr);
				}
				else {
					return sqrt(template_funcs::DSQR(betaH[i]) - k_sqr_ncl_sqr);
				}
			}

		}
		else {
			throw std::invalid_argument("Error: double slab::q(int i, bool t)\nNo modes have been computed\n");
			return 0;
		}
	}
	catch (std::range_error &e) {
		std::cerr << e.what();
		return 0;
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab::r(int i, bool t)
{
	//Wavenumber in Rib

	try {

		if (nbeta(t) > 0) {

			if (i<0 || i > nbeta(t)) {
				throw std::range_error("Error: double slab_tl::r(int i, bool t)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {
				if (t) {
					return sqrt( template_funcs::DSQR(betaE[i]) - k_sqr_nr_sqr);
				}
				else {
					return sqrt( template_funcs::DSQR(betaH[i]) - k_sqr_nr_sqr);
				}
			}

		}
		else {
			throw std::invalid_argument("Error: double slab::r(int i, bool t)\nNo modes have been computed\n");
			return 0;
		}
	}
	catch (std::range_error &e) {
		std::cerr << e.what();
		return 0;
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

// Definitions for the three layer slab derived class
slab_tl::slab_tl(double width, double lambda, double ncore, double nsub, double nclad)
{
	set_params(width, lambda, ncore, nsub, nclad);
}

void slab_tl::set_params(double width, double lambda, double ncore, double nsub, double nclad)
{
	// assign values to the parameters for the three layer slab waveguide
	// throw exception if one of the inputs is not correct

	try {

		bool c1 = width > 0.0 ? true : false;
		bool c2 = lambda > 0.0 ? true : false;
		bool c3 = nclad >= 1.0 ? true : false;
		bool c4 = nsub >= 1.0 ? true : false;
		bool c5 = ncore > std::max(nsub, nclad) ? true : false;

		if (c1 && c2 && c3 && c4 && c5) {

			d = width;

			l = lambda;

			nc = ncore;
			nc_sqr = template_funcs::DSQR(nc);

			ns = std::max(nsub, nclad);
			ns_sqr = template_funcs::DSQR(ns);

			ncl = std::min(nsub, nclad);
			ncl_sqr = template_funcs::DSQR(ncl);

			aa = (nc_sqr / ns_sqr);

			bb = (nc_sqr / ncl_sqr);

			k = Two_PI / l;
			k_sqr = template_funcs::DSQR(k); // k_{0}^{2}

			k_sqr_nc_sqr = k_sqr * nc_sqr; // k_{0}^{2} n_{c}^{2}
			k_sqr_ns_sqr = k_sqr * ns_sqr; // k_{0}^{2} n_{c}^{2}
			k_sqr_ncl_sqr = k_sqr * ncl_sqr; // k_{0}^{2} n_{c}^{2}

			double x = nc_sqr - ns_sqr;
			double y = ns_sqr - ncl_sqr;

			// WG asymmetry paramater
			// g = ( template_funcs::DSQR(ns) - template_funcs::DSQR(ncl) )/( template_funcs::DSQR(nc) - template_funcs::DSQR(ns) );
			g = y > 0.0 ? (y / x) : 0.0;

			// na = sqrt( template_funcs::DSQR(nc) - template_funcs::DSQR(ns) );
			na = sqrt(x); // numerical aperture

			//V = ( PI*d*sqrt( template_funcs::DSQR(nc) - template_funcs::DSQR(ns) ) ) / l;
			V = (PI*d*na) / l; // V-parameter

			// predicted number of modes
			M = y > 0 ? static_cast<int>( std::max(1.0, ceil((2.0*V / PI) - atan(g) / PI)) ) : static_cast<int>(std::max(1.0, ceil((2.0*V / PI))));

			lower = k * ns; // lower bound of search space

			upper = k * nc; // upper bound of search space

			w = k * SPEED_OF_LIGHT;

			efieldint = 0.5*EPSILON*SPEED_OF_LIGHT;

			hfieldint = 0.5*MU*SPEED_OF_LIGHT;

			// Empty the std::vector each time a new instance of the class is called
			betaE.clear();
			betaH.clear();
		}
		else {
			std::string reason = "Error: void slab_tl::set_params(double width,double lambda,double ncore,double nsub,double nclad)\n";
			if (!c1) reason += "WG width = " + template_funcs::toString(width, 3) + " is negative\n";
			if (!c2) reason += "Wavelength = " + template_funcs::toString(lambda, 3) + " is negative\n";
			if (!c3) reason += "Cladding Index = " + template_funcs::toString(nclad, 3) + " is less than one\n";
			if (!c4) reason += "Substrate Index = " + template_funcs::toString(nsub, 3) + " is less than one\n";
			if (!c5) reason += "Core Index = " + template_funcs::toString(ncore, 3) + " is less than cladding / substrate index\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab_tl::g1(int i, bool t)
{
	try {

		if (nbeta(t)) {

			if (i<0 || i > nbeta(t)) {
				throw std::range_error("Error: double slab_tl::g1(int i, bool t)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {
				if (t) {
					double beta_sqr = template_funcs::DSQR(betaE[i]);
					return (beta_sqr / k_sqr_nc_sqr) + (beta_sqr / k_sqr_ns_sqr) - 1.0;
				}
				else {
					double beta_sqr = template_funcs::DSQR(betaH[i]);
					return (beta_sqr / k_sqr_nc_sqr) + (beta_sqr / k_sqr_ns_sqr) - 1.0;
				}
			}
		}
		else {
			throw std::invalid_argument("Error: double slab_tl::g1(int i, bool t)\nNo modes have been computed\n");
			return 0;
		}
	}
	catch (std::range_error &e) {
		std::cerr << e.what();
		return 0;
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab_tl::g2(int i, bool t)
{
	try {

		if (nbeta(t) > 0) {

			if (i<0 || i > nbeta(t)) {
				throw std::range_error("Error: double slab_tl::g2(int i, bool t)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {
				if (t) {
					double beta_sqr = template_funcs::DSQR(betaE[i]);
					return (beta_sqr / k_sqr_nc_sqr) + (beta_sqr / k_sqr_ncl_sqr) - 1.0;
				}
				else {
					double beta_sqr = template_funcs::DSQR(betaH[i]);
					return (beta_sqr / k_sqr_nc_sqr) + (beta_sqr / k_sqr_ncl_sqr) - 1.0;
				}
			}
		}
		else {
			throw std::invalid_argument("Error: double slab_tl::g2(int i, bool t)\nNo modes have been computed\n");
			return 0;
		}
	}
	catch (std::range_error &e) {
		std::cerr << e.what();
		return 0;
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab_tl::deff(int i, bool t)
{
	// Effective width of waveguide mode

	try {

		if (nbeta(t) > 0) {
			if (i<0 || i > nbeta(t)) {
				throw std::range_error("Error: double slab_tl::deff(int i, bool t)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {
				if (t) {
					return d + (1.0 / p(i, t)) + (1.0 / q(i, t));
				}
				else {
					return d + (1.0 / (g1(i, t) * p(i, t))) + (1.0 / (g2(i, t) * q(i, t)));
				}
			}
		}
		else {
			throw std::invalid_argument("Error: double slab_tl::deff(int i, bool t)\nNo modes have been computed\n");
			return 0;
		}
	}
	catch (std::range_error &e) {
		std::cerr << e.what();
		return 0;
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab_tl::phase(int i, bool t)
{
	// Phase of mode in slab waveguide

	try {

		if (nbeta(t) > 0) {
			if (i<0 || i > nbeta(t)) {
				throw std::range_error("Error: double slab_tl::phase(int i, bool t)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {
				if (t) {
					double hh = h(i, t);
					double pp = p(i, t);
					double qq = q(i, t);
					return i * (PI / 2) + 0.5 * atan((pp / hh)) - 0.5 * atan((qq / hh));
				}
				else {
					double hh = h(i, t);
					double pp = p(i, t);
					double qq = q(i, t);
					return  i * (PI / 2) + 0.5*atan(aa*(pp / hh)) - 0.5*atan(bb * (qq / hh));
				}
			}
		}
		else {
			throw std::invalid_argument("Error: double slab_tl::phase(int i, bool t)\nNo modes have been computed\n");
			return 0;
		}
	}
	catch (std::range_error &e) {
		std::cerr << e.what();
		return 0;
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab_tl::norm_const(int i, bool t)
{
	// Normalisation constant for mode in a slab waveguide

	try {

		if (nbeta(t) > 0) {
			if (i<0 || i > nbeta(t)) {
				throw std::range_error("Error: double slab_tl::norm_const(int i, bool t)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {
				if (t) {
					return sqrt((w * MU) / (betaE[i] * deff(i, t)));
				}
				else {
					return sqrt((w * EPSILON * nc_sqr) / (betaH[i] * deff(i, t)));
				}
			}
		}
		else {
			throw std::invalid_argument("Error: double slab_tl::norm_const(int i, bool t)\nNo modes have been computed\n");
			return 0;
		}
	}
	catch (std::range_error &e) {
		std::cerr << e.what();
		return 0;
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab_tl::conf_fact(int i, bool t)
{
	// Confinement factor of mode in slab waveguide

	try {

		if (nbeta(t) > 0) {
			if (i<0 || i > nbeta(t)) {
				throw std::range_error("Error: double slab_tl::conf_fact(int i, bool t)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {
				if (t) {
					double hh_sqr = template_funcs::DSQR(h(i, t));
					double pp = p(i, t);
					double pp_sqr = template_funcs::DSQR(pp);
					double qq = q(i, t);
					double qq_sqr = template_funcs::DSQR(qq);

					return (d + (1.0 / pp)*(1.0 / (1.0 + hh_sqr / pp_sqr)) + (1.0 / qq)*(1.0 / (1.0 + hh_sqr / qq_sqr))) / (deff(i, t));
				}
				else {
					double hh_sqr = template_funcs::DSQR(h(i, t));
					double pp = p(i, t);
					double pp_sqr = template_funcs::DSQR(pp);
					double qq = q(i, t);
					double qq_sqr = template_funcs::DSQR(qq);
					return (d + (1.0 / (g1(i, t) * pp))*(1.0 / (1.0 + hh_sqr / pp_sqr)) + (1.0 / (g2(i, t) * qq))*(1.0 / (1.0 + hh_sqr / qq_sqr))) / (deff(i, t));
				}
			}
		}
		else {
			throw std::invalid_argument("Error: double slab_tl::conf_fact(int i, bool t)\nNo modes have been computed\n");
			return 0;
		}
	}
	catch (std::range_error &e) {
		std::cerr << e.what();
		return 0;
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab_tl::eigeneqn(double x, bool t)
{
	// Dispersion equation for slab waveguide written as a sum of sine and cosine functions
	// The version in eigeneqn_3, written as a sum over atan(*) is preferred for fast root finding

	try {

		if (x >= lower && x <= upper) {

			double x_sqr = template_funcs::DSQR(x);
			double h = sqrt(k_sqr_nc_sqr - x_sqr);
			double p = sqrt(x_sqr - k_sqr_ns_sqr);
			double q = sqrt(x_sqr - k_sqr_ncl_sqr);

			if (t) {//TE modes
				return (template_funcs::DSQR(h) - (p*q)) * sin((d * h)) - h * (p + q) * cos((d * h));
			}
			else {//TM modes
				return (template_funcs::DSQR(h) - aa * bb * (p * q)) * sin((d * h)) - (h * aa * p + h * bb * q) * cos((d * h));
			}
		}
		else {
			return 0;
			throw std::range_error("Error: double slab_tl::eigeneqn(double x,bool t)\n Attempting to compute dispersion equation outside valid range\n");
		}
	}
	catch (std::range_error &e) {
		std::cerr << e.what();
		return 0;
	}
}

double slab_tl::TE_TM(double x, int i, bool mode)
{
	// Function which defines the shape of the modes
	// This method uses the stored computed propagation constants
	// Assumes that the slab wg has core region defined on -t < x < t, t = d/2

	try {

		if (nbeta(mode) > 0) {

			if (i<0 || i>nbeta(mode)) {
				throw std::range_error("Error: double slab_tl::TE_TM(double x,int i,bool mode)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {
				double t = d / 2;

				if (x > t) {
					return norm_const(i, mode)*cos(h(i, mode)*t - phase(i, mode))*exp(-q(i, mode)*(x - t)); // Cladding
				}
				else if (x < -t) {
					return norm_const(i, mode)*cos(h(i, mode)*t + phase(i, mode))*exp(p(i, mode)*(x + t)); // Substrate
				}
				else {
					return norm_const(i, mode)*cos(h(i, mode)*x - phase(i, mode)); // Core
				}
			}
		}
		else {
			throw std::invalid_argument("Error: double slab_tl::TE_TM(int i, bool t)\nNo modes have been computed\n");
			return 0;
		}
	}
	catch (std::range_error &e) {
		std::cerr << e.what();
		return 0;
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab_tl::eigeneqn_3(double x, bool t, int mm)
{
	// Using this version of the dispersion equation means that you don't have to do a bracketing step
	// which you would have to do if you used eigeneqn

	try {

		if (k_sqr_nc_sqr > k_sqr_ns_sqr) {

			double h, p, q, tmp;

			double x_sqr = template_funcs::DSQR(x);

			tmp = k_sqr_nc_sqr - x_sqr;
			h = (tmp>0 ? sqrt(tmp) : 0.0);

			tmp = x_sqr - k_sqr_ns_sqr;
			p = (tmp>0 ? sqrt(tmp) : 0.0);

			tmp = x_sqr - k_sqr_ncl_sqr;
			q = (tmp>0 ? sqrt(tmp) : 0.0);

			if (t) {//TE modes
				return d * h - mm * PI - atan(q / h) - atan(p / h);
			}
			else {//TM modes
				return d * h - mm * PI - atan(bb*(q / h)) - atan(aa*(p / h));
			}
		}
		else {
			std::string reason = "Error: double slab_tl::eigeneqn_3(double x,bool t,int mm)\n";
			reason += "Input parameters not correctly defined\n";
			reason += "k_{0}^{2} n_{c}^{2} = " + template_funcs::toString(k_sqr_nc_sqr) + ", k_{0}^{2} n_{s}^{2} = " + template_funcs::toString(k_sqr_ns_sqr) + "\n";
			return 0;
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab_tl::zbrent(double x1, double x2, double tol, bool t, int mm)
{
	//Find the root of the function between x1 and x2 using brent's method
	//The root is refined to +/- tol
	//Seemingly one of the best methods to use

	//This will be used to compute the roots of eigeneq_3
	//R. Sheehan 28 - 5 - 2010

	try {

		bool c1 = fabs(x1 - x2) > 1.0e-9 ? true : false; // cannot have x1 == x2
		bool c2 = tol > 1.0e-16 ? true : false;
		bool c3 = mm < M ? true : false;

		if (c1 && c2 && c3) {

			int iter;

			static const int ITMAX = 100;//Used in rtbis, rtflsp, rtsec, zriddr

			double a = std::min(x1, x2), b = std::max(x1, x2), c = std::max(x1, x2), d, e, min1, min2;
			double fc, p, q, r, s, tol1, xm;
			double fa = eigeneqn_3(a, t, mm), fb = eigeneqn_3(b, t, mm);

			if ((fa>0.0 && fb>0.0) || (fa<0.0 && fb<0.0)) {
				std::cerr << "Root must be bracketed in zbrent\n";
			}
			fc = fb;
			for (iter = 1; iter <= ITMAX; iter++) {
				if ((fb>0.0 && fc>0.0) || (fb<0.0 && fc<0.0)) {
					c = a;
					fc = fa;
					e = d = b - a;
				}
				if (fabs(fc)<fabs(fb)) {
					a = b;
					b = c;
					c = a;
					fa = fb;
					fb = fc;
					fc = fa;
				}
				tol1 = 2.0*EPS*fabs(b) + 0.5*tol;
				xm = 0.5*(c - b);
				if (fabs(xm) <= tol1 || fb == 0.0) return b;
				/*if(fabs(xm)<=tol1 || fb==0.0){
				std::cout<<"Brent's Method converged in "<<iter<<" iterations\n";
				return b;
				}*/
				if (fabs(e) >= tol1 && fabs(fa)>fabs(fb)) {
					s = fb / fa;
					if (a == c) {
						p = 2.0*xm*s;
						q = 1.0 - s;
					}
					else {
						q = fa / fc;
						r = fb / fc;
						p = s * (2.0*xm*q*(q - r) - (b - a)*(r - 1.0));
						q = (q - 1.0)*(r - 1.0)*(s - 1.0);
					}
					if (p>0.0) q = -q;
					p = fabs(p);
					min1 = 3.0*xm*q - fabs(tol1*q);
					min2 = fabs(e*q);
					if (2.0*p<std::min(min1, min2)) {
						e = d;
						d = p / q;
					}
					else {
						d = xm;
						e = d;
					}
				}
				else {
					d = xm;
					e = d;
				}
				a = b;
				fa = fb;
				if (fabs(d)>tol1) {
					b += d;
				}
				else {
					b += template_funcs::SIGN(tol1, xm);
				}
				fb = eigeneqn_3(b, t, mm);
			}
			std::cerr << "Maximum number of iterations exceeded in zbrent\n";
			return 0.0;
		}
		else {
			std::string reason = "Error: double slab_tl::zbrent(double x1,double x2,double tol,bool t,int mm)\n";
			if (!c1) reason += "Cannot have x1 = x2\nx1 = " + template_funcs::toString(x1) + ", x2 = " + template_funcs::toString(x1) + "\n";
			if (!c2) reason += "Desired tolerance is less than smallest allowable EPS\n";
			if (!c3) reason += "mm >= M\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

// Definitions for the four layer slab derived class
// Type 1



// Definitions for the four layer slab derived class
// Type 2