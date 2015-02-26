/* Basic adaptive Runge-Kutta ODE algorithm implementation
 * Program solves differential equations of the form dx/dt = F(x(t),t)
 */

/* Given an initial x0, y0, 
 * (all doubles), as well as the funcion F(y(x), x) of the form shown above,
 * we will calculate k1, k2, k3, k4, and then finally y1, the next value of y */

#include<functional>
#include<cmath>

double exx1(double h, double x, double y, std::function<double(double, double)> F){
	double kay1 = h*F(x, y);
	double kay2 = h*F(x + .5*h, y + .5*kay1);
	double kay3 = h*F(x + .5*h, y + .5*kay2);
	double kay4 = h*F(x + h, y + kay3);	
	return  y + kay1/6 + kay2/3 + kay3/3 + kay4/6;
}

bool rungeKutta1(std::function<double(double, double)> F, double y0, std::vector<double>& x, \
		std::vector<double>& y, double left, double right, double e) {
	// Debugging variable
	unsigned int evals = 0;

	// Setting a few parameters
	const double h0 = (right-left) / 100000.;
	const double hlow = h0/10;
	const double minDel = .0000001;
	const double rMax = 2.;
	
	double h = h0;

	// Core of the runge-Kutta evaluation process
	for(double x0=left; x0<right; x0+=h){
		x.push_back(x0);
		y.push_back(y0);
		
		double del = std::max( std::abs( exx1(h, x0, y0, F)-exx1(h/2, x0, y0, F) ) , minDel ); evals +=3;
		double ratio = (15./16.)*std::pow(e/del, 1./5.);
		// Don't want the new h to be more than rMax larger
		ratio = std::min(rMax, ratio);
		// If the ratio is zero, well we're kind of hopeless.  Abandon ship!
		if (ratio == 0) { return false; }

		h *= ratio;
		if (h < hlow) {h = hlow;}
		y0 = exx1(h, x0, y0, F); ++evals;
		if (evals > 1000000) break;
	}

	// We made it here, means we haven't failed!  Debug messages below.
	// std::cout << "We have done " << evals << " exx1 evaluations." << std::endl;
	// std::cout << "We have made it to x=" << x.at(x.size()-1) << std::endl;
	return true;
}

