/* Basic adaptive Runge-Kutta ODE algorithm implementation
 * Given an initial x0, y0, 
 * (all doubles), as well as the funcion F(x,y(x)) of the form shown above,
 * we will calculate k1, k2, k3, k4, and then finally y1, the next value of y */

#include<functional>
#include<cmath>
#include<algorithm>
#include<iterator>
#include<vector>

// Overloading vector addition and scalar multiplication.  First one taken from stackoverflow question 3376124
template <typename T>
std::vector<T> operator+(const std::vector<T> &a, const std::vector<T> &b) {
	assert(a.size() == b.size());
	std::vector<T> result;
	result.reserve(a.size());
	std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::plus<T>());
	return result;
}
template <typename T>
std::vector<T> operator*(T c, const std::vector<T> &v) {
	std::vector<T> res;
	res.reserve(v.size());
	std::transform(v.begin(), v.end(), std::back_inserter(res), [c](const T& x)->T{ return c*x; });
	return res;
}
	

double exx1(double h, double x, double y, const std::function<double(double, double)> &F){
	double kay1 = h*F(x, y);
	double kay2 = h*F(x + .5*h, y + .5*kay1);
	double kay3 = h*F(x + .5*h, y + .5*kay2);
	double kay4 = h*F(x + h, y + kay3);	
	return  y + kay1/6 + kay2/3 + kay3/3 + kay4/6;
}

// Overload of exx1 to take and return vector-valued y
std::vector<double> exx1(double h, double x, const std::vector<double> &y, \
	   	const std::function<std::vector<double>(double, std::vector<double>) > &F) {
	// Calculate all the kays, just like before thanks to operator overloading
	std::vector<double> kay1 = h*F(x, y);
	std::vector<double> kay2 = h*F(x + .5*h, y + .5*kay1);
	std::vector<double> kay3 = h*F(x + .5*h, y + .5*kay2);
	std::vector<double> kay4 = h*F(x + h, y + kay3);
	return y + (1./6.)*kay1 + (1./3.)*kay2 + (1./3.)*kay3 + (1./6.)*kay4;
}
	
bool rungeKutta1(const std::function<double(double, double)> &F, double y0, std::vector<double> &x, \
		std::vector<double> &y, double left, double right, double e) {
	// Debugging variable
	unsigned int evals = 0;

	// Setting a few parameters
	const double h0 = (right-left) / 1000.;
	const double hlow = h0/10;
	const double minDel = .0000001;
	const double rMax = 2.;
	
	double h = h0;

	// Core of the runge-Kutta evaluation process
	for(double x0=left; x0<right; x0+=h){
		x.push_back(x0);
		y.push_back(y0);
		
		double del = std::max( std::abs( exx1(h, x0, y0, F)-exx1(h/2, x0, y0, F) ) , minDel ); evals +=2;
		double ratio = (15./16.)*std::pow(e/del, 1./5.);
		// Don't want the new h to be more than rMax larger
		ratio = std::min(rMax, ratio);
		// If the ratio is zero, well we're kind of hopeless.  Abandon ship!
		if (ratio == 0) { return false; }

		h *= ratio;
		if (h < hlow) {h = hlow;}
		y0 = exx1(h, x0, y0, F); ++evals;
		if (evals > 5000000) break;
	}

	// We made it here, means we haven't failed!  Debug messages below.
	// std::cout << "We have done " << evals << " exx1 evaluations." << std::endl;
	std::cout << "We have made it to x=" << x.at(x.size()-1);
	return true;
}

// Overload of runge-kutta where it and F have a dependant variable that is a
// vector of doubles instead of just a scalar double.
bool rungeKutta1(const std::function<std::vector<double>(double, std::vector<double>)> &F, \
		std::vector<double> y0, std::vector<double> &x, \
		std::vector< std::vector<double> > &y, double left, double right, double e) {

	// Dimension of vectors being used.
	unsigned int dim = y0.size();

	// Debugging variable
	unsigned int evals = 0;

	// Setting a few necessary parameters
	const double h0 = (right - left) / 100000.;
	const double hlow = h0/1000;
	const double minDel = .00000001;
	const double rMax = 2.;

	double h = h0;

	// Core of the runge-Kutta evaluation process
	for(double x0=left; x0<right; x0+=h){
		x.push_back(x0);
		y.push_back(y0);

		std::vector<double> toDel = exx1(h, x0, y0, F) + (-1.)*exx1(h/2, x0, y0, F) ; evals+=2*dim;

		// I want to the greatest del to set my next h
		double del = minDel;
		std::for_each(toDel.begin(), toDel.end(), [&del](double &dub){ del=std::max(del, std::abs(dub)); });
		
		double ratio = std::min( (15./16.)*std::pow(e/del, 1./5.) , rMax);

		h *= ratio;
		if (h < hlow) {return false;}
		// Sometimes h reaches out of the range (left,right) because the function is so ridiculously easy to
		// evaluate without error.  This keeps it in there!
		if ( (right - x0) < h){ 
			h = right - x0 - hlow;
			if (h < hlow) {break;} // you're already so close to right.  We can stop now.
		}

		y0 = exx1(h, x0, y0, F); evals+=dim;
		if (evals > 4000000) break;
	}
		// We made it here!  This means catastrophic failure has been averted, hopefully.  Debug msgs below.
		// std::cout << "We have done " << evals << " exx1 evaluations." << std::endl;
		std::cout << "We have made it to x=" << x.at(x.size()-1) << " using " << evals \
			<< " evaluations!" << std::endl;
		std::cout << "The last h was " << h << std::endl;
		return true;
}
