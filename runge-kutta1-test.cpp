/* Tests of the functionality of the functions in runge-kutta1.cpp */

#include<iostream>
#include<cmath>
#include"runge-kutta1.cpp"
#include<vector>
#include<iomanip>

void double_exx1_test_1(){
	/* Let's first test this on sin(t) */
	std::vector<double> kutta_values;
	std::vector<double> actual_values;

	double h0 = .01;
	double t0 = -5.0;
	double x0 = std::sin(-5.0);

	/* sin(t) satisifies d(sin t) = cos t , thus F = cos(t) */
	auto eff = [](double x, double t)->double{ return std::cos(t); };
	double(* F)(double, double) = eff;

	kutta_values.push_back(x0);
	
	for (int i = 0; i < 100; i++){
		x0 = exx1(h0, x0, t0, eff);
		t0 += h0;
		kutta_values.push_back(x0);
	}
	
	t0 = -5.0;
	x0 = std::sin(-5.0);
	actual_values.push_back(x0);
	for (int i = 0; i < 101; i++){
		actual_values.push_back(std::sin(t0));
		t0 += h0;
	}

	std::cout << "actual values | kutta values | diff" << std::endl;
	int theSmallerSize = (actual_values.size() < kutta_values.size() ? actual_values.size() : kutta_values.size() );
	for (int i = 0; i < theSmallerSize; i++){
		std::cout << std::setw(13) << actual_values.at(i) << " | " << std::setw(12) << kutta_values.at(i);
		std::cout << " | " << actual_values.at(i) - kutta_values.at(i) << std::endl;
	}

	return;
}

void double_exx1_test_2(){
	/* Let's try it on the exponential function now!*/
	std::vector<double> kutta_values;
	std::vector<double> actual_values;

	double h0 = .01;
	double t0 = -5.0;
	double x0 = std::exp(-5.0);

	/* exp(t) satisifies d(exp t) = exp t , thus F = exp(t) */
	auto eff = [](double x, double t)->double{ return std::exp(t); };
	double(* F)(double, double) = eff;

	kutta_values.push_back(x0);
	
	for (int i = 0; i < 100; i++){
		for(int j = 0; j < 100; j++){
			x0 = exx1(h0, x0, t0, eff);
			t0 += h0;
		}
		kutta_values.push_back(x0);
	}
	
	t0 = -5.0;
	x0 = std::exp(-5.0);
	actual_values.push_back(x0);
	for (int i = 0; i < 101; i++){
		actual_values.push_back(std::exp(t0));
		t0 += h0*10;
	}

	std::cout << "actual values | kutta values | diff" << std::endl;
	int theSmallerSize = (actual_values.size() < kutta_values.size() ? actual_values.size() : kutta_values.size() );
	for (int i = 0; i < theSmallerSize; i++){
		std::cout << std::setw(13) << actual_values.at(i) << " | " << std::setw(12) << kutta_values.at(i);
		std::cout << " | " << actual_values.at(i) - kutta_values.at(i) << std::endl;
	}

	return;
}

int main(){
	double_exx1_test_1();
	double_exx1_test_2();
	std::cout << "Welcome!";
	std::cin;
	return 0;
}
