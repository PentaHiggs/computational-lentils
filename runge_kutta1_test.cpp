/* Tests of the functionality of the functions in runge-kutta1.cpp */

// NOTE: Currently not made to be used outside of ROOT, will NOT compile in g++
//
#include<iostream>
#include<cmath>
#include"runge-kutta1.cpp"
#include<vector>
#include<iomanip>

// ROOT dependencies
#include "TCanvas.h"
#include "TROOT.h"
#include "TF1.h"
#include "TApplication.h"
#include "TGraph.h"

// Compile with g++ and additional option 'root-config --cflags --glibs'

void double_exx1_test_1(){
	/* Let's first test this on sin(t)
	 * We're going to graph a computer-evaluated sin(t) and a kutta-evaluated 
	 * x=sin(t) from t=-5 to t=5.  (d/dt)sin(t) = cos(t)
	 */
	
	std::vector<double> kutta_t;
	std::vector<double> kutta_x;

	const float tBounds[2] = {-5.,5.};

	auto cosine = [](double x, double t)->double { return std::cos(t); };
	double(*F)(double, double) = cosine;

	const double h0 = 1;

	double x = std::sin(tBounds[0]);
	for(double t = tBounds[0]; t < tBounds[1] ; t += h0) {
		kutta_t.push_back(t);
		kutta_x.push_back(x);
		x = exx1(h0, x, t, F);
	}
	TGraph rungeGraph(int(kutta_t.size()), kutta_t.data(), kutta_x.data());
	
	rungeGraph.SetTitle("Runge-Kutta approximation of cosine");
	rungeGraph.SetLineColor(kBlue);

	// The canvas on which the runge-kutta and actual functions will be plotted
	TCanvas* canvas = new TCanvas();

	rungeGraph.DrawClone("AL");

	
	// Drawing an actual cosine
	TF1* sine = new TF1("Fsine", "sin(x)", tBounds[0], tBounds[1]);
	sine->SetLineColor(kRed);
	sine->DrawClone("Same");	

	return;
}


void double_exx1_test_2(){
	/* Let's try it on the exponential function now!*/
	std::vector<double> kutta_t;
	std::vector<double> kutta_x;

	const float tBounds[2] = {-5.,20.};

	auto exponential = [](double x, double t)->double { return std::exp(t); };
	double(*F)(double, double) = exponential;

	const double h0 = 1;

	double x = std::exp(tBounds[0]);
	for(double t = tBounds[0]; t < tBounds[1] ; t += h0) {
		kutta_t.push_back(t);
		kutta_x.push_back(x);
		x = exx1(h0, x, t, F);
	}
	TGraph rungeGraph(int(kutta_t.size()), kutta_t.data(), kutta_x.data());
	
	rungeGraph.SetTitle("Runge-Kutta approximation of exponential function");
	rungeGraph.SetLineColor(kBlue);

	// The canvas on which the runge-kutta and actual functions will be plotted
	TCanvas* canvas = new TCanvas();

	rungeGraph.DrawClone("AL");

	
	// Drawing an actual cosine
	TF1* sine = new TF1("Fsine", "exp(x)", tBounds[0], tBounds[1]);
	sine->SetLineColor(kRed);
	sine->DrawClone("Same");	

	return;
}

void runge_kutta1_test() {
	std::cout << "Run test 1? (y/n): ";
	char yesno;
	std::cin >> yesno;
	if (yesno=='y') {
		double_exx1_test_1();
	}
	std::cout << std::endl;
	std::cout << "Run test 2? (y/n): ";
	std::cin >> yesno;
	if (yesno=='y') {
		double_exx1_test_2();
	}
}
