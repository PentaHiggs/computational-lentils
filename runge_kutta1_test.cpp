/* Tests of the functionality of the functions in runge-kutta1.cpp */

// NOTE: Currently not made to be used outside of ROOT, will NOT compile in g++
//
#include<iostream>
#include<cmath>
#include"runge_kutta1.cpp"
#include<vector>
#include<iomanip>

// ROOT dependencies
#include "TCanvas.h"
#include "TROOT.h"
#include "TF1.h"
//#include "TApplication.h"
#include "TGraph.h"
#include "TMultiGraph.h"

// Compile with g++ and additional option 'root-config --cflags --glibs'

void double_rungeKutta1_test_1() {
	/* Let's try this on a few functions.
	 * (d/dt)sin(3t) = 3cos(3t) = F(x,t) 
	 * (d/dt)cos(2t) = -2sin(2t) = F(x,t)
	 * Both with desired error e < 10^-4, on the
	 * interval [-5,5] (tested at endpoint and midpoint)
	 */

	std::vector<double> kutta_t;
	std::vector<double> kutta_x;
	const double tBounds[2] = {-5., 5};
	const double e = .0001;

	// Defining the two Fs and the intial values for the two testing functions
	std::function<double(double)> F1([](double){return 3.*std::cos(3*t);});
	const double x01 = std::sin(3*tBounds[0]);
	std::function<double(double)> F2([](double){return -2.*std::sin(2*t);});
	const double x02 = std::cos(2*tBounds[0]);

	TCanvas* canvas = new TCanvas();
	TMultiGraph* mg = new TMultiGraph();
	mg->setTitle("Runge-kutta approximation of sin(3t), cos(2t)");
	if (rungeKutta1(F1 , x01, &kutta_t, &kutta_x, tBounds[0], tBounds[1], e);) {
		std::cout << "Runge-kutta sin(3t) Successful!" << std::endl;
		TGraph* rungeGraph(kutta_t.size(), kutta_t.data(), kutta_x.data());
		mg->Add(rungeGraph);
	} else {
		std::cout << "Runge-kutta sin(3t) Failed" << std::endl;
	} 
	// Empty out the vectors first of old data
	kutta_t.clear(); kutta_x.clear();

	if (rungeKutta1(F2 , x02, &kutta_t, &kutta_x, tBounds[0], tBounds[1], e);) {
		std::cout << "Runge-kutta cos(2t) Successful!" << std::endl;
		TGraph* rungeGraph2(kutta_t.size(), kutta_t.data(), kutta_x.data());
		mg->Add(rungeGraph2);
	} else {
		std::cout << "Runge-kutta cos(2t) Failed" << std::endl;
	}
	return;
}




void double_exx1_test_1(){
	/* Let's first test this on sin(t)
	 * We're going to graph a computer-evaluated sin(t) and a kutta-evaluated 
	 * x=sin(t) from t=-5 to t=5.  (d/dt)sin(t) = cos(t)
	 */
	
	std::vector<double> kutta_t;
	std::vector<double> kutta_x;

	const float tBounds[2] = {0,100};

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
	sine->SetNpx(1000);
	sine->DrawClone("Same");	

	return;
}

void double_exx1_test_3() {
	std::vector<double> kutta_t;
	std::vector<double> kutta_x;

	const double tBounds[2] = {-5.,5.};

	auto poly = [](double x, double t)->double{return t*t*3.;};
	double(*F)(double, double) = poly;
	const double h0 = .1;
	double x = -125.;
	for(double t = tBounds[0]; t < tBounds[1] ; t += h0) {
		kutta_t.push_back(t);
		kutta_x.push_back(x);
		x = exx1(h0, x, t, F);
	}
	TGraph rungeGraph(int(kutta_t.size()), kutta_t.data(), kutta_x.data());

	TCanvas* canvas = new TCanvas();

	rungeGraph.SetTitle("Runge-Kutta approximation of t^3");
	rungeGraph.SetLineColor(kBlue);
	rungeGraph.DrawClone("AL");
	std::cout << "At four, we have " << rungeGraph.Eval(4,0,"");
	return;
}

void double_exx1_test_2(){
	/* Let's try it on the exponential function now!*/
	std::vector<double> kutta_t;
	std::vector<double> kutta_x;

	const float tBounds[2] = {-5.,30.};

	auto exponential = [](double x, double t)->double { return std::exp(t); };
	double(*F)(double, double) = exponential;

	const double h0 = .1;

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
	std::cout << std::endl;
	std::cout << "Run test 3? (y/n): ";
	std::cin >> yesno;
	if (yesno=='y') {
		double_exx1_test_3();
	}
}
