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

void bool_rungeKutta1_test_1() {
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
	std::function<double(double, double)> F1([](double t, double x){ \
			return 3.*std::cos(3*t);});
	const double x01 = std::sin(3*tBounds[0]);
	std::function<double(double, double)> F2([](double t, double x){ \
			return -2.*std::sin(2*t);});
	const double x02 = std::cos(2*tBounds[0]);

	TCanvas* canvas = new TCanvas();
	TMultiGraph* mg = new TMultiGraph();
	//mg->setTitle("Runge-kutta approximation of sin(3t)[red], cos(2t)[blue];x;y");
	if (rungeKutta1(F1 , x01, kutta_t, kutta_x, tBounds[0], tBounds[1], e)) {
		std::cout << "Runge-kutta sin(3t) Successful!" << std::endl;
		TGraph* rungeGraph = new TGraph(kutta_x.size(), kutta_t.data(), kutta_x.data());
		rungeGraph->SetLineColor(kRed);
		mg->Add(rungeGraph);
	} else {
		std::cout << "Runge-kutta sin(3t) Failed" << std::endl;
	} 
	// Empty out the vectors first of old data
	kutta_t.clear(); kutta_x.clear();

	if (rungeKutta1(F2 , x02, kutta_t, kutta_x, tBounds[0], tBounds[1], e)) {
		std::cout << "Runge-kutta cos(2t) Successful!" << std::endl;
		TGraph* rungeGraph2 = new TGraph(kutta_x.size(), kutta_t.data(), kutta_x.data());
		rungeGraph2->SetLineColor(kBlue);
		mg->Add(rungeGraph2);
	} else {
		std::cout << "Runge-kutta cos(2t) Failed" << std::endl;
	}
	mg->Draw("AL");
	return;
}

void bool_rungeKutta1_test_2(){
	/* Here, we test the overloaded version of runge-kutta that takes in and spits out vectors.
	 * (d/dt){x_1,x_2} = {x_2,-x_1} This should produce a pair of trigs, x_1 being A*cos(t+o) and
	 * x_2 being -A*sin(t+o), where o = arcsin(-x_1/A) = arccos(x_2/A).
	 * Lets hope for the best!
	 */  
	typedef std::function<std::vector<double>(double, std::vector<double>)> myFunc;
	typedef std::vector<std::vector<double> > doubleVec;

	std::vector<double> kutta_t;
	doubleVec kutta_x;

	TCanvas* canvas = new TCanvas();
	TMultiGraph* mg = new TMultiGraph;
	mg->SetTitle("Runge-kutta(vector) approximations of trig functions");
	
	myFunc F = [](double t, const std::vector<double> &x) {
		std::vector<double> retVal;
		retVal.push_back(x[1]);
		retVal.push_back(-x[0]);
		return retVal;
	};

	// This starting condition, x=x_dot=0 should give us the trivial solution where x_1(t)=x_2(t)=0
	std::vector<double> x0_trivial = {0.,0.};
	// This starting condition, x_dot=1, x=0, should give us x_1=sin(t)
	std::vector<double> x0_sine = {std::sin(-5.), std::cos(-5.)};
	// This starting condition should give us x_1 = cos(t)
	std::vector<double> x0_cosine = {std::cos(-5.),-std::sin(-5)};
	// Lets stick them all in one place for ease of execution
	doubleVec all_x0;
	all_x0.push_back(x0_trivial);
  	all_x0.push_back(x0_sine); 
	all_x0.push_back(x0_cosine);

	double theBounds[2] = {-5., 5.};
	
	unsigned int lineColor = 0;

	// God knows what type kBlue, kRed, kGreen are (probably strings/structs/int constants or something)
	std::vector<decltype(kBlue)> colors_woo = {kBlue, kRed, kGreen};

	for (std::vector<double> x0 : all_x0){
		kutta_t.clear();
		kutta_x.clear();
		if(rungeKutta1(F, x0, kutta_t, kutta_x, theBounds[0], theBounds[1], .0001)){
			// Extract the x first
			std::vector<double> ecks;
			ecks.reserve(kutta_x.size());
			for(std::vector<double> x : kutta_x){ ecks.push_back(x[0]); }

			// Now graph it
			TGraph* rungeGraph = new TGraph(kutta_x.size(), kutta_t.data(), ecks.data());
			rungeGraph->SetLineColor(colors_woo[lineColor]);
			mg->Add(rungeGraph);
		} else {
			std::cout << "Runge-kutta has failed!" << std::endl;
			return;
		}
		lineColor++;
	}
	mg->Draw("AL");
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

	auto cosine = [](double t, double x)->double { return std::cos(t); };
	double(*F)(double, double) = cosine;

	const double h0 = 1;

	double x = std::sin(tBounds[0]);
	for(double t = tBounds[0]; t < tBounds[1] ; t += h0) {
		kutta_t.push_back(t);
		kutta_x.push_back(x);
		x = exx1(h0, t, x, F);
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

	auto poly = [](double t, double x)->double{return t*t*3.;};
	double(*F)(double, double) = poly;
	const double h0 = .1;
	double x = -125.;
	for(double t = tBounds[0]; t < tBounds[1] ; t += h0) {
		kutta_t.push_back(t);
		kutta_x.push_back(x);
		x = exx1(h0, t, x, F);
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

	auto exponential = [](double t, double x)->double { return std::exp(t); };
	double(*F)(double, double) = exponential;

	const double h0 = .1;

	double x = std::exp(tBounds[0]);
	for(double t = tBounds[0]; t < tBounds[1] ; t += h0) {
		kutta_t.push_back(t);
		kutta_x.push_back(x);
		x = exx1(h0, t, x, F);
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

	std::cout << std::endl;
	std::cout << "Run rungeKutta test 1? (y/n): ";
	std::cin >> yesno;
	if (yesno=='y') {
		bool_rungeKutta1_test_1();
	}

	std::cout << std::endl;
	std::cout << "Run rungeKutta test 2? (y/n): ";
	std::cin >> yesno;
	if (yesno=='y') {
		bool_rungeKutta1_test_2();
	}
}
