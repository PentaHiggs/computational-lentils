/* Basic program to explore the chaotic properties of the Duffing Oscillator.
 * Uses cern's ROOT library/software in order to draw figures/graphs.
 * Program CANNOT compile on it's own without modification.  However, all code used is 
 * valid C++ code, the file is just not an independent C program.
 */

// C library includes
#include<functional>
#include<cmath>
#include<iostream>
#include<vector>
#include<algorithm>

// ROOT includes
#include"TROOT.h"
#include"TF1.h"
#include"TCanvas.h"
#include"TGraph.h"
#include"TMultiGraph.h"

// My integration routine
#include"runge_kutta1.cpp"

std::vector<double> duffing_eq(double t, const std::vector<double> &x, double gamma, double f0, double w){
	// Duffing equation with m=1, a=1/4, b=1/2
	std::vector<double> ret;
	ret.push_back(x[1]);
	ret.push_back(-gamma*x[1] + .5*x[0] - 2*std::pow(x[0],3) + f0*std::cos(w*t));
	return ret;
}

void duffing_oscillator(){
	typedef std::function<std::vector<double>(double, std::vector<double>)> doubleFunc;
	typedef std::vector<std::vector<double> > doubleVec;

	// Recreation/graphing of the duff oscillator shown in the book
	
	const double tBounds[2] = {0., 100.};
	doubleFunc duffing_example = [](double t, const std::vector<double> &x) { 
		return duffing_eq(t, x, .01, 2.0, 2.4);};
	std::vector<double> kutta_t;
	doubleVec kutta_x;
	std::vector<double> x0 = {0.5,0};

	TCanvas* canvas = new TCanvas();
	TMultiGraph* mg = new TMultiGraph();

	if(rungeKutta1(duffing_example, x0, kutta_t, kutta_x, tBounds[0], tBounds[1], .001)){
		std::cout << "Runge-kutta on textbook duff example a success!" << std::endl;
		// Here, I have to extract the first coordinate out of every 2-vector in kutta_x, that is
		// x (as opposed to its derivative, the second coordinates).  Can probably be done with a custom
		// iterator that gets passed to TGraph instead, but not good enough at C++ to do that yet
		std::vector<double> actual_x; actual_x.reserve(kutta_x.size());
		std::vector<double> actual_dx; actual_dx.reserve(kutta_x.size());
		for(const std::vector<double> &vec : kutta_x) {
			actual_x.push_back(vec[0]);
			actual_dx.push_back(vec[1]);
		}
		/*
		TGraph* g1 = new TGraph(kutta_t.size(), kutta_t.data(), actual_dx.data());
		g1->SetLineColor(kRed);
		mg.add(g1);
		*/
		TGraph* g2 = new TGraph(kutta_t.size(), kutta_t.data(), actual_x.data());
		g2->SetLineColor(kBlue);
		mg->Add(g2);
		mg->Draw("AL");
	} else { 
		std::cout << "Runge-kutta has Failed : <" << std::endl;
	}
	return;
}
