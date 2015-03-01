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
	ret.push_back(-gamma*x[1] + .5*x[0] - 2.*std::pow(x[0],3.) + f0*std::cos(w*t));
	return ret;
}

void duffing_oscillator(){
	typedef std::function<std::vector<double>(double, std::vector<double>)> doubleFunc;
	typedef std::vector<std::vector<double> > doubleVec;

	// Recreation/graphing of the duff oscillator shown in the book
	
	const double tBounds[2] = {0., 200.};
	doubleFunc duffing_example = [](double t, const std::vector<double> &x) { 
		return duffing_eq(t, x, 0.1, 2.0, 2.4);};
	std::vector<double> kutta_t;
	doubleVec kutta_x;
	std::vector<double> x0 = {0.5,0};
	std::vector<double> x1 = {0.5,0};
	TCanvas* canvas = new TCanvas();
	TMultiGraph* mg = new TMultiGraph();

	if(rungeKutta1(duffing_example, x0, kutta_t, kutta_x, tBounds[0], tBounds[1], .0000001)){
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
		mg->Add(g1);
		*/
		TGraph* g2 = new TGraph(kutta_t.size(), kutta_t.data(), actual_x.data());
		g2->SetLineColor(kBlue);
		mg->Add(g2);
	} else { 
		std::cout << "Runge-kutta has Failed : <" << std::endl;
	}
	
	kutta_t.clear();
	kutta_x.clear();
	
	std::vector<double> actual_x; actual_x.reserve(kutta_x.size());
	std::vector<double> actual_dx; actual_dx.reserve(kutta_x.size());

	if(rungeKutta1(duffing_example, x1, kutta_t, kutta_x, tBounds[0], tBounds[1], .0001)) {

		for(const std::vector<double> &vec : kutta_x) {
			actual_x.push_back(vec[0]);
			actual_dx.push_back(vec[1]);
		}
		TGraph* g3 = new TGraph(kutta_t.size(), kutta_t.data(), actual_x.data());
		g3->SetLineColor(kGreen);
		mg->Add(g3);
	} else { std::cout << "Runge-kutta has failed for the second one" << std::endl; }
	
	mg->Draw("AL");

	const double pie = std::atan(1);

	// Routine to graph the strange attractor of the duffing oscillator
	TCanvas* canvas2 = new TCanvas("cStrange1", "Strange Attractor for duffing oscillator");
	double tee0 = 2.*pie/2.4;
	double tee = 0;

	std::vector<double> eks; eks.reserve(std::ceil(100./tee));
	std::vector<double> dee_eks; dee_eks.reserve(std::ceil(100./tee));

	for (int i = 0; i < kutta_t.size(); ++i){
		if (kutta_t.at(i) >= tee) {
			eks.push_back(kutta_x.at(i)[0]);
			dee_eks.push_back(kutta_x.at(i)[1]);
			tee += tee0;
		}
	}

	TGraph* g4 = new TGraph(eks.size(), eks.data(), dee_eks.data());
	g4->GetXaxis()->SetTitle("x");
	g4->GetYaxis()->SetTitle("x_dot");
	g4->GetXaxis()->CenterTitle();
	g4->GetYaxis()->CenterTitle();
	g4->DrawClone("A*");

	// Draw the phase space diagram
	TCanvas* canvas3 = new TCanvas("cPhase1", "Phase diagram for duffing oscillator");
	TGraph* g5 = new TGraph(actual_x.size(), actual_x.data(), actual_dx.data());
	g5->GetXaxis()->SetTitle("x");
	g5->GetYaxis()->SetTitle("x_dot");
	g5->GetXaxis()->CenterTitle();
	g5->GetYaxis()->CenterTitle();
	g5->DrawClone("AL");
	return;
}
