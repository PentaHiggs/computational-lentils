/* Interface for the ROOT graphincs package, coded in order to reduce the amount of boilerplate code
 * required for the simple graphing of results obtained from numerical algorithms I code
 */

#include<vector>
#include<array>
#include<string>

// ROOT dependencies
#include "TCanvas.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TMultiGraph.h"

#include "ROOT_interface.h"

MyGraph::MyGraph(std::string name) : current_graph_color(2), name(name+";t;x"){
	canvas = new TCanvas();
	graphs = new TMultiGraph();
	graphs->SetTitle(name.c_str());
	current_graph_color = 2;
}
void MyGraph::AddGraph(const std::vector<double>& x, const std::vector<double>& y){
	TGraph* graph = new TGraph(x.size(), x.data(), y.data());
	++current_graph_color;
	graph->SetLineColor(current_graph_color);
	graphs->Add(graph);
	graphs->Draw("AL");
}

void MyGraph::AddGraph(const std::vector<double>& x, const std::vector<std::vector<double> >& y, int n) {
	std::vector<double> newY;
	newY.reserve(y.size());
	for(std::vector<double> doub : y) { newY.push_back(doub[n]); }
	AddGraph(x, newY);
}

void MyGraph::AddGraph(const std::vector<double>& x, const std::vector<std::vector<double> >& y) {
	size_t size = y[0].size();
	for (size_t i = 0; i < size; ++i) { AddGraph(x, y, i); }
}
