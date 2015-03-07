/* ROOT_interface.h */

#ifndef INC_ROOT_INTERFACE_H
#define INC_ROOT_INTERFACE_H

#include<vector>
#include<array>
#include<string>

#include"TCanvas.h"
#include"TMultiGraph.h"

class MyGraph {
	public: 
		TCanvas *canvas;
		TMultiGraph *graphs;
		std::string name;
		MyGraph(std::string);
		void AddGraph(const std::vector<double>&, const std::vector<double>&);
		void AddGraph(const std::vector<double>&, const std::vector<std::vector<double> >&);
		void AddGraph(const std::vector<double>&, const std::vector<std::vector<double> >&, int);
		template <int n>
		void AddGraph(const std::vector<double>&, const std::vector<std::array<double,n> >&, int);
		template <int n>
		void AddGraph(const std::vector<double>&, const std::vector<std::array<double,n> >&);
	private:
		int current_graph_color;
};

// Template functions

template<int n>
void MyGraph::AddGraph(const std::vector<double>& x,const std::vector<std::array<double,n> >& y, int k) {
	std::vector<double> newY;
	newY.reserve(y.size());
	for(std::array<double, n> doub : y) { newY.push_back(doub[n]); }
	AddGraph(x, newY);
}

template<int n>
void MyGraph::AddGraph(const std::vector<double>& x,const std::vector<std::array<double,n> >& y) {
	size_t size = y[0].size();
	for (size_t i = 0; i < size; ++i) { AddGraph(x, y, i); }
}

#endif /* INC_ROOT_INTERFACE_H */
