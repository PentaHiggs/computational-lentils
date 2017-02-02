#include<iostream>
#include<cfloat>
#include<typeinfo>

using namespace std;

int main() {
	cout << DBL_EPSILON << " is the amount of epsilon" << endl;
	double boop = 1.+DBL_EPSILON;
	cout << "Type of literal is: " << typeid(313.).name() << endl;
	if ( boop == 1. ){
		cout << " WOMBO COMBO !!! " << endl;
	}
	return 0;
}
