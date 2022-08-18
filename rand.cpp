#include <iostream>
#include <ctime>
#include <random>
#include "tissue.h"
using namespace std;

int Tissue::unifRandInt(int a,int b) {
	//std::random_device rd;
	//std::mt19937 gen(rd());
	std::uniform_int_distribution<> dis(a,b);
	return dis(this->gen);	
}
double Tissue::unifRand() {
	std::uniform_real_distribution<> dis(0,1);
	return dis(this->gen);
}
double Tissue::unifRand(double a, double b) {
	std::uniform_real_distribution<> dis(a,b);
	//return (b-a)*unifRand() + a;
	return dis(this->gen);
}

