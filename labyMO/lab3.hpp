#include <iostream>
#include <math.h>
#include <cmath>
#include <iomanip>


namespace lab3
{
	using namespace std;

	double picard(double (*fun)(double), double (*phi)(double), double xStart, int maxN, double tolx, double tolf);
	void picardTest();

	double bisection(double (*fun)(double), double jumpAmountOfAB, int numberOfJumps, int maxN, double tolx, double tolf);
	void bisectionTest();

	double newton(double (*fun)(double), double (*funDev)(double), double xStart, int maxN, double tolx, double tolf);
	void newtonTest();

	double secant(double (*fun)(double), double xStart1, double xStart2, int maxN, double tolx, double tolf);
	void secantTest();

	void test();
	

}