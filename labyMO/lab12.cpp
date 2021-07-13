#include "lab12.hpp"

namespace lab12
{
#define PI 3.141592653589793238462643383279502884L

#define leftMax  -1.0
#define rightMax 1.0
#define sizeOfNodeGrid 20
#define analytic(x) 1.0 / (1.0 + 10.0 * x * x * x * x * x * x)


	double** createMatrix(int numberOfRows, int numberOfColls)
	{
		double** A = new double* [numberOfRows];

		for (int i = 0; i < numberOfRows; i++)
		{
			A[i] = new double[numberOfColls];

			for (int j = 0; j < numberOfColls; j++)
				A[i][j] = 0.0;
		}
		return A;
	}

	double* createVector(int size)
	{
		double* v = new double[size];
		for (int i = 0; i < size; i++)
			v[i] = 0.0;

		return v;
	}


	double interpolateFromToVec(double* v1, double* v2,int vecSize, double interpolationParameter)
	{
		
		double interResult = v2[vecSize - 1] * (interpolationParameter - v1[vecSize - 2]) + v2[vecSize - 2]; // 9 8 8

		for (int i = vecSize - 2; i >= 1 ; i--) //8 7 6 5...
			interResult = interResult * (interpolationParameter - v1[i - 1]) + v2[i-1];


		return interResult;

	}


	double* calculateEqualDist(double l, double r, int sizeOfGrid)
	{
		double* results = createVector(sizeOfGrid);

		double xd = l;
		double h = (r - l) / (sizeOfGrid - 1);
		for (int i = 0; i < sizeOfGrid; i++, xd += h)
			results[i] = xd;
			
		return results;
	}


	double* calculateChebyshevNodes(double l, double r, int gridSize)
	{
		double* results = createVector(gridSize + 1);
		for (int i = 0; i < gridSize+1; i++)
		{
			results[i] = (l + r) / 2.0;
		}

		double xi;
		for (int i = 0; i < gridSize + 1; i++)
		{
			xi = cos((2.0 * i + 1.0) / (2.0 * gridSize + 2.0) * PI);
			results[i] += (r - l) / 2.0 * xi;
		}
		return results;
	}


	double* calculateNewtonNodes(double* xNow,double* vals,int vecSize)
	{
		double* temp = createVector(vecSize);
		double* calculatedVector = createVector(vecSize);
		int iterCount = 0;

		for (int i = 0; i < vecSize; i++)
			calculatedVector[i] = vals[i];


		for (int i = 1; i < vecSize; i++)
		{
			for (int i = 0; i < vecSize; i++)
				temp[i] = calculatedVector[i];

			iterCount = 0;
			for (int j = i; j < vecSize; j++,iterCount++)
				temp[j] = (calculatedVector[j] - calculatedVector[j - 1]) / (xNow[j] - xNow[iterCount]);


			for (int i = 0; i < vecSize; i++)
				calculatedVector[i] = temp[i];
		}

		return calculatedVector;

	}


	void test()
	{
		fstream wyniki;
		wyniki.open("C:\\Users\\CZOLG\\Desktop\\lab12\\wyniki.txt");

		if (wyniki.bad())
			throw exception("could not open file");

		//setup
		double* xEqualDist = calculateEqualDist(leftMax, rightMax, sizeOfNodeGrid);

		//analitic conv
		double* yEqualDist = createVector(sizeOfNodeGrid);
		for (int i = 0; i < sizeOfNodeGrid - 1; i++)
			yEqualDist[i] = analytic(xEqualDist[i]);
		

		//calculate
		double* calculatedNewton = calculateNewtonNodes(xEqualDist, yEqualDist,sizeOfNodeGrid);
		
		double* calculatedChebyshev = calculateChebyshevNodes(leftMax,rightMax,sizeOfNodeGrid);
		
		//analitic conv
		double* analyticChebyshev = createVector(sizeOfNodeGrid);
		for (int i = 0; i < sizeOfNodeGrid - 1; i++)
			analyticChebyshev[i] = analytic(calculatedChebyshev[i]);


		double* calculatedNewtonWithChebyshev = calculateNewtonNodes(calculatedChebyshev, analyticChebyshev,sizeOfNodeGrid);


		//get results
		double x = leftMax;
		x += 0.01;
		while (x < rightMax)
		{
			wyniki << scientific << setprecision(16)
				<< x << " "
				<< analytic(x) << " "
				<< interpolateFromToVec(xEqualDist, calculatedNewton,  sizeOfNodeGrid,x) << " "
				<< interpolateFromToVec(calculatedChebyshev, calculatedNewtonWithChebyshev, sizeOfNodeGrid, x) << endl;

			x += 0.02;
			
		}

		wyniki.close();


	}
}
