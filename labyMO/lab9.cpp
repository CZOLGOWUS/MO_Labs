#include "lab9.hpp"

namespace lab9
{
	using namespace std;

#define func(x)  ( exp( 2.0-2.0*x ) - 4.0*exp( 4.0 - x * 2.0 )+ 4.0*exp( x * 2.0 ) - exp(2.0 + 2.0 * x) - x + x * exp(4.0)) / (4.0 - 4 * exp(4.0))

#define p 1.0
#define q 0.0
#define r -4.0

#define alpha 0.0
#define beta 1.0
#define gamma -1.0

#define fi 0.0
#define psi 1.0
#define theta 0.0

#define xStart 0.0
#define xEnd 1.0


	void printMatrix(double** a, int row, int c)
	{
		for (int i = 0; i < row; i++)
		{
			for (int j = 0; j < c; j++)
				cout << setw(7) << a[i][j];

			cout << endl;
		}
	}

	void printVector(double* v,  int size)
	{
		for (int i = 0; i < size; i++)
			cout << setw(10) << v[i] << endl;
	}

	double** createMatrix( int numberOfRows,  int numberOfColls)
	{
		double** A = new double* [numberOfRows];

		for (int i = 0; i < numberOfRows; i++)
			A[i] = new double[numberOfColls];

		return A;
	}

	double* createVector( int size)
	{
		return new double[size];
	}

	void applyThomasAlgo(double* low, double* mid, double* up, double* B, double* x,  int size)
	{
		double* tempMid = new double[size];
		double* tempB = new double[size];
		tempMid[0] = mid[0];
		tempB[0] = B[0];


		for (int i = 1; i < size; i++)
		{
			tempMid[i] = mid[i] - low[i - 1] * (up[i - 1] / tempMid[i - 1]);
			tempB[i] = B[i] - (low[i - 1] * tempB[i - 1]) / tempMid[i - 1];
		}

		//for (int i = 1; i < size; i++)
		//	tempB[i] = B[i] - low[i - 1] * tempB[i - 1] / tempMid[i - 1];

		x[size - 1] = tempB[size - 1] / tempMid[size - 1];

		for (int i = size - 2; i >= 0; i--)
			x[i] = (tempB[i] - up[i] * x[i + 1]) / tempMid[i];

		delete[] tempMid;
		delete[] tempB;
	}


	double numerowMethod( double h , int n)
	{
		double xn = xStart, xi = xStart;
		double maxError = 0.0;

		double* lower = new double[n];
		double* mid = new double[n];
		double* upper = new double[n];
		double* B = new double[n];
		double* x = new double[n];
		double* error = new double[n];

		//first iter
		upper[0] = alpha / h;
		mid[0] = beta - alpha / h;
		B[0] = -gamma;

		double xnh;

		//iterate over every row in matrix:
		for (int i = 1; i < n - 1; i++)
		{
			xnh = i * h;
			//lower and upper have n-1 size in reality
			lower[i - 1] = p / (h * h) + r / (12.0);
			mid[i] = (-2.0 * p) / (h * h) + r * (10.0/12.0);
			upper[i] = p / (h * h) + r / (12.0);
			B[i] = (xn + xnh - h) / 12.0 + (10.0 / 12.0) * (xn + xnh) + (xn + xnh + h) / 12.0;
		}
		lower[n - 2] = -fi / h;
		mid[n - 1] = -fi / h + psi;
		B[n - 1] = -theta;

		applyThomasAlgo(lower,mid,upper,B,x,n);

		//calculate Error
		for (int i = 0; i < n; i++)
		{
			error[i] = abs(x[i] - func(xn));

			xn = xn + h;


			if (abs(error[i]) > maxError)
				maxError = error[i];
		}

		//todo: save to file
		if (n == 30)
		{
			fstream file;
			file.open("C:\\Users\\CZOLG\\Desktop\\numerowMethod.txt");
			
			cout << "Numerow Method\n";
			cout << setw(16) << "xi" << setw(16) << "x[i]" << setw(16) << "function(xi)" << endl;

			if (file.bad()) throw exception("file error (bad)");
			
			for (int i = 0; i < n; i++)
			{
				file << xi << " " << x[i] << " " << func(xi) << endl;
				cout <<setw(3) << i << setw(16) << scientific << xi << setw(16) << x[i] << setw(16) << func(xi) << endl;
				xi += h;
			}

			file.close();
		}

		delete[] lower;
		delete[] mid;
		delete[] upper;
		delete[] B;
		delete[] x;
		delete[] error;


		return maxError;
	}


	double conventionalMethod( double h, int n)
	{
		double xn = xStart, xi = xStart;
		double maxError = 0.0;

		double* lower = new double[n];
		double* mid = new double[n];
		double* upper = new double[n];
		double* B = new double[n];
		double* x = new double[n];
		double* error = new double[n];

		//first iter
		upper[0] = alpha / h;
		mid[0] = beta - alpha / h;
		B[0] = -gamma;

		double xnh;

		//iterate over every row in matrix:
		for (int i = 1; i < n - 1; i++)
		{
			xnh = h * i;
			//lower and upper have n-1 size in reality
			lower[i - 1] = p / (h * h) - q / (2.0 * h);
			mid[i] = r + (-2.0*p) / (h*h);
			upper[i] = p / (h * h) + q / (2.0*h);
			B[i] = (xn + xnh);
		}
		lower[n - 2] = -fi / h;
		mid[n - 1] = -fi / h + psi;
		B[n - 1] = -theta;

		applyThomasAlgo(lower, mid, upper, B, x, n);

		//calculate Error
		for (int i = 0; i < n; i++)
		{
			error[i] = abs(x[i] - func(xn));
			if (abs(error[i]) > maxError)
				maxError = error[i];

			xn = xn + h;
		}

		//todo: save to file
		if (n == 30)
		{
			fstream file;
			file.open("C:\\Users\\CZOLG\\Desktop\\conventionalMethod.txt",  std::ios::out);
			if (file.bad()) throw exception("file bad");

			cout << "conventionalMethod\n";
			cout << setw(16) << "xi" << setw(16) << "x[i]" << setw(16) << "function(xi)" << endl;

			for (int i = 0; i < n; i++)
			{
				file << xi << " " << x[i] << " " << func(xi) << endl;
				cout << setw(3) << i << setw(16) << scientific << xi << setw(16) << x[i] << setw(16) << func(xi) << endl;
				xi += h;
			}

			file.close();
		}

		delete[] lower;
		delete[] mid;
		delete[] upper;
		delete[] B;
		delete[] x;
		delete[] error;


		return maxError;
	}


	void test()
	{
		
		fstream file;
		file.open("C:\\Users\\CZOLG\\Desktop\\errors.txt",std::ios::out);
		file << scientific;

		double h;

		for (int i = 10; i < 200000; i += 10)
		{
			h = (xEnd - xStart) / (i - 1);
			file << setw(16) << log10(h) << " " << setw(16) << log10(conventionalMethod(h, i)) << " " << setw(16) << log10(numerowMethod(h, i)) << endl;
		}

		file.close();

	}

}