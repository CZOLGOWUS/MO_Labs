#include "lab7.hpp"
#include "Matrix.hpp"
#include <typeinfo>

namespace lab7
{
#define maxRows 4
#define maxColls 4

	void printMatrix(double** a,int r,int c)
	{
		for (int i = 0; i < r; i++)
		{
			for (int j = 0; j < c; j++)
				cout << setw(7) << a[i][j];

			cout << endl;
		}
	}

	void printVector(double* v, unsigned short size)
	{
		for (int i = 0; i < size; i++)
			cout << setw(10) << v[i] << endl;
	}

	double** createMatrix(unsigned short numberOfRows,unsigned short numberOfColls)
	{
		double** A = new double* [numberOfRows];

		for (int i = 0; i < numberOfRows; i++)
			A[i] = new double[numberOfColls];

		return A;
	}

	void fillMatrix(double** a, unsigned short numberOfRows, unsigned short numberOfColls, double filler)
	{
		for (int i = 0; i < numberOfRows; i++)
			for (int j = 0; j < numberOfColls; j++)
				a[i][j] = filler;
	}

	void fillVector(double* v,unsigned short size,double filler)
	{
		for (int j = 0; j < size; j++)
			v[j] = filler;
	}


	void jacobMethod(double** A, double* B, double* x0, double toly,double tolx, unsigned int n)
	{
		double rightSide = 0.0;

		double* tempX;
		double* xNext = new double[maxRows];
		double* xNow = new double[maxRows];
		xNow[0] = x0[0];
		xNow[1] = x0[1];
		xNow[2] = x0[2];
		xNow[3] = x0[3];


		double residuum[maxRows] = {0.0,0.0,0.0,0.0};
		bool found = false;
		int ns = 0;

		while (ns < n)
		{

			//-(L+U)*x0 + b
			for (int i = 0; i < maxColls; i++)
			{
				rightSide = 0.0;
				for (int j = 0; j < maxRows; j++)
				{
					if(j != i)
						rightSide += A[i][j] * xNow[j];
				}

				rightSide = -rightSide + B[i];

				//xNext aka. : kazdy z elementów (-(L+U)*x0 + b) dzielimy przez kazdy z elementów D
				xNext[i] = rightSide / A[i][i];
			}

			for (int i = 0; i < maxColls; i++)
			{
				residuum[i] = 0.0;

				for (int j = 0; j < maxRows; j++)
					residuum[i] += A[i][j] * xNow[j];

				residuum[i] -= B[i];
			}



			cout << setw(10) << "X" << "   " << setw(10) << "residuum" << "   " << setw(10) << "tolx"<< endl;
			cout << setw(10) << xNext[0]<< "   " << setw(10) << residuum[0] << "   " << setw(10) << abs(xNext[0] - xNow[0]) << endl;
			cout << setw(10) << xNext[1]<< "   " << setw(10) << residuum[1] << "   " << setw(10) << abs(xNext[1] - xNow[1]) << endl;
			cout << setw(10) << xNext[2]<< "   " << setw(10) << residuum[2] << "   " << setw(10) << abs(xNext[2] - xNow[2]) << endl;
			cout << setw(10) << xNext[3]<< "   " << setw(10) << residuum[3] << "   " << setw(10) << abs(xNext[3] - xNow[3]) << endl;
			cout << endl;


			if (
				residuum[0] <= toly
				&& residuum[1] <= toly
				&& residuum[2] <= toly
				&& residuum[3] <= toly
				&& abs(xNext[0] - xNow[0]) <= tolx
				&& abs(xNext[1] - xNow[1]) <= tolx
				&& abs(xNext[2] - xNow[2]) <= tolx
				&& abs(xNext[3] - xNow[3]) <= tolx
				)
			{
				found = true;
				break;
			}

			ns++;

			tempX = xNow;
			xNow = xNext;
			xNext = tempX;

		}

		if (found)
		{
			cout << "przyblizenie jest w granicach, numer iteracji: " << ns << endl;

			printVector(xNow, maxRows);
		}
		else
		{
			cout << "limit iteracji przekroczony\n";
		}

		delete[] xNow;
		delete[] xNext;
	}

	void GSmethod(double** A, double* B, double* x0, double toly, double tolx, unsigned int n)
	{

		double rightSide = 0;
		double sum = 0.0;

		double* tempX;
		double* xNext = new double[maxRows];
		double* xNow = new double[maxRows];
		xNow[0] = x0[0];
		xNow[1] = x0[1];
		xNow[2] = x0[2];
		xNow[3] = x0[3];


		double residuum[maxRows] = { 0.0,0.0,0.0,0.0 };
		bool found = false;
		int ns = 0;


		while (ns < n)
		{
			for (int i = 0; i < maxRows; i++)
			{
				rightSide = 0.0;

				for (int j = 0; j < maxRows; j++)
					if(j > i)
						rightSide += A[i][j] * xNow[j];

				rightSide = -rightSide + B[i];

				sum = 0.0;

				for (int j = 0; j < i; j++)
					sum += A[i][j] * xNext[j];

				xNext[i] = (rightSide - sum) / A[i][i];
			}

			for (int i = 0; i < maxColls; i++)
			{
				residuum[i] = 0;
				for (int j = 0; j < maxRows; j++)
					residuum[i] += A[i][j] * xNow[j];

				residuum[i] -= B[i];
			}


			cout << setw(10) << "X" << "   " << setw(10) << "residuum" << "   " << setw(10) << "tolx" << endl;
			cout << setw(10) << xNext[0] << "   " << setw(10) << residuum[0] << "   " << setw(10) << abs(xNext[0] - xNow[0]) << endl;
			cout << setw(10) << xNext[1] << "   " << setw(10) << residuum[1] << "   " << setw(10) << abs(xNext[1] - xNow[1]) << endl;
			cout << setw(10) << xNext[2] << "   " << setw(10) << residuum[2] << "   " << setw(10) << abs(xNext[2] - xNow[2]) << endl;
			cout << setw(10) << xNext[3] << "   " << setw(10) << residuum[3] << "   " << setw(10) << abs(xNext[3] - xNow[3]) << endl;
			cout << endl;


			if (
				residuum[0] <= toly
				&& residuum[1] <= toly
				&& residuum[2] <= toly
				&& residuum[3] <= toly
				&& abs(xNext[0] - xNow[0]) <= tolx
				&& abs(xNext[1] - xNow[1]) <= tolx
				&& abs(xNext[2] - xNow[2]) <= tolx
				&& abs(xNext[3] - xNow[3]) <= tolx
				)
			{
				found = true;
				break;
			}


			ns++;

			tempX = xNow;
			xNow = xNext;
			xNext = tempX;

		}

		if (found)
		{
			cout << "przyblizenie jest w granicach, numer iteracji: " << ns << endl;

			printVector(xNow, maxRows);
		} 
		else
		{
			cout << "limit iteracji przekroczony\n";
		}

		delete[] xNow;
		delete[] xNext;
	}


	void SORmethod(double** A, double* B, double* x0,double w, double toly, double tolx, unsigned int n)
	{

		double rightSide = 0;
		double sum = 0.0;

		double* tempX;
		double* xNext = new double[maxRows];
		double* xNow = new double[maxRows];
		xNow[0] = x0[0];
		xNow[1] = x0[1];
		xNow[2] = x0[2];
		xNow[3] = x0[3];


		double residuum[maxRows] = { 0.0,0.0,0.0,0.0 };
		bool found = false;
		int ns = 0;


		while (ns < n)
		{

			for (int i = 0; i < maxRows; i++)
			{
				rightSide = 0.0;
				for (int j = 0; j < maxRows; j++)
					if (j > i)
						rightSide += A[i][j] * xNow[j];
					else if(i == j)
						rightSide += (1 - (1/w))*A[i][j] * xNow[j];

				rightSide = -rightSide + B[i];

				sum = 0.0;

				for (int j = 0; j < i; j++)
						sum += A[i][j] * xNext[j];

				xNext[i] = (rightSide - sum) / (A[i][i]/w);

			}

			for (int i = 0; i < maxColls; i++)
			{
				residuum[i] = 0;
				for (int j = 0; j < maxRows; j++)
					residuum[i] += A[i][j] * xNow[j];

				residuum[i] -= B[i];
			}


			cout << setw(10) << "X" << "   " << setw(10) << "residuum" << "   " << setw(10) << "tolx" << endl;
			cout << setw(10) << xNext[0] << "   " << setw(10) << residuum[0] << "   " << setw(10) << abs(xNext[0] - xNow[0]) << endl;
			cout << setw(10) << xNext[1] << "   " << setw(10) << residuum[1] << "   " << setw(10) << abs(xNext[1] - xNow[1]) << endl;
			cout << setw(10) << xNext[2] << "   " << setw(10) << residuum[2] << "   " << setw(10) << abs(xNext[2] - xNow[2]) << endl;
			cout << setw(10) << xNext[3] << "   " << setw(10) << residuum[3] << "   " << setw(10) << abs(xNext[3] - xNow[3]) << endl;
			cout << endl;


			if (
				residuum[0] <= toly
				&& residuum[1] <= toly
				&& residuum[2] <= toly
				&& residuum[3] <= toly
				&& abs(xNext[0] - xNow[0]) <= tolx
				&& abs(xNext[1] - xNow[1]) <= tolx
				&& abs(xNext[2] - xNow[2]) <= tolx
				&& abs(xNext[3] - xNow[3]) <= tolx
				)
			{
				found = true;
				break;
			}


			ns++;

			tempX = xNow;
			xNow = xNext;
			xNext = tempX;

		}

		if (found)
		{
			cout << "przyblizenie jest w granicach, numer iteracji: " << ns << endl;

			printVector(xNow, maxRows);
		}
		else
		{
			cout << "limit iteracji przekroczony\n";
		}

		delete[] xNow;
		delete[] xNext;
	}



	void test()
	{

		double** A = createMatrix(maxRows, maxColls);
		double* B = new double[maxRows];
		double* x = new double[maxRows];

		B[0] = 116.0;
		B[1] = -226.0;
		B[2] = 912.0;
		B[3] = -1174.0;

		x[0] = 2.0;
		x[1] = 2.0;
		x[2] = 2.0;
		x[3] = 2.0;

		double temp[4][4] =
		{
			{100.0 , -1.0 , 2.0 , -3.0 },
			{1.0 , 200.0 , -4.0 , 5.0},
			{-2.0 , 4.0 , 300.0 , -6.0},
			{3.0 , -5.0 , 6.0 , 400.0}
		};

		for (int i = 0; i < maxRows; i++)
			for (int j = 0; j < maxColls; j++)
				A[i][j] = temp[i][j];

		cout << "A : \n";
		printMatrix(A,maxRows,maxColls);
		cout << endl;

		cout << "B : \n";
		printVector(B,maxRows);
		cout << endl;
		
		cout << "metoda Jacobiego:\n";
		jacobMethod(A,B,x,0.0001,0.0001,100);
		cout << endl << endl;
		
		cout << "metoda Gaussa-Seidela:\n";
		GSmethod(A, B, x, 0.0001, 0.0001, 100);
		cout << endl << endl;
		
		cout << "metoda SOR:\n";
		SORmethod(A, B, x, 0.5, 0.0001, 0.0001, 100);

	}

}