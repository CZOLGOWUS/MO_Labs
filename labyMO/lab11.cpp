#include "lab11.hpp"
#include "calerf.h"

namespace lab11
{

#define maxT 2.0
#define parameterD 1.0
#define h 0.07
#define xMax 6.0*sqrt(parameterD*maxT)
#define lambda 0.4

	using namespace std;

	double analitic(double x,double t)
	{
		return calerf::ERFCL(x / (2 * sqrt(parameterD * t)));
	}

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



	//Thomas + LDU ----------------------------
	void ThomasMethod(double* upper,double* mid,double* low, double* x,double* B, int size)
	{
		double* tb = createVector(size);
		double* tmid = createVector(size);

		tb[0] = B[0];
		tmid[0] = mid[0];

		for (int i = 1; i < size; i++)
		{
			tmid[i] = mid[i] - low[i - 1] * (upper[i - 1] / tmid[i - 1]);
			tb[i] = B[i] - low[i - 1] * tb[i - 1] / tmid[i - 1];
		}

		x[size - 1] = tb[size - 1] / tmid[size - 1];

		for (int i = size - 2; i >= 0; i--)
			x[i] = (tb[i] - upper[i] * x[i + 1]) / tmid[i];

		delete[] tb, tmid;
	}

	void LDUAlgo(double** org, double** u, double** m, double** l, int size)
	{
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				if (i == j)
					m[i][i] = org[i][i];
				else if (j < i)
					l[i][j] = org[i][j];
				else if (j > i)
					u[i][j] = org[i][j];
			}
		}
	}


	//new thomas
	
	void thomasReadyMatrix(double** org, double** result, int size)
	{
		result[0][0] = org[0][1];
		result[0][1] = org[0][2];
		for (int i = 1; i < size; i++)
		{
			result[i][0] = org[i][1] - (org[i][0] * (1.0 / result[i - 1][0]) * org[i - 1][2]);
			result[i][1] = org[i][2];
		}
	}
	
	
	void thomasReadyVector(double* orgVec, double* result, double** orgMatrix, double** resultMatrix, int size)
	{
		result[0] = orgVec[0];
		for (int i = 1; i < size; i++)
			result[i] = orgVec[i] - (orgMatrix[i][0] * (1.0 / resultMatrix[i - 1][0]) * result[i - 1]);
	}

	void ThomasAlgo(double** thomasReadyMatrix,double* thomasReadyVec,double* result,int size)
	{
		result[size - 1] = (1.0 / thomasReadyMatrix[size - 1][0]) * thomasReadyVec[size - 1];
		for (int i = size -2; i >= 0; i--)
		{
			result[i] = (1.0 / thomasReadyMatrix[i][0]) * (thomasReadyVec[i] - thomasReadyMatrix[i][1]*result[i+1]);
		}
	}



	//JACOBS ALGO:---------------

	


	void multiplyMatrixiesSameShape(double** a, double** b, int size)
	{
		double** result = createMatrix(size, size);
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				for (int k = 0; k < size; k++)
				{
					result[i][j] += a[i][k] * b[k][j];
				}
			}
		}
		for (int i = 0; i < size; i++)
		{
			delete[] result[i];
		}
		delete[] result;
	}

	double** multiplyMatrixiesSameShapeReturn(double** a, double** b, int size)
	{
		double** result = createMatrix(size, size);
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				for (int k = 0; k < size; k++)
				{
					result[i][j] += a[i][k] * b[k][j];
				}
			}
		}

		return result;
	}
	
	void multiVecMat(double* vec, double** m, int size)
	{
		double* result = createVector(size);
		double sum = 0.0;
	
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				sum += m[i][j] * vec[j];
			}
			result[i] = sum;
		}
		for (int i = 0; i < size; i++)
		{
			vec[i] = result[i];
		}
		delete[] result;
	}

	double* multiVecMatReturn(double* vec, double** m, int size)
	{
		double* result = createVector(size);
		double sum = 0.0;

		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				sum += m[i][j] * vec[j];
			}
			result[i] = sum;
		}

		return result;
	}
	
	void addVec(double* a, double* b,int size)
	{
		for (int i = 0; i < size; i++)
			a[i] += b[i];
	}
	
	bool testRes(double** A, double* b, double* x, int size)
	{
		double* resTemp = createVector(size);
		bool res = true;
		multiVecMat(resTemp, A, size);
		for (int i = 0; i < size; i++)
		{
			if (resTemp[i] > 0.0001)
			{
				res = false;
				break;
			}
		}
		delete[] resTemp;
		return res;
	}
	
	bool testEs(double* vecPrev, double* x, int size)
	{
		for (int i = 0; i < size; i++)
		{
			vecPrev[i] -= x[i];
			if (vecPrev[i] > 0.0001)
			{
				return false;
			}
		}
		return true;
	}
	
	void JacobiMethod(double** A, double* b, double** M, double* C, double* x, int size,int maxIter)
	{
		double* vecP = createVector(size);
	
		for (int i = 0; i < size; i++)
		{
			vecP[i] = x[i];
		}
	
		for (int i = 0; i < maxIter; i++)
		{
			multiVecMat(x,M,size);
			for (int j = 0; j < size; j++)
			{
				x[j] += C[j];
			}
			if (testRes(A,b,x,size) && testEs(vecP,x,size))
			{
				break;
			}
			else
			{
				for (int i = 0; i < size; i++)
				{
					vecP[i] = x[i];
				}
			}
	
		}
	
		delete[] vecP;
	
	}

	void CN(double** A, int size)
	{
		A[0][0] = 0.0;
		A[0][1] = 1.0;
		A[0][2] = 0.0;
		for (int i = 1; i < size - 1; i++)
		{
			A[i][0] = lambda * 0.5;
			A[i][1] = -lambda - 1.0;
			A[i][2] = lambda * 0.5;
		}


		A[size - 1][0] = 0.0;
		A[size - 1][1] = 1.0;
		A[size - 1][2] = 0.0;
	}

	void classicDirectMethod()
	{
		fstream errXtime,timeDiff,logStepErr;
		errXtime.open("C:\\Users\\CZOLG\\Desktop\\lab11\\errorsXtime.txt", std::ios::out);
		logStepErr.open("C:\\Users\\CZOLG\\Desktop\\lab11\\errorChangeStepC.txt", std::ios::out);

		if (errXtime.bad() && logStepErr.bad())
			throw new exception("cant open file stream errXtime/timeDiff");
		else
		{
			errXtime << scientific << setprecision(16);
			timeDiff << scientific << setprecision(16);
		}


		double dt;
		double tNow;
		double maxErr = -INFINITY;
		int xNumberOfSteps;
		int tNumberOfSteps;
		double* temp;
		double analiticResult;
		double xNowErr;


		dt = (lambda * (h * h)) / parameterD;
		tNow = dt;

		xNumberOfSteps = (int)(xMax / h);
		double* xNow = createVector(xNumberOfSteps);
		double* xPrev = createVector(xNumberOfSteps);

		tNumberOfSteps = (int)(maxT / dt);

		for (int i = 0; i < xNumberOfSteps; i++)
			xPrev[i] = 1.0;


		for (int timeStep = 0 ; timeStep < tNumberOfSteps; timeStep++)
		{
			xNow[0] = 0.0;
			xNow[xNumberOfSteps - 1] = 1.0;

			for (int i = 1; i < xNumberOfSteps-1; i++)
			{

				xNow[i] = lambda * xPrev[i - 1] + (1.0 - (2.0 * lambda)) * xPrev[i] + lambda * xPrev[i+1];
			}

//save depending on timeStep
#pragma region filesave
			if (timeStep == tNumberOfSteps - 1)
			{
				timeDiff.open("C:\\Users\\CZOLG\\Desktop\\lab11\\timeMax.txt", std::ios::out);
			}
			else if (timeStep == (int)tNumberOfSteps * 0.05)
			{
				timeDiff.open("C:\\Users\\CZOLG\\Desktop\\lab11\\time005.txt", std::ios::out);
			}
			else if (timeStep == (int)tNumberOfSteps * 0.5)
			{
				timeDiff.open("C:\\Users\\CZOLG\\Desktop\\lab11\\time05.txt", std::ios::out);
			}
			else if (timeStep == (int)(tNumberOfSteps * 0.8))
			{
				timeDiff.open("C:\\Users\\CZOLG\\Desktop\\lab11\\time08.txt", std::ios::out);
			}

			double x = 0.0;
			for (int j = 0; j < xNumberOfSteps; j++,x += h)
			{
				timeDiff << setw(18) << " " << x << setw(18) << " " << xNow[j] << setw(18) << " " << analitic(x, tNow) << endl;
			}
			
			timeDiff.close();
#pragma endregion

			//error:
			double xd = 0.0;
			maxErr = abs(analitic(xd, tNow) - xNow[0]);
			analiticResult = 0.0;
			xNowErr = 0.0;

			for (int i = 0; i < xNumberOfSteps; i++)
			{
				xNowErr = abs(analitic(xd, tNow) - xNow[i]);
				if (xNowErr > maxErr)
					maxErr = xNowErr;

				xd += h;
			}

			//save
			errXtime << setw(18) << tNow << " " << setw(18) << maxErr << endl;


			temp = xPrev;
			xPrev = xNow;
			xNow = temp;

		
			tNow += dt;
		}
		errXtime.close();

		//error after done

		double xd = 0.0;
		tNow = tNow - dt;
		analiticResult = 0.0;
		xNowErr = 0.0;
		for (int i = 0; i < xNumberOfSteps; i++)
		{
			xNowErr = abs(analitic(xd, tNow) - xNow[i]);
			if (xNowErr > maxErr)
				maxErr = xNowErr;

			xd += h;
		}

		logStepErr << "h ="<< h << "  log(h) = " << scientific << log10(h) << "  log(err) = " << log10(maxErr);
		cout << "log(h) = " << scientific << log10(h) << "  log(err) = " << log10(maxErr);

		delete[] xPrev, xNow;

	}


	void directCrankNicolsonWithThomas()
	{
		fstream errXtime, timeDiff, logStepErr;
		errXtime.open("C:\\Users\\CZOLG\\Desktop\\lab11\\errorsXtimeCNT.txt", std::ios::out);
		logStepErr.open("C:\\Users\\CZOLG\\Desktop\\lab11\\errorChangeStepCN.txt", std::ios::out);

		if (errXtime.bad() && timeDiff.bad())
			throw new exception("cant open file stream errXtime");
		else
		{
			errXtime << scientific << setprecision(16);
			timeDiff << scientific << setprecision(16);
		}


		double dt;
		double timeStep;
		double maxErr = -INFINITY;
		int xNumberOfSteps;
		int tNumberOfSteps;

		//helping vars
		double analiticResult;
		double xNowErr;

		//helping variables
		double xd = 0.0;

		dt = (lambda * (h * h)) / parameterD;
		double tNow = dt;

		xNumberOfSteps = (int)(xMax / h);

		double* b = createVector(xNumberOfSteps);
		double* xFinal = createVector(xNumberOfSteps);
		double* thomasReshapedVector = createVector(xNumberOfSteps);

		//double** upper = createMatrix(xNumberOfSteps, xNumberOfSteps);
		//double** mid = createMatrix(xNumberOfSteps, xNumberOfSteps);
		//double** low = createMatrix(xNumberOfSteps, xNumberOfSteps);

		double** A = createMatrix(xNumberOfSteps,3);
		double** thomasReshapedMatrix = createMatrix(xNumberOfSteps,2);


		tNumberOfSteps = (int)(maxT / dt);


		for (int i = 0; i < xNumberOfSteps; i++)
			xFinal[i] = 1.0;

		b[0] = 0.0;
		b[xNumberOfSteps - 1] = 1.0;
		CN(A,xNumberOfSteps);

		thomasReadyMatrix(A, thomasReshapedMatrix,xNumberOfSteps);
		
		for (int timeStep = 0; timeStep < tNumberOfSteps; timeStep++)
		{

			for (int i = 1; i < xNumberOfSteps-1; i++)
				b[i] = -((lambda * 0.5) * xFinal[i - 1] + (1.0 - lambda) * xFinal[i] + (lambda * 0.5) * xFinal[i + 1]);

			thomasReadyVector(b, thomasReshapedVector, A, thomasReshapedMatrix,xNumberOfSteps);
			ThomasAlgo(thomasReshapedMatrix, thomasReshapedVector, xFinal,xNumberOfSteps);

//save depending on timeStep
#pragma region filesave
			if (timeStep == tNumberOfSteps - 1.0)
			{
				timeDiff.open("C:\\Users\\CZOLG\\Desktop\\lab11\\timeMaxCNT.txt", std::ios::out);
			}
			else if (timeStep == tNumberOfSteps * 0.05)
			{
				timeDiff.open("C:\\Users\\CZOLG\\Desktop\\lab11\\time005CNT.txt", std::ios::out);
			}
			else if (timeStep == tNumberOfSteps * 0.5)
			{
				timeDiff.open("C:\\Users\\CZOLG\\Desktop\\lab11\\time05CNT.txt", std::ios::out);
			}
			else if (timeStep == (int)(tNumberOfSteps * 0.8))
			{
				timeDiff.open("C:\\Users\\CZOLG\\Desktop\\lab11\\time08CNT.txt", std::ios::out);
			}

			xd = 0.0;
			for (int j = 0; j < xNumberOfSteps; j++, xd += h)
			{
				timeDiff << setw(18) << " " << xd << setw(18) << " " << xFinal[j] << setw(18) << " " << analitic(xd, tNow) << endl;
			}

			timeDiff.close();
#pragma endregion

			//error:
			xd = 0.0;
			maxErr = abs(analitic(xd, tNow) - xFinal[0]);

			analiticResult = 0.0;
			xNowErr = 0.0;
			for (int i = 0; i < xNumberOfSteps; i++)
			{
				xNowErr = abs(analitic(xd, tNow) - xFinal[i]);
				if (xNowErr > maxErr)
					maxErr = xNowErr;

				xd += h;
			}

			//save
			errXtime << setw(18) << tNow << " " << setw(18) << maxErr << endl;

			tNow = tNow + dt;
		}


		errXtime.close();
		tNow = tNow - dt;

		//max error after done
		xd = 0.0;
		analiticResult = 0.0;
		xNowErr = 0.0;
		for (int i = 0; i < xNumberOfSteps; i++)
		{
			xNowErr = abs(analitic(xd, tNow) - xFinal[i]);
			if (xNowErr > maxErr)
				maxErr = xNowErr;

			xd += h;
		}

		logStepErr << "h =" << h << "  log(h) = " << scientific << log10(h) << "  log(err) = " << log10(maxErr);
		cout << "log(h) = " << scientific << log10(h) << "  log(err) = " << log10(maxErr);


		for (int i = 0; i < xNumberOfSteps; i++)
		{
			delete[]  A[i], thomasReshapedMatrix[i];
		}
		delete[]  A, thomasReshapedVector, b, xFinal;

	}

#pragma region CrankNicolJacMeth
	void directCrankNicolsonWithJacob()
	{
		fstream errXtime, timeDiff;
		errXtime.open("C:\\Users\\CZOLG\\Desktop\\lab11\\errorsXtimeCNJ.txt", std::ios::out);

		if (errXtime.bad() && timeDiff.bad())
			throw new exception("cant open file stream errXtime");
		else
		{
			errXtime << scientific << setprecision(16);
			timeDiff << scientific << setprecision(16);
		}


		double dt;
		double t;
		double maxErr = -INFINITY;
		int xNumberOfSteps;
		int tNumberOfSteps;

		//helping variables
		double x = 0.0;

		dt = (lambda * (h * h)) / parameterD;
		t = dt;

		xNumberOfSteps = (int)(xMax / h);
		tNumberOfSteps = (int)(maxT / dt);

		double* b = createVector(xNumberOfSteps);
		double* thomasReshapedVector = createVector(xNumberOfSteps);

		double** M = createMatrix(xNumberOfSteps, xNumberOfSteps);
		double** D = createMatrix(xNumberOfSteps, xNumberOfSteps);
		double** LU = createMatrix(xNumberOfSteps, xNumberOfSteps);
		double** squareA = createMatrix(xNumberOfSteps, xNumberOfSteps);

		double* C = createVector(xNumberOfSteps);
		double* u = createVector(xNumberOfSteps);

		double** A = createMatrix(xNumberOfSteps, 3);
		double** thomasReshapedMatrix = createMatrix(xNumberOfSteps, 2);





		for (int i = 0; i < xNumberOfSteps; i++)
			u[i] = 1.0;

		CN(A, xNumberOfSteps);

		for (int i = 0; i < xNumberOfSteps; i++)
		{
			for (int j = -1; j <= 1; j++)
			{
				if (j ==0)
				{
					D[i][i] = -1 / A[i][j + 1];
				}
				else
				{
					if (i+j>= 0 && i+j<xNumberOfSteps)
					{
						LU[i][i + j] = A[i][j + 1];
					}
				}
				squareA[i][i + j] = A[i][j + 1];
			}
		}
		M = multiplyMatrixiesSameShapeReturn(D, LU, xNumberOfSteps);

		for (int k = 0; k < tNumberOfSteps; k++)
		{
			b[0] = 0.0;
			b[xNumberOfSteps - 1] = 1.0;
			for (int i = 1; i < xNumberOfSteps - 1; i++)
				b[i] = -((lambda * 0.5) * u[i - 1] + (1.0 - lambda) * u[i] + (lambda * 0.5) * u[i + 1]);

			b[xNumberOfSteps - 1] = 1.0;
			C = multiVecMatReturn(b, D, xNumberOfSteps);

			for (int i = 0; i < xNumberOfSteps; i++)
				C[i] = C[i];

			JacobiMethod(squareA,b,M,C,u,xNumberOfSteps,100);

#pragma region filesave
			if (k == tNumberOfSteps - 1.0)
			{
				timeDiff.open("C:\\Users\\CZOLG\\Desktop\\lab11\\timeMaxCNJ.txt", std::ios::out);
			}
			else if (k == tNumberOfSteps *0.05)
			{
				timeDiff.open("C:\\Users\\CZOLG\\Desktop\\lab11\\time-20CNJ.txt", std::ios::out);
			}
			else if (k == tNumberOfSteps *0.5)
			{
				timeDiff.open("C:\\Users\\CZOLG\\Desktop\\lab11\\timeDiff-2CNJ.txt", std::ios::out);
			}
			else if (k == (int)(tNumberOfSteps * 0.8))
			{
				timeDiff.open("C:\\Users\\CZOLG\\Desktop\\lab11\\timeDiff-3/4CNJ.txt", std::ios::out);
			}

			x = 0.0;
			for (int j = 0; j < xNumberOfSteps; j++, x += h)
			{
				timeDiff << setw(18) << " " << x << setw(18) << " " << u[j] << setw(18) << " " << analitic(x, t) << endl;
			}

			timeDiff.close();
#pragma endregion

			//error:
			x = 0.0;
			maxErr = abs(analitic(x, t) - u[0]);

			double analiticResult = 0.0;
			double xNowErr = 0.0;

			for (int i = 0; i < xNumberOfSteps; i++)
			{
				analiticResult = analitic(x, t);
				xNowErr = abs(analitic(x, t) - u[i]);

				if (xNowErr > maxErr && i < xNumberOfSteps)
					maxErr = xNowErr;

				x += h;
			}

			//save
			errXtime << setw(18) << t << " " << setw(18) << maxErr << endl;

			t = t + dt;
		}


		errXtime.close();

		t = t - dt;

		//error after done
		double analiticResult;
		double xNowErr;

		x = 0.0;
		for (int i = 0; i < xNumberOfSteps; i++)
		{
			analiticResult = analitic(x, t);
			xNowErr = abs(analitic(x, t) - u[i]);

			if (xNowErr > maxErr && i < xNumberOfSteps)
				maxErr = xNowErr;

			x += h;
		}

		cout << "log(h) = " << scientific << log10(h) << "  log(err) = " << log10(maxErr);
	}

#pragma endregion

	void test()
	{
		classicDirectMethod();
		directCrankNicolsonWithThomas();
		//directCrankNicolsonWithJacob();
	}
}