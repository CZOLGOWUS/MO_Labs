#include "lab4.hpp"


using namespace std;

namespace lab4
{



	void test()
	{

		std::function<double(double, double, double)> Jrev[3][3];
		std::function<double(double, double, double)> fXnArr[3];
		std::function<double(double, double, double)> lambda2[3];


		Jrev[0][0] = ([](double x, double y, double z) -> double {return 0.0; });
		Jrev[0][1] = ([](double x, double y, double z) -> double {return 1.0 / (2.0 * x + 4.0 * x * y); });
		Jrev[0][2] = ([](double x, double y, double z) -> double {return y / (x + 2.0 * x * y); });
		Jrev[1][0] = ([](double x, double y, double z) -> double {return 0.0; });
		Jrev[1][1] = ([](double x, double y, double z) -> double {return 1.0 / (2.0 * y + 1.0); });
		Jrev[1][2] = ([](double x, double y, double z) -> double {return -1.0 / (2.0 * y + 1.0); });
		Jrev[2][0] = ([](double x, double y, double z) -> double {return 1.0 / (2.0 * z); });
		Jrev[2][1] = ([](double x, double y, double z) -> double {return -1.0 / (2.0 * z); });
		Jrev[2][2] = ([](double x, double y, double z) -> double {return 0.0; });

		fXnArr[0] = ([](double x, double y, double z) -> double {return f1(x, y, z); });
		fXnArr[1] = ([](double x, double y, double z) -> double {return f2(x, y, z); });
		fXnArr[2] = ([](double x, double y, double z) -> double {return f3(x, y, z); });

		 //all lambda functions pre calculated
		lambda2[0] = ([](double x, double y, double z) -> double {return (x * x + 3.0 * (y * y) + 2.0 * (x * x) * y - 1.0) / (2.0 * x + 4.0 * x * y); });
		lambda2[1] = ([](double x, double y, double z) -> double {return (y * y - y - 1.0) / (2.0 * y + 1.0); });
		lambda2[2] = ([](double x, double y, double z) -> double {return (z*z - 1.0) / (2.0 * z); });
		


		double xn[3] = { 0.8,0.6,2.0 };
		double Xnext[3] = { 0.0 };

		double JrevValues[3][3] = { 0.0 };
		double fXnArrValues[3] = { 0.0 };

		double lambdaV[3] = { 0.0 };

		double reziduum[3] = { 0.0 };
		double error[3] = { 0.0 };

		double tolx = 0.0001;
		double tolf = 0.0001;
		int maxN = 100;
		int n = 0;


		cout << "enter starting x,y,z:\n";
		cin >> xn[0] >> xn[1] >> xn[2];
		cout << "enter max iterations:\n";
		cin >> maxN;
		cout << "enter tolx:\n";
		cin >> tolx;
		cout << "enter tolf:\n";
		cin >> tolf;
		cout << endl;

		while (true)
		{
			for (int i = 0; i < 3; i++)
			{
				fXnArrValues[i] = fXnArr[i](xn[0], xn[1], xn[2]);

				for (int j = 0; j < 3; j++)
				{
					JrevValues[i][j] = Jrev[i][j](xn[0], xn[1], xn[2]);
				}
			}


			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
				{
					lambdaV[i] += JrevValues[i][j] * fXnArrValues[j];
				}

			/*
			for (int i = 0; i < 3 ; i++)
				lambdaV[i] = lambda2[i](xn[0], xn[1], xn[2]);
			*/
			
			for (int i = 0; i < 3 ;i++)
			{
				Xnext[i] = xn[i] - lambdaV[i];
			}


			reziduum[0] = abs(f1(Xnext[0], Xnext[1], Xnext[2]));
			reziduum[1] = abs(f2(Xnext[0], Xnext[1], Xnext[2]));
			reziduum[2] = abs(f3(Xnext[0], Xnext[1], Xnext[2]));
			
			error[0] = abs(Xnext[0] - xn[0]);
			error[1] = abs(Xnext[1] - xn[1]);
			error[2] = abs(Xnext[2] - xn[2]);


			cout << setw(4) << "Xn = " << setw(15) << Xnext[0] << "  lambda for x:" << setw(10) << lambdaV[0] 
				 << "  reziduum = "<< setw(10) << reziduum[0] << "  error = " << setw(10) << error[0] << endl;

			cout << setw(4) << "Yn = " << setw(15) << Xnext[1] << "  lambda for y:" << setw(10) << lambdaV[1] 
				<< "  reziduum = " << setw(10) << reziduum[1] << "  error = " << setw(10) << error[1] << endl;

			cout << setw(4) << "Zn = " << setw(15) << Xnext[2] << "  lambda for z:" << setw(10) << lambdaV[2] 
				<< "  reziduum = " << setw(10) << reziduum[2] << "  error = " << setw(10) << error[2] << endl;
			cout << endl;


			if (n >= maxN)
			{
				cout << "iteration over extended\n";
				break;
			}
			else if (error[0] <= tolx && error[1] <= tolx && error[2] <= tolx)
			{
				cout << "result is in given error tolerance\n";
				break;
			}
			else if (reziduum[0] < tolf && reziduum[1] < tolf && reziduum[2] < tolf)
			{
				cout << "result is in given reziduum tolerance\n";
				break;
			}
			else if (reziduum == 0)
			{
				cout << "perfect approximation has been achived!\n";
				break;
			}

			//clear before next iteration
			for (int i = 0; i < 3; i++)
			{
				xn[i] = Xnext[i];
				Xnext[i] = 0.0;
				lambdaV[i] = 0.0;
			}
			n++;
		}

		cout << "approximated values of [x,y,z] are : x = " << xn[0] << " y = " << xn[1] << " z = " << xn[2] << endl;

	}

}