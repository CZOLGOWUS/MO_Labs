#include "lab2.hpp"


namespace lab2
{
    using namespace std;
#pragma region lab2

    double exponential(double x, int n)
    {
        double sum = 1.0f; // initialize sum of series  

        for (int i = n - 1; i > 0; --i)
            sum = 1 + x * sum / i;

        return sum;
    }

    double eToPow(double x, int n)
    {
        static double p = 1, f = 1;
        double r;

        if (n == 0)
            return 1;

        r = eToPow(x, n - 1);

        p = p * x;
        f = f * n;

        return (r + p / f);
    }

    double NormalFunction(double x)
    {
        return (1 - exp(-x)) / x;
    }

    double TaylorFunction(double x)
    {

        if (fabs(x) < 0.5)
        {
            //1-exp(x):
            int n = 100;
            double ex = 0.0;

            double facN, powX;

            for (int i = 1; i < n; i++)
            {
                facN = 1.0;
                for (int j = 1; j <= i; j++)
                {
                    facN *= j;
                }
                powX = pow(-x, i);
                ex += powX / facN;
            }

            return -ex / x;
        }

        return (1 - exp(-x)) / x;

    }

    void test()
    {
        stringstream sstream;
        string xString;
        string fString;

        double x = 0.0;
        double fx = 0.0;

        ifstream xData("E:\\studia\\semestr 4\\Metody obliczeniowe\\laby\\plots\\lab2\\XdataRef.txt");
        ifstream fData("E:\\studia\\semestr 4\\Metody obliczeniowe\\laby\\plots\\lab2\\fdataRef.txt");
        //ifstream logData("E:\\studia\\semestr 4\\Metody obliczeniowe\\laby\\plots\\lab2\\logdataRef.txt");
        ofstream dataNormalFun("E:\\studia\\semestr 4\\Metody obliczeniowe\\laby\\plots\\lab2\\normalFun.txt");
        ofstream dataTaylorFun("E:\\studia\\semestr 4\\Metody obliczeniowe\\laby\\plots\\lab2\\TaylorFun.txt");


        while (getline(xData, xString) && getline(fData, fString))
        {
            sstream << xString;
            sstream >> x;
            sstream.clear();

            sstream << fString;
            sstream >> fx;
            sstream.clear();

            dataNormalFun << setprecision(20) << scientific << log10(x) << " " << setprecision(20) << scientific << log10(fabs(lab2::NormalFunction(x) - fx)) << "\n";
            dataTaylorFun << setprecision(20) << scientific << log10(x) << " " << setprecision(20) << scientific << log10(fabs(lab2::TaylorFunction(x) - fx)) << "\n";
        }


        xData.close();
        fData.close();
        dataNormalFun.close();
        dataTaylorFun.close();
    }

#pragma endregion
}