#include <iostream>
#include <fstream>
#include <cstdio>
#include <iomanip>
#include <math.h>
#include <string>
#include <sstream>

using namespace std;

namespace lab1 
{
#pragma region zaj1zad1

    template <typename T>
    void zad1Temp()
    {
        T epsilon = (T)0;
        T x = (T)1;
        int count = 0;

        for (int i = 1; x < (T)2; i++)
        {
            x += (T)pow(2, -i);

            if (x >= 2)
                break;

            epsilon = (T)pow(2, -i);
            count++;

            //cout << "x=" << setprecision(30) << x << " epsilon=" << epsilon << " mantise bits=" << count << endl;
        }

        // liczba bitów zgadza siê z liczb¹ przejsc fora z prostego powodu ustawiania
        // kolejnych bitów w mantysie i gdy wszsytkie bity mantysy s¹ ustawione na 1
        //  nastêpuje "zaokr¹glenie"(overflow) w góre a jako mniejsze liczby niz najmniejszy bit mantysy
        //  nie mog¹ byæ przedstawione w danym typie to ta liczba siê niepowieksza

        cout << "x=" << setprecision(30) << x << " epsilon=" << epsilon << " mantise bits=" << count << endl;
    }

    void test()
    {
        zad1Temp<float>();
        zad1Temp<double>();
    }

#pragma endregion
}
namespace lab2
{

#pragma region zaj1zad2

    double exponential(double x,int n)
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

        //1-exp(x):
        int n = 100;
        double ex = 0.0;

        double facN,powX;

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

        //cout << ex << endl;
        return -ex / x;

    }


#pragma endregion
}

int main()
{

#pragma region zaj1zad1
    //lab1::test();
#pragma endregion

#pragma region zaj1zad2

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


    while (getline(xData, xString) && getline(fData,fString))
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
    dataNormalFun.close();
    dataTaylorFun.close();


#pragma endregion



}