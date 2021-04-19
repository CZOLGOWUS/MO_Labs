#include "lab3.hpp"

namespace lab3
{

using namespace std;

#pragma region Picard

    double picard(double (*fun)(double),double (*phi)(double),double xStart,int maxN , double tolx, double tolf)
    {
        double reziduum = abs(fun(xStart));
        double xNew = NULL;
        double xOld = xStart;
        double error = NULL;
        int n = 1;


        cout << setw(20) << "new X"
            << setw(20) << "reziduum"
            << setw(20) << "error estimator" << "\n";

        while (1)
        {
            xNew = phi(xOld);
            reziduum = abs(fun(xNew));
            error = abs(xOld - xNew);

            cout << setw(20) << setprecision(10) << xNew
                << setw(20) << setprecision(10) << reziduum
                << setw(20) << setprecision(10) << error << endl;

            if (n >= maxN)
            {
                cout << "iteration over extended";
                break;
            }
            else if (error <= tolx)
            {
                cout << "result is in given error tolerance";
                break;
            }
            else if (reziduum < tolf)
            {
                cout << "result is in given reziduum tolerance";
                break;
            }
            else if (reziduum == 0)
            {
                cout << "perfect approximation has been achived!\n";
                break;
            }

            xOld = xNew;
            n++;
        }

        return xNew;
    }
    void picardTest()
    {
        cout << "\naproximation for picards method for sin(x/4) - x = 0 is: " << picard(
            [](double x) -> double {return sin(x / 4) - x; },
            [](double x) -> double {return sin(x / 4); },
            5.0,
            1000,
            0.0001,
            0.0001) << endl << endl;

        cout << "\naproximation for picards method for tan(2*x) - 1 - x = 0 is: " << picard(
            [](double x) -> double {return tan(2 * x) - 1 - x; },
            [](double x) -> double {return tan(2 * x) - 1; },
            1.5,
            100,
            0.0001,
            0.0001) << endl<<endl<<endl;
    }
#pragma endregion

#pragma region bisection

    double bisection(double (*fun)(double),double jumpAmountOfAB,int numberOfJumps, int maxN, double tolx, double tolf)
    {
        double a = -2, b = 1;
        double aNew = NULL, bNew = NULL;
        double x = NULL;
        int n = 0;
        double reziduum = NULL;
        double error = NULL;

        while (signbit(fun(a)) == signbit(fun(b)) || n <= numberOfJumps)
        {
            n++;
            a -= jumpAmountOfAB;
            b += jumpAmountOfAB;
        }
        n = 1;

        cout << "[a,b] = [" << a << ","<< b << "]\n";

        cout << setw(20) << setprecision(10) << "new X"
            << setw(20) << setprecision(10) << "reziduum"
            << setw(20) << setprecision(10) << "error" << endl;


        while (1)
        {
            x = (a + b) / 2.0;

            if (signbit(fun(x)) == signbit(fun(a)))
                a = x;
            else
                b = x;

            error = abs((b - a) / 2);

            reziduum = abs(fun(x));


            cout << setw(20) << setprecision(10) << x
                << setw(20) << setprecision(10) << reziduum
                << setw(20) << setprecision(10) << error << endl;



            if (reziduum < tolf)
            {
                cout << "result is in given reziduum tolerance";
                break;
            }
            else if (error <= tolx)
            {
                cout << "result is in given error tolerance";
                break;
            }
            else if (reziduum == 0)
            {
                cout << "perfect approximation has been achived!\n";
                break;
            }
            else if (n >= maxN)
            {
                cout << "iteration over extended";
                break;
            }

            n++;
        }

        return x;

    }
    void bisectionTest()
    {
        cout << "\naproximation for bisection method for sin(x/4) - x = 0 is: " << bisection(
            [](double y) -> double {return sin(y / 4) - y; },
            0.01,
            1000,
            1000,
            0.000001,
            0.000001) << endl << endl;

        cout << "\naproximation for bisection method for tan(2*x) - 1  - x = 0 is: " << bisection(
            [](double x) -> double {return tan(2 * x) - 1 - x; },
            0.01,
            1000,
            1000,
            0.000001,
            0.000001) << endl << endl;
    }

#pragma endregion

#pragma region Newton
    double newton(double (*fun)(double),double (*funDev)(double),double xStart, int maxN, double tolx, double tolf)
    {
        double xPrev = xStart;
        double x = xStart;
        double xPrevOld = NULL;
        int n = 1;
        double reziduum = NULL;
        double error = NULL;


        cout << setw(20) << setprecision(10) << "new X"
            << setw(20) << setprecision(10) << "reziduum"
            << setw(20) << setprecision(10) << "error" << endl;


        while (1)
        {
            if (funDev(xPrev) == 0)
            {
                cout << "deveritive of function in new X is 0\n";
                break;
            }

            x = xPrev - (fun(xPrev) / funDev(xPrev));
            error = abs(x - xPrev);
            reziduum = abs(fun(x));


            cout << setw(20) << setprecision(10) << x
                << setw(20) << setprecision(10) << reziduum
                << setw(20) << setprecision(10) << error << endl;



            if (reziduum < tolf)
            {
                cout << "result is in given reziduum tolerance";
                break;
            }
            else if (error <= tolx)
            {
                cout << "result is in given error tolerance";
                break;
            }
            else if (reziduum == 0)
            {
                cout << "perfect approximation has been achived!\n";
                break;
            }
            else if (n >= maxN)
            {
                cout << "iteration over extended";
                break;
            }
            else if (xPrevOld == x)
            {
                cout << "iteration are in a loop";
                break;
            }

            xPrevOld = xPrev;
            xPrev = x;
            n++;

        }

        return x;

    }

    void newtonTest()
    {
        cout << "\naproximation for Newton method for sin(x/4)  - x = 0 is: " << newton(
            [](double x) -> double {return sin(x / 4) - x; },
            [](double x) -> double {return cos(x / 4) / 4 - 1.0; },
            5.2516,
            1000,
            0.000001,
            0.0001) << endl << endl;


        cout << "\naproximation for Newton method for tan(2*x) - 1  - x = 0 is: " << newton(
            [](double x) -> double {return tan(2.0 * x) - 1 - x; },
            [](double x) -> double {return 2.0 * 1 / pow(cos(2.0 * x), 2) - 1.0; },
            0.2516,
            1000,
            0.0001,
            0.0001) << endl << endl;
    }

#pragma endregion

#pragma region secant

    double secant(double (*fun)(double), double xStart1,double xStart2, int maxN, double tolx, double tolf)
    {
        double xPrevOld = xStart2;
        double xPrev = xStart1;
        double x = NULL;

        int n = 1;
        double reziduum = NULL;
        double error = NULL;


        cout << setw(20) << setprecision(10) << "new X"
            << setw(20) << setprecision(10) << "reziduum"
            << setw(20) << setprecision(10) << "error" << endl;


        while (1)
        {
            if ((xPrev - xPrevOld) == 0 || (fun(xPrev) - fun(xPrevOld) == 0))
            {
                cout << "cant divide by zero\n";
                break;
            }

            x = xPrev - fun(xPrev) / ((fun(xPrev) - fun(xPrevOld)) / (xPrev - xPrevOld));
            reziduum = abs(fun(x));
            error = abs(x - xPrev);

            cout << setw(20) << setprecision(10) << x
                << setw(20) << setprecision(10) << reziduum
                << setw(20) << setprecision(10) << error << endl;


            if (reziduum < tolf)
            {
                cout << "result is in given reziduum tolerance";
                break;
            }
            else if (error <= tolx)
            {
                cout << "result is in given error tolerance";
                break;
            }
            else if (reziduum == 0)
            {
                cout << "perfect approximation has been achived!\n";
                break;
            }
            else if (n >= maxN)
            {
                cout << "iteration over extended";
                break;
            }


            xPrevOld = xPrev;
            xPrev = x;
            n++;
        }
        return x;
    }

    void secantTest()
    {
        cout << "\naproximation for secant method for sin(x/4)  - x = 0 is: " << secant(
            [](double x) -> double {return sin(x / 4) - x; },
            5.2516,
            1.24,
            1000,
            0.0001,
            0.0001) << endl << endl;


        cout << "\naproximation for secant method for tan(2*x) - 1  - x = 0 is: " << secant(
            [](double x) -> double {return tan(2.0 * x) - 1 - x; },
            0.2516,
            1.24,
            1000,
            0.0001,
            0.0001) << endl << endl;
    }

#pragma endregion

    void test()
    {
        char choice = '0';

        while (choice != 'e')
        {
            
            cout << "tolf = 0.0001 | tolx = 0.0001\n";
            cout << "menu:\n";
            cout << "1. Picard method\n";
            cout << "2. Bisection method\n";
            cout << "3. Newton method\n";
            cout << "4. Secant method\n";
            cout << "'e' to exit\n";
            cout << "choose from 1 to 4:";
            cin >> choice;

            switch (choice)
            {
                case '1':
                {
                    cout << "Picard method:\n";
                    picardTest();
                    break;
                }
                case '2':
                {
                    cout << "Bisection Method:\n";
                    bisectionTest();
                    break;
                }
                case '3':
                {
                    cout << "Netwon Method test:\n";
                    newtonTest();
                    break;
                }
                case '4':
                {
                    cout << "Secant Method test:\n";
                    secantTest();
                    break;
                }
                case 'e':
                {
                    break;
                }
            }
        }

    }

}

