#include <iostream>
#include <iomanip>


namespace lab1
{
    using namespace std;
#pragma region lab1

    template <typename T>
    static void zad1Temp()
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

    static void test()
    {
        zad1Temp<float>();
        zad1Temp<double>();
    }

#pragma endregion
}