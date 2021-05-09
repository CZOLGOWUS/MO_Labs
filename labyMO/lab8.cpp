#include "lab8.hpp"

namespace lab8
{
	using namespace std;

#define PI (T)3.14159265358979323846
#define fun(x) (T)sin(x)
#define derivativeFun(x) (T)cos(x)
#define numberOfCols 8
#define numberOfIterrations 200

//defining 3 main points in graph
#define start (T)0.0
#define mid (T)PI/4.0
#define end (T)PI/2.0

//two pointed differences
	template <typename T> 
	T twoProgressive(T x,T h)
	{
		return (fun(x + h) - fun(x)) / (h);
	}

	template <typename T>
	T twoBackward(T x,T h)
	{
		return (fun(x) - fun(x - h)) / (h);
	}

	template <typename T>
	T twoCentral(T x,T h)
	{
		return (fun(x + h) - fun(x - h)) / (2 * h);
	}

//three pointed differences of h-net at the beggining and the end of the h-net
	template <typename T>
	T threeProgressive(T x, T h)
	{
		return ( -((T)3.f / (T)2.f) * fun(x) + (T)2.f * fun(x + h) - ((T)1.f / (T)2.f) * fun(x + h + h) ) / (h);
	}

	template <typename T>
	T ThreeBackward(T x, T h)
	{
		return ( ((T)1.f / (T)2.f) * fun(x-h-h) - (T)2.f * fun(x - h) + ((T)3.f / (T)2.f) * fun(x)) / (h);
	}

	template <typename T>
	T** calculateError(double toly,double tolx)
	{
		cout.precision(12);
		cout.scientific;

		T h = (T)0.2;
		T** error = new T * [numberOfIterrations];
		for (int i = 0; i < numberOfIterrations; i++)
			error[i] = new T[numberOfCols];

		for (int i = 0; i < numberOfIterrations; i++)
		{
			error[i][0] = abs(twoProgressive<T>(start, h) - derivativeFun(start));
			error[i][1] = abs(twoCentral<T>(mid, h) - derivativeFun(mid));
			error[i][2] = abs(twoBackward<T>(end, h) - derivativeFun(end));

			error[i][3] = abs(threeProgressive<T>(start, h) - derivativeFun(start));
			error[i][4] = abs(threeProgressive<T>(mid, h) - derivativeFun(mid));
						  
			error[i][5] = abs(ThreeBackward<T>(mid, h) - derivativeFun(mid));
			error[i][6] = abs(ThreeBackward<T>(end, h) - derivativeFun(end));

			error[i][7] = h;
			

			cout <<setw(15) << error[i][7] << " ";
			for (int j = 0; j < numberOfCols - 1; j++) 
				cout <<setw(15)<< error[i][j] << " ";
			 cout << endl;


			 h = h * (T)0.7;

		}

		cout << "\npresision levels: \n";
		std::cout <<setw(35)<< "2 point progresive start   : " << (log10(error[1][0]) - log10(error[0][0])) / (log10(error[1][7]) - log10(error[0][7])) << std::endl;
		std::cout <<setw(35)<< "2 point progresive mid   : " << (log10(error[1][1]) - log10(error[0][1])) / (log10(error[1][7]) - log10(error[0][7])) << std::endl;
		std::cout <<setw(35)<< "2 point progresive end   : " << (log10(error[1][2]) - log10(error[0][2])) / (log10(error[1][7]) - log10(error[0][7])) << std::endl;
		std::cout <<setw(35)<< "3 point progressive start   : " << (log10(error[1][3]) - log10(error[0][3])) / (log10(error[1][7]) - log10(error[0][7])) << std::endl;
		std::cout <<setw(35)<< "3 point progressive mid   : " << (log10(error[1][4]) - log10(error[0][4])) / (log10(error[1][7]) - log10(error[0][7])) << std::endl;
		std::cout <<setw(35)<< "3 point backwards mid   : " << (log10(error[1][5]) - log10(error[0][5])) / (log10(error[1][7]) - log10(error[0][7])) << std::endl;
		std::cout <<setw(35)<< "3 point backwards end   : " << (log10(error[1][6]) - log10(error[0][6])) / (log10(error[1][7]) - log10(error[0][7])) << std::endl;

		return error;
	}

	template <typename T> bool saveTable(T** arr,string path)
	{
		ofstream file;
		file.open(path);

		file.scientific;
		file.precision(12);


		if (file.good())
		{
			for (int i = 0; i < numberOfIterrations;i++)
			{
				file << setw(15) << log10(arr[i][7]) << " ";
				for (int j = 0; j < numberOfCols - 1 ; j++)
					file << setw(15)<< scientific << log10(arr[i][j]) << " ";
				file << endl;
			}
		}
		else
		{
			return false;
		}

		file.close();
		return true;
	}

	void test()
	{

		cout.precision(12);
		cout.scientific;

		float** derivitiveF = calculateError<float>(0.0001,0.0001);
		saveTable(derivitiveF, "C:\\Users\\CZOLG\\Desktop\\testFloats.txt");

		double** derivitiveD = calculateError<double>(0.0001, 0.0001);
		saveTable(derivitiveD, "C:\\Users\\CZOLG\\Desktop\\testDoubles.txt");

	}

}