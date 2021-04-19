#include "lab6.hpp"


namespace lab6
{
	
	vector<double>& mThomas(const vector<double>& L, const vector<double>& D, const vector<double>& U)
	{
		vector<double>* DThom = new vector<double>;
		vector<double>& dRef = *DThom;

		dRef = D;
		
		for (int i = 1, j = 0, k = 1; i < (int)D.size(); i++)
		{
			dRef[i] -= L[j++] * U[k++ - 1] / dRef[i - 1];
		}

		return dRef;

	}

	vector<double>& vThomas(const vector<double>& L, const vector<double>& D,const vector<double> B)
	{
		vector<double>* v = new vector<double>;
		vector<double>& vRef = *v;

		vRef = B;

		for (int i = 1, j = 0; i < (int)D.size(); i++)
		{
			vRef[i] -= L[j++] * B[i - 1] / D[i - 1];
		}

		return vRef;
	}


	vector<double>& getVectorX(const vector<double>& U, const vector<double>& D, const vector<double> B)
	{

		vector<double>* x = new vector<double>;
		vector<double>& xRef = *x;

		for (int i = 0; i < (int)B.size(); i++)
			xRef.push_back(0.0);

		xRef[D.size() - 1] = B[D.size() - 1] / D[D.size() - 1];

		for (int i = (int)xRef.size() - 2; i >= 0; i--)
		{
			xRef[i] = (B[i] - (U[i] * xRef[i + 1])) / D[i];
		}


		return xRef;
	}



	void test()
	{

		vector<vector<double>> A;
		vector<double> U;
		vector<double> D;
		vector<double> L;

		vector<double> B;

		vector<double> X;


		A =
		{
			{0.5 , 0.25 , 1.0 / 6.0 , 1.0 / 8.0 , 1.0 / 10.0},
			{10.0 , 20.0 , 30.0 , 30.0 , 20.0 , 10.0 } ,
			{1.0/3.0 , 1.0/5.0 , 1.0/7.0 , 1.0/9.0 , 1.0/11.0}
		};
		
		U = { 0.5 , 0.25 , 1.0 / 6.0 , 1.0 / 8.0 , 1.0 / 10.0 };
		D = { 10.0 , 20.0 , 30.0 , 30.0 , 20.0 , 10.0 };
		L = { 1.0 / 3.0 , 1.0 / 5.0 , 1.0 / 7.0 , 1.0 / 9.0 , 1.0 / 11.0 };
		
		B = { 31.0, 165.0 / 4.0 , 917.0 / 30.0 , 851.0 / 28.0 , 3637.0 / 90.0 , 332.0 / 11.0 };
		
		X = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };


		vector<double> DThom = mThomas(L, D, U);
		vector<double> VectorThom = vThomas(L, DThom, B);
		vector<double> VectorX = getVectorX(U, DThom, VectorThom);


		cout << "Diagonal Vector after Thomas function D:\n";
		for (int i = 0; i < (int)DThom.size(); i++)
		{
			cout << DThom[i] << "  ";
		}
		cout << endl;
		cout << endl;

		cout << "B vector after Thomas function B:\n";
		for (int i = 0; i < (int)VectorThom.size(); i++)
		{
			cout << VectorThom[i] << "  ";
		}
		cout << endl;
		cout << endl;


		cout << "Solution Vector X:\n";
		for (int i = 0; i < (int)VectorX.size(); i++)
		{
			cout << VectorX[i] << "  ";
		}
		cout << endl;

	}
}