#include "lab5.hpp"


namespace lab5
{

	int getIndexRowOfMaxValueInCollumn(vector<vector<float>> A, int col)
	{
		float max = -INFINITY;
		int index = 0;
		for (unsigned int i = col; i < A.size(); i++)
		{
			if (fabs(A[i][col]) > max)
			{
				max = fabs(A[i][col]);
				index = i;
			}
		}

		return index;
	}

	MatrixesLU& getLU(vector<vector<float>> A)
	{
		MatrixesLU* result = new MatrixesLU;

		result->L = 
		{
			{1,0,0,0},
			{0,1,0,0},
			{0,0,1,0},
			{0,0,0,1}
		};

		result->U = A;
		result->changes = { 0,1,2,3 };

		float ratio = 0;

		for (unsigned int i = 0; i < A[0].size(); i++)
		{
			if (fabs(result->U[i][i]) < 0.0001f)
			{
				int swapRow = getIndexRowOfMaxValueInCollumn(result->U, i);

				result->U[i].swap(result->U[swapRow]);
				
				int temp = result->changes[i];
				result->changes[i] = result->changes[swapRow];
				result->changes[swapRow] = temp;
			}
			for (unsigned int y = 1+i; y < A.size(); y++)
			{
				ratio = result->U[y][i] / result->U[i][i];

				result->L[A.size()-y+i][i] = ratio;

				for (unsigned int x = i; x < A[0].size(); x++)
					result->U[y][x] -= result->U[i][x] * ratio;
			}
		}
#pragma region testing

		/*
		
		for (unsigned int i = 0; i < A.size(); i++)
		{
			for (unsigned int j = 0; j < A[0].size(); j++)
				cout << setw(10) << result->U[i][j] << "  ";
			cout << endl;
		}

		cout << endl; 
		cout << endl;

		for (unsigned int i = 0; i < A.size(); i++)
		{
			for (unsigned int j = 0; j < A[0].size(); j++)
				cout << setw(10) << result->L[i][j] << "  ";
			cout << endl;
		}
		
		cout << endl;
		cout << endl;

		vector<vector<float>> multiTest = 
		{
			{ 0,0,0,0 },
			{ 0,0,0,0 },
			{ 0,0,0,0 },
			{ 0,0,0,0 }
		};

		
		cout << "checking with multiplication:\n";
		for (unsigned int i = 0; i < A.size(); i++)
		{
			for (unsigned int j = 0; j < A[0].size(); j++)
			{
				for (int k = 0; k < A.size(); k++)
					multiTest[i][j] += result->L[i][k] * result->U[k][j];
			}
		}

		for (unsigned int i = 0; i < A.size(); i++)
		{
			for (unsigned int j = 0; j < A[0].size(); j++)
				cout << setw(10) << multiTest[i][j] << "  ";
			cout << endl;
		}
		cout << endl; cout << endl;
		for (unsigned int i = 0; i < result->changes.size(); i++)
		{
			cout << result->changes[i] << "  ";
		}

		cout << endl; cout << endl;
		
		*/
#pragma endregion

		return *result;
	}


	vector<float>& getB(const vector<float>& b ,const MatrixesLU& LU)
	{
		vector<float>* result = new vector<float>;
		(*result) = b;

		vector<int> changed;

		for (unsigned int i = 0; i < b.size(); i++)
		{
			if (LU.changes[i] != i && 
				find(changed.begin(),changed.end(),i) == changed.end() &&
				find(changed.begin(), changed.end(), LU.changes[i]) == changed.end()
				)
			{
				changed.push_back(i);
				changed.push_back(LU.changes[i]);

				float temp = (*result)[i];
				(*result)[i] = (*result)[LU.changes[i]];
				(*result)[LU.changes[i]] = temp;
			}
		}

#pragma region test
		/*
		for (unsigned int i = 0; i < b.size(); i++)
			cout << (*result)[i] << endl;
		cout << endl;
		*/
#pragma endregion

		return (*result);
	}

	vector<float>& solveX(MatrixesLU& LU, vector<float>& b)
	{

		vector<float>* x = new vector<float>;
		vector<float>& xRef = *x;
		xRef = { 0,0,0,0 };

		vector<float> y = {0,0,0,0};
		float temp = 0;

#pragma region showingResults

		cout << endl;
		cout << endl;
		cout << "Macierz U:" << endl;

		for (unsigned int i = 0; i < LU.U.size(); i++)
		{
			for (unsigned int j = 0; j < LU.U[0].size(); j++)
				cout << setw(10) << LU.U[i][j] << "  ";
			cout << endl;
		}

		cout << endl;
		cout << endl;
		cout << "Macierz L:" << endl;

		for (unsigned int i = 0; i < LU.L.size(); i++)
		{
			for (unsigned int j = 0; j < LU.L[0].size(); j++)
				cout << setw(10) << LU.L[i][j] << "  ";
			cout << endl;
		}

		cout << endl;
		cout << endl;
#pragma endregion

		//solve Y
		y[0] = b[0];

		for (unsigned int i = 1; i < b.size(); i++)
		{
			temp = 0.0f;
			for (unsigned int j = 0; j < i; j++)
			{
				temp += LU.L[i][j] * y[j];
			}

			y[i] = b[i] - temp;
		}

		//solve X

		for (int i = b.size() - 1; i >= 0; i--)
		{
			temp = 0.0f;
			for (unsigned int j = i+1; j < b.size(); j++)
			{
				temp += LU.U[i][j] * xRef[j];
			}

			xRef[i] = (y[i] - temp) / LU.U[i][i];
		}


		return xRef;
	}

	void test()
	{

		vector<vector<float>> A =
		{
			{1.0f, -20.0f, 30.0f, -4.0f},
			{2.0f, -40.0f, -6.0f, 50.0f },
			{9.0f, -180.0f, 11.0f, -12.0f },
			{-16.0f, 15.0f, -140.0f, 13.0f}
		};

		vector<float> b = { 35.0f , 104.0f , -366.0f , -354.0f };

		MatrixesLU LU = getLU(A);
		vector<float> UB = getB(b, LU);

		vector<float> x = { 0,0,0,0 };


		x = solveX(LU,UB);

		cout << "X is equal to:\n";
		for (unsigned int i = 0; i < x.size(); i++)
		{
			cout<< setw(5) << x[i] << endl;
		}
		
	}

}