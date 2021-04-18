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
		
		float ratio = 0;

		for (unsigned int i = 0; i < A[0].size(); i++)
		{
			if(fabs(result->U[i][i]) < 0.0001f)
				result->U[i].swap(result->U[getIndexRowOfMaxValueInCollumn(result->U, i)]);

			for (unsigned int y = 1+i; y < A.size(); y++)
			{

				cout << result->U[i][i];
				ratio = result->U[y][i] / result->U[i][i];
				

				result->L[y][i] = ratio;

				for (unsigned int x = i; x < A[0].size(); x++)
				{

					result->U[y][x] -= result->U[i][x] * ratio;

				}
				cout << endl;
			}
			
		}


		
		for (unsigned int i = 0; i < A.size(); i++)
		{
			for (unsigned int j = 0; j < A[0].size(); j++)
				cout << result->U[i][j] << "  ";
			cout << endl;
		}

		cout << endl; 
		cout << endl;

		for (unsigned int i = 0; i < A.size(); i++)
		{
			for (unsigned int j = 0; j < A[0].size(); j++)
				cout << result->L[i][j] << "  ";
			cout << endl;
		}
		
		cout << endl;
		cout << endl;

		vector<vector<float>> multiTest = {
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
				cout << multiTest[i][j] << "  ";
			cout << endl;
		}


		return *result;
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

		int rowChanges[4] = { 1,2,3,4 };

		MatrixesLU LU = getLU(A);

	}

}