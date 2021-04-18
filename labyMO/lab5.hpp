#include <iostream>
#include <iomanip>
#include <vector>

#pragma once

namespace lab5
{
	using namespace std;

	struct MatrixesLU
	{
		vector<vector<float>> L;
		vector<vector<float>> U;
	};



	MatrixesLU& getLU(vector<vector<float>> A);
	void test();

}