#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>

#pragma once

namespace lab5
{
	using namespace std;

	struct MatrixesLU
	{
		vector<vector<float>> L;
		vector<vector<float>> U;
		vector<int> changes;
	};



	MatrixesLU& getLU(vector<vector<float>> A);
	void test();

}