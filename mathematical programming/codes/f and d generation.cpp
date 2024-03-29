#include "Common.h"
#include <iostream>
#include <random>
#include <unordered_set>
#include <boost/random/mersenne_twister.hpp>
#include <boost/graph/random.hpp>

int main()
{
	const int N = 1432;
	boost::mt19937 gen(3);
	boost::random::uniform_int_distribution<int> fdist(50000, 60000);
	boost::random::uniform_int_distribution<int> Ddist(100, 200);

	//write-to-a file
	ofstream ftofile("f.txt");
	ofstream dtofile("d.txt");

	int* f = new int[N];
	int* d = new int[N];

	for (int i = 0; i < N; ++i)
	{
		d[i] = Ddist(gen);
		dtofile << d[i] << "\n";
		f[i] = fdist(gen);
		ftofile << f[i] << "\n";
	}
	
	return 0;
}