#include "Common.h"
#include <iostream>
#include <cstdlib>
#include <random>
#include <cmath>
#include <vector>
#include <stack>
#include <unordered_set>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>

using namespace std;
//int p[maxx][maxx];
const int k = 10;
const int N = 1304;
int edgenum = 0;
const double p = 0.30;

int main()
{
	int** adj = new int*[N];
	int** cij = new int*[N];

	for (int i = 0; i < N; ++i)
	{
		adj[i] = new int[N];
		cij[i] = new int[N];
	}
		

	boost::mt19937 gen;
	boost::random::uniform_real_distribution<> dist(0, 1);
	boost::random::uniform_int_distribution<int> edgedist(30, 50);

	for (int i = 0; i < N; ++i)
		for (int j = i + 1; j < N; ++j)
			if (dist(gen) < p)
			{
				cij[i][j] = edgedist(gen);
				adj[i][j] = 1;
				++edgenum;
			}

	//write-to-a file
	ofstream edgetofile("er-" + to_string(N) +  ".txt");

	edgetofile << N << "\t" << edgenum << "\t" << k << "\n";

	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j)
		{
			if (adj[i][j] == 1)
			{
				edgetofile << i << "\t" << j << "\t" << cij[i][j] << "\n";
			}
				
		}
		
	return 0;
}