//#include <bits/stdc++.h> 
using namespace std;
#include <iostream>
#include <windows.h>
#include <cstdlib>
#include <fstream>
#include <filesystem>
#include <random>
#include <unordered_set>
#include <boost/random/mersenne_twister.hpp>
#include <boost/graph/random.hpp>

int main()
{
	string pathloop = "C:\\Users\\murat\\Desktop\\sensornetwork\\networkoptim\\input\\otherdatalist.txt";
	string line;
	int seed = 1500;

	ifstream fileloop(pathloop);
	if (fileloop.is_open())
	{
		while (getline(fileloop, line))
		{
			int N;
			string filename;
			istringstream Streamloop(line);
			Streamloop >> N >> filename;

			boost::mt19937 gen(seed);
			boost::random::uniform_int_distribution<int> fdist(300, 400);

			ofstream ftofile("C:\\Users\\murat\\Desktop\\sensornetwork\\networkoptim\\input\\sensitivity\\" + filename + "\\f.txt");

			int* f = new int[N];

			for (int i = 0; i < N; ++i)
			{
				f[i] = fdist(gen);
				ftofile << f[i] << "\n";
			}

			delete f;
			++seed;
		}
	}
		
	return 0;
}