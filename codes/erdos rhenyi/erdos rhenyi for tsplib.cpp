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
int k, N, edgenum;
const double p = 0.30;

struct Vertex
{
	int I;
	int J;

	Vertex(int Iin, int Jin)
	{
		I = Iin;
		J = Jin;
	}
};

int main()
{
	//tsp225
	//string path = "C:\\Users\\hp\\Desktop\\son tubitak baþvurusu\\slovenya kod\\cplex tez\\tsp libraries\\d657.tsp\\d657.txt"; //upper triangular

	//1.8.23-for loop
	string pathloop = "C:\\Users\\hp\\Desktop\\son tubitak baþvurusu\\slovenya kod\\cplex tez\\input\\tsplist.txt";
	string line;

	ifstream fileloop(pathloop);
	if (fileloop.is_open())
	{
		while (getline(fileloop, line))
		{
			string filename;
			istringstream Streamloop(line);
			Streamloop >> filename;
			string path = "C:\\Users\\hp\\Desktop\\son tubitak baþvurusu\\slovenya kod\\cplex tez\\tsp libraries\\" + filename +  ".tsp\\" + filename + ".txt"; //upper triangular

			ifstream file(path);
			int edgenum = 0;
			double** cij;
			int** adj;
			string line;
			vector<Vertex> nodes;

			if (file.is_open())
			{
				int count = 0;
				int i = 0;
				int j = 0;

				while (getline(file, line))
				{
					istringstream Stream(line); //ID pj rj dj Ymax
					if (count == 0)
					{
						Stream >> N;

						cij = new double*[N];

						for (int i = 0; i < N; i++)
							cij[i] = new double[N];

						adj = new int*[N];

						for (int i = 0; i < N; i++)
							adj[i] = new int[N];

					}

					else
					{
						//tsp225
						//euclid distance					
						int i = -1;
						int j = -1;

						Stream >> i >> j;
						nodes.emplace_back(i, j);
						//tsp225
					}

					++count;
				}

				file.close();
			}

			edgenum = 0;
			boost::mt19937 gen;
			boost::random::uniform_real_distribution<> dist(0, 1);

			for (int i = 0; i < N; ++i)
				for (int j = i + 1; j < N; ++j)
					if (dist(gen) < p)
					{
						//// Calculation for distance matrix. Gets the euclidean distance
						cij[i][j] = sqrt(pow(nodes[i].I - nodes[j].I, 2) + pow(nodes[i].J - nodes[j].J, 2));
						++edgenum;
						adj[i][j] = 1;
					}

			//write-to-a file
			ofstream edgetofile(filename + ".txt");

			edgetofile << N << "\t" << edgenum << "\t" << k << "\n";

			for (int i = 0; i < N; ++i)
				for (int j = 0; j < N; ++j)
				{
					if (adj[i][j] == 1)
						edgetofile << i << "\t" << j << "\t" << cij[i][j] << "\n";
				}

		}
	}


		
	return 0;
}