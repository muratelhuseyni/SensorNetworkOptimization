#include "Common.h"
#include <iostream>
#include <vector>

using namespace std;

// Function to perform DFS traversal
void dfs(int node, vector<vector<int>>& graph, vector<bool>& visited) {
	visited[node] = true;
	for (int adjacent : graph[node]) {
		if (!visited[adjacent]) {
			dfs(adjacent, graph, visited);
		}
	}
}

// Function to check if the graph is connected
bool isConnectedGraph(vector<vector<int>>& graph, int numNodes) 
{
	vector<bool> visited(numNodes, false);

	// Perform DFS traversal starting from node 0
	dfs(0, graph, visited);

	// Check if all nodes were visited
	for (bool visit : visited) 
	{
		if (!visit) 
			return false;  // If any node was not visited, graph is not connected
	}

	return true;  // If all nodes were visited, graph is connected
}

int main() 
{
	/////7
	int N;

	//string conlist = "datalist";
	string pathloop = "C:\\Users\\murat\\Desktop\\sensornetwork\\networkoptim\\datalist.txt";

	string line;
	string upperfolder = "albert-barabasi";
	ifstream fileloop(pathloop);
	ofstream resultfile("conresults_" + upperfolder + "_.txt");

	if (fileloop.is_open())
	{
		//loopun içini düzenle
		while (getline(fileloop, line))
		{
			//vector<vector<int>> graph;
			istringstream Streamloop(line); //ID pj rj dj Ymax
			string filename;
			Streamloop >> filename >> N;
			vector<vector<int>> graph(N);
			//graph = zraph;

			string path = "C:\\Users\\murat\\Desktop\\sensornetwork\\networkoptim\\input\\" +
				upperfolder + "\\" + filename + "\\" + filename + ".txt";
			
			ifstream file(path);
			if (file.is_open())
			{				

				while (getline(file, line))
				{
					istringstream Stream(line); //ID pj rj dj Ymax
					int u, v, val;
					Stream >> u >> v >> val;
					//graph[u-1].push_back(v);
					//graph[v-1].push_back(u);  // For undirected graph, add both edges
					graph[u].push_back(v);
					graph[v].push_back(u);  // For undirected graph, add both edges

				}

				if (isConnectedGraph(graph, N))
					resultfile << filename << "connected" << endl;
				else
					resultfile << filename << "disconnected" << endl;

				file.close();
			}

			//graph.empty();
		}
		
		fileloop.close();
	}
	
	return 0;
}
