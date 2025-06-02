#include "Common.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <stack>
#include <chrono>
#include <unordered_set>
#include<array> 
#include <set>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>

using namespace boost;

typedef adjacency_list<vecS, vecS, directedS,
	no_property, property<edge_weight_t, int >> GraphDijkstra;

typedef graph_traits<GraphDijkstra>::vertex_descriptor Vertex;

typedef std::pair<int, int> Edge;

struct DijkstraResults 
{
	std::vector<Vertex> path;
	int weight;
};

struct Path
{
	int i;
	int j;
	int val;
	int compi;
	int compj;

	Path(int iIn, int jIn, int valIn, int compiIn, int compjIn) //constructor
	{
		i = iIn;
		j = jIn;
		val = valIn;
		compi = compiIn;
		compj = compjIn;
	}
};

struct RSPath
{
	int i;
	int j;
	vector<int> ps;
	vector<int> Cps;

	RSPath(int iIn, int jIn, vector<int> psIn, vector<int> CpsIn) //constructor
	{
		i = iIn;
		j = jIn;
		ps = psIn;
		Cps = CpsIn;
	}
};

vector<int> FindCps(vector<int> R, vector<int> PS, vector<vector<int>>& S)
{
	vector<int> SR;
	vector<int> SPS;

	for (int i = 0; i < R.size(); ++i)
	{
		if (i == 0)
			SR = S[R[i]];
		else
		{
			int facil = R[i];
			for (int client : S[facil])
			{
				std::vector < int > ::iterator itr;
				itr = find(SR.begin(), SR.end(), client);
				if (itr != SR.end()) //found
					continue;
				SR.push_back(client);
			}
		}
	}

	for (int i = 0; i < PS.size(); ++i)
	{
		if (i == 0)
			SPS = S[PS[i]];
		else
		{
			int facil = PS[i];
			for (int client : S[facil])
			{
				std::vector < int > ::iterator itr;
				itr = find(SPS.begin(), SPS.end(), client);
				if (itr != SPS.end()) //found
					continue;
				SPS.push_back(client);
			}
		}
	}

	vector<int> CPS;
	//S_PS - S_R
	for (int clientps : SPS)
	{
		std::vector < int > ::iterator itr;
		itr = find(SR.begin(), SR.end(), clientps);
		if (itr == SR.end()) //notfound
			CPS.push_back(clientps);
	}

	return CPS;
}

vector<int> FindPs(DijkstraResults&route)
{
	vector<int> PS;
	for (int i = 1; i < route.path.size(); ++i)
		PS.push_back(route.path[i]);

	return PS;
}


void Repair(vector<pair<int, int>>& edges, int N, vector<int>& vertices, int** cij, int p)
{
	while (vertices.size() < p)
	{
		int min = INT_MAX;
		int insel, outsel;

		for (int i = 0; i < N; ++i)
		{
			if (find(vertices.begin(), vertices.end(), i)
				!= vertices.end())
				continue;
			//then i is adjacent-outside
			for (int inside : vertices)
			{
				if (cij[inside][i] < min)
				{
					min = cij[inside][i];
					insel = inside;
					outsel = i;
				}
			}
		}

		vertices.push_back(outsel);
		edges.push_back({ insel, outsel });
	}

}

bool FindIntersection(vector<int>& arr1, vector<int>& arr2)
{
	bool intersect = false;
	for (int i : arr1)
		for (int j : arr2)
			if (i == j)
			{
				intersect = true;
				return intersect;
			}


	return intersect;
}

vector<vector<int>> MinSpanTree(int** graph, int numNodes)
{
	// Create arrays to keep track of visited nodes and minimum edge weights
	bool* visited = new bool[numNodes];
	int* minWeights = new int[numNodes];

	// Initialize the visited array and minimum weights
	for (int i = 0; i < numNodes; i++) {
		visited[i] = false;
		minWeights[i] = INT_MAX;
	}

	// Start with the first node
	minWeights[0] = 0;

	// Variables to keep track of the minimum spanning tree
	int* parent = new int[numNodes];
	parent[0] = -1;

	// Find the minimum spanning tree
	for (int i = 0; i < numNodes - 1; i++) {
		int minNode = -1;
		for (int j = 0; j < numNodes; j++) {
			if (!visited[j] && (minNode == -1 || minWeights[j] < minWeights[minNode])) {
				minNode = j;
			}
		}

		visited[minNode] = true;

		for (int j = 0; j < numNodes; j++) {
			if (graph[minNode][j] && !visited[j] && graph[minNode][j] < minWeights[j]) {
				parent[j] = minNode;
				minWeights[j] = graph[minNode][j];
			}
		}
	}

	vector<vector<int>> Edges;
	Edges.resize(numNodes);

	// Print the minimum spanning tree
	//std::cout << "Edges in the Minimum Spanning Tree:" << std::endl;
	for (int i = 1; i < numNodes; i++)
	{
		if (parent[i]<0 || parent[i]>numNodes)
		{
			Edges.clear();
			Edges.resize(0);
			return Edges;
		}

		//std::cout << "Edge " << parent[i] << " - " << i << " with weight " << graph[i][parent[i]] << std::endl;
		//Edges[i].emplace_back(parent[i], i, graph[i][parent[i]]);
		Edges[i].push_back(parent[i]);
		Edges[i].push_back(i);
		Edges[i].push_back(graph[i][parent[i]]);
	}

	return Edges;
}

DijkstraResults dijkstra_shortest_path(const GraphDijkstra& g, Vertex start, Vertex goal)
{
	// Set up Dijkstra's algorithm
	std::vector<Vertex> parents(num_vertices(g), graph_traits<GraphDijkstra>::null_vertex());
	std::vector<int> distances(num_vertices(g), std::numeric_limits<int>::max());

	// Apply Dijkstra's algorithm
	dijkstra_shortest_paths(
		g, start,
		predecessor_map(make_iterator_property_map(parents.begin(), get(vertex_index, g))).
		distance_map(make_iterator_property_map(distances.begin(), get(vertex_index, g)))
		);

	// Check if the goal is reachable
	if (distances[goal] == std::numeric_limits<int>::max())
	{
		// Goal is unreachable
		return{ {}, -1 }; // Return empty path and invalid distance
	}

	// Reconstruct the path
	std::vector<Vertex> path;
	for (Vertex v = goal; v != graph_traits<GraphDijkstra>::null_vertex(); v = parents[v])
	{
		path.push_back(v);
		if (v == start)
			break;
	}

	std::reverse(path.begin(), path.end());
	return{ path, distances[goal] };
}


using namespace std;
//int p[maxx][maxx];
int k, N, edgenum;

//k_5, N=250, after .. sec, addition of D makes it harder to solve - 4% after 10 mins

// Function to find all subsets of given set. Any repeated
// subset is considered only once in the output
vector<vector<int> > findPowerSet(vector<int>& nums)
{
	// Size of array to set bit
	int bits = nums.size();

	// Total number of subsets = pow(2,
	// sizeof(arr))
	unsigned int pow_set_size = pow(2, bits);

	// Sort to avoid adding permutation of
	// the substring
	sort(nums.begin(), nums.end());
	vector<vector<int> > ans;

	// To store subset as a list to
	// avoid adding exact duplicates
	vector<string> list;

	// Counter 000..0 to 111..1
	for (int counter = 0; counter < pow_set_size;
		counter++) {
		vector<int> subset;
		string temp = "";

		// Check for the current bit in the counter
		for (int j = 0; j < bits; j++) {
			if (counter & (1 << j)) {
				subset.push_back(nums[j]);

				// Add special character to separate
				// integers
				temp += to_string(nums[j]) + '$';
			}
		}

		if (find(list.begin(), list.end(), temp)
			== list.end()) {
			ans.push_back(subset);
			list.push_back(temp);
		}
	}

	return ans;
}

void find_k_element_subsets(vector<int>& set, int k, int start, vector<int>& subset, vector<vector<int>>& result) {
	if (subset.size() == k) {
		result.push_back(subset);
		return;
	}

	for (int i = start; i < set.size(); i++) {
		subset.push_back(set[i]);
		find_k_element_subsets(set, k, i + 1, subset, result);
		subset.pop_back();
	}
}

vector<vector<int>> get_k_element_subsets(vector<int>& set, int k) {
	vector<vector<int>> result;
	vector<int> subset;
	find_k_element_subsets(set, k, 0, subset, result);
	return result;
}

//floydmarshal
void allpairshort(int** a, int n)
{
	/*if (dist[i][j] > (dist[i][k] + dist[k][j])
		&& (dist[k][j] != INF
		&& dist[i][k] != INF))*/

	int k, i, j;
	for (k = 0; k < n; k++)
	{
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				if (a[i][k] + a[k][j] < a[i][j])
				{
					a[i][j] = a[i][k] + a[k][j];
				}
			}
		}
	}

	/*for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			cout << a[i][j] << endl;*/

}

// Define a custom comparison function to sort pairs by value in descending order
bool comparePairs(const pair<int, int>& a, const pair<int, int>& b) {
	return a.second > b.second;
}

bool compareBySize(const std::vector<int>& a, const std::vector<int>& b)
{
	return a.size() > b.size();
}

vector<int> FindCover(vector<vector<int>>& S, vector<pair<int, double>> sortedchead)
{
	vector<int> U;
	vector<int> V;
	int N = S.size();

	for (int i = 0; i < N; ++i)
		V.push_back(i);
	
	for (const auto& pair : sortedchead)
	{
		int member = pair.first;
		int count = 0;
		for (int i : S[member])
			if (binary_search(V.begin(), V.end(), i))
			{
				if (count == 0)
					U.push_back(member);
				++count;
				V.erase(std::remove(V.begin(), V.end(), i), V.end());
				if (V.empty())
					return U;
			}

	}

	return U;
}

vector<int> GetCandidateFacil(vector<vector<int>>& S, vector<int> V)
{
	vector<int> U;
	int N = S.size();

	//sort(V.begin(), V.end());

	for (int j = 0; j < N; ++j)
	{
		if (S[j].size() == 0)
			continue;

		int count = 0;

		for (int i : S[j])
			if (binary_search(V.begin(), V.end(), i))
			{
				++count;
				if (count == 1)
					U.push_back(j);
				V.erase(std::remove(V.begin(), V.end(), i), V.end());
				if (V.empty())
					return U;
			}
	}

	return U;
}

void FillNoncoveredclients(vector<vector<int>>& S, vector<int>& vertices, int** cij)
{
	vector<int> U;
	int N = S.size();

	vector<int> V;
	for (int i = 0; i < N; ++i)
		V.push_back(i);

	//sort(V.begin(), V.end());

	for (int j : vertices)
	{
		if (S[j].size() == 0)
			continue;

		for (int i : S[j])
			if (binary_search(V.begin(), V.end(), i))
			{
				V.erase(std::remove(V.begin(), V.end(), i), V.end());
				if (V.empty())
					return;
			}
	}

	if (!V.empty())
	{
		//find minimum cij from noncovered client to candidate facility
		for (int i : V)
		{
			int min = INT_MAX;
			int node = -1;
			for (int j : vertices)
			{
				if (cij[i][j] < min)
				{
					min = cij[i][j];
					node = j;
				}
			}
			S[node].push_back(i);
		}
	}
}

void ConnectedSetCover(vector<pair<int, int>>& edges, vector<int>& vertices, vector<vector<int>>& S)
{
	int N = S.size();
	vector<int> V;
	for (int i = 0; i < N; ++i)
		V.push_back(i);

	vector<int> R;
	vector<int> U;

	vector<pair<int, int>> Redges; //edges.push_back({ u, v });

	//find adjacents: cover adjacent union graph-adjacent
	int** adj;
	adj = new int*[N];

	for (int i = 0; i < N; i++)
		adj[i] = new int[N];

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			if (i == j)
				adj[i][j] = 0;
			else
				adj[i][j] = INT_MAX;

	for (const auto& pair : edges)
	{
		adj[pair.first][pair.second] = 1;
		adj[pair.second][pair.first] = 1;
	}

	//Initialitzation

	//find max s_j
	int max = INT_MIN;
	int node = -1;

	for (int j : vertices)
	{
		int size = S[j].size();
		if (size > max)
		{
			max = size;
			node = j;
		}
	}

	R.push_back(node);
	U = S[node];
	//deduct S_node from V
	for (int i : U)
		if (binary_search(V.begin(), V.end(), i))
			V.erase(std::remove(V.begin(), V.end(), i), V.end());

	//step 2.1
	vertices.erase(std::remove(vertices.begin(), vertices.end(), node), vertices.end()); //vertices=S-R	

	while (V.size() > 0)
	{
		//Adjacents
		//get dist 1 ones
		vector<int> adjacent;

		for (int inside : R)
			for (int outside : vertices)
			{
				if (adj[inside][outside] == 1)
					adjacent.push_back(outside);
			}

		//foreach S_j and S_i set

		for (int inside : R)
			for (int outside : vertices)
				if (adj[inside][outside] > 1)
				{
					bool find = FindIntersection(S[inside], S[outside]);
					if (find)
						adjacent.push_back(outside);
				}

		//find shortest R-S path - use dijkstra shortest path
		GraphDijkstra gd(N);

		for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j)
			if (adj[i][j] == 1)
				add_edge(i, j, 1, gd); //1 for unit graph, O.W use c_ij

		vector<RSPath> RSpaths;
		for (int outside : adjacent)
		{
			vector<RSPath> endpoints;
			for (int inside : R)
			{
				DijkstraResults route = dijkstra_shortest_path(gd, inside, outside);
				vector<int> PS = FindPs(route);
				vector<int> CPS = FindCps(R, PS, S);
				endpoints.emplace_back(inside, outside, PS, CPS);
			}

			int min = INT_MAX;
			int index = -1;
			for (int i = 0; i < endpoints.size(); ++i)
			{
				RSPath rspath = endpoints[i];
				if (rspath.ps.size() < min)
				{
					min = rspath.ps.size();
					index = i;
				}

			}

			RSpaths.push_back(endpoints[index]);
		}

		//step 2.2
		//find min weight e(PS) ratio
		double min = INT_MAX;
		int index = -1;
		for (int i = 0; i < RSpaths.size(); ++i)
		{
			RSPath rspath = RSpaths[i];
			if (rspath.Cps.size() == 0)
				continue;
			double eps = rspath.ps.size() / ((double)rspath.Cps.size());
			if (eps < min)
			{
				min = eps;
				index = i;
			}
		}

		//add min eps as PS and add it to R and U
		if (index <0)
			break;
		RSPath Pspath = RSpaths[index]; //Redges;
		vector<int> PS = Pspath.ps;
		int source = Pspath.i;

		Redges.push_back({ source, PS[0] });


		if (PS.size() == 1)
		{
			R.push_back(PS[0]);
			vertices.erase(std::remove(vertices.begin(), vertices.end(), PS[0]), vertices.end());
		}

		//add this path to Redges
		for (int i = 0; i < PS.size() - 1; ++i)
		{
			int u = PS[i];
			int v = PS[i + 1];
			Redges.push_back({ u, v });
			R.push_back(u);
			vertices.erase(std::remove(vertices.begin(), vertices.end(), u), vertices.end());

			if (v == PS[PS.size() - 1])
			{
				R.push_back(v);
				vertices.erase(std::remove(vertices.begin(), vertices.end(), v), vertices.end());
			}
		}

		for (int i : Pspath.Cps)
			if (binary_search(V.begin(), V.end(), i))
			{
				V.erase(std::remove(V.begin(), V.end(), i), V.end());
				if (V.empty())
				{
					edges = Redges;
					vertices = R;
					return;
				}

			}

	}

	edges = Redges;
	vertices = R;
	return;
}

class Graph {
public:
	Graph(int vertices) : V(vertices), adj(vertices) {}

	void addEdge(int u, int v)
	{
		if (adj[u].find(v) == adj[u].end()) //checks whether u and v are adjacent
		{
			adj[u].insert(v);
			adj[v].insert(u);
		}
	}

	/*void addEdgeCharikar(int u, int v)
	{
	adj[u].insert(v);
	adj[v].insert(u);
	}*/

	set<int> findDominatingSet(const set<int>& requiredVertices) {
		set<int> dominatingSet = requiredVertices;  // Start with required vertices
		vector<bool> dominated(V, false);

		// Mark vertices dominated by required vertices
		for (int v : requiredVertices) {
			dominated[v] = true;
			for (int neighbor : adj[v]) {
				dominated[neighbor] = true;
			}
		}

		// Get list of undominated vertices
		vector<int> undominated;
		for (int i = 0; i < V; ++i) {
			if (!dominated[i] && !adj[i].empty()) {  // Only consider vertices with edges
				undominated.push_back(i);
			}
		}

		// While there are undominated vertices
		while (!undominated.empty()) {
			// Find vertex that dominates maximum number of undominated vertices
			int bestVertex = -1;
			int maxNewDominated = -1;

			for (int i = 0; i < V; ++i) {
				if (dominatingSet.count(i) > 0) continue;  // Skip if already in dominating set

				// Count how many undominated vertices this vertex would dominate
				int newDominated = 0;
				for (int neighbor : adj[i]) {
					if (!dominated[neighbor]) {
						newDominated++;
					}
				}
				if (!dominated[i]) newDominated++;  // Count self if undominated

				if (newDominated > maxNewDominated) {
					maxNewDominated = newDominated;
					bestVertex = i;
				}
			}

			if (bestVertex == -1) break;  // No improvement possible

			// Add best vertex to dominating set
			dominatingSet.insert(bestVertex);

			// Mark newly dominated vertices
			dominated[bestVertex] = true;
			for (int neighbor : adj[bestVertex]) {
				dominated[neighbor] = true;
			}

			// Update undominated list
			undominated.clear();
			for (int i = 0; i < V; ++i) {
				if (!dominated[i] && !adj[i].empty()) {
					undominated.push_back(i);
				}
			}
		}

		return dominatingSet;
	}


	bool DoesExist(int u, int v)
	{
		return adj[u].find(v) != adj[u].end();
	}

	vector<set<int>> findConnectedComponents() {
		vector<bool> visited(V, false);
		vector<set<int>> components;

		for (int v = 0; v < V; ++v) {
			if (!visited[v]) {
				set<int> component;
				dfs(v, visited, component);
				//new
				if (component.size() == 1)
					continue;
				components.push_back(component);
			}
		}

		return components;
	}

	vector<pair<int, int>> getAllEdges()
	{
		vector<pair<int, int>> edges;
		for (int u = 0; u < V; ++u)
		{
			for (int v : adj[u])
			{
				if (u < v)
					edges.push_back({ u, v });
			}
		}
		return edges;
	}

	vector<int> getAllVertices()
	{
		vector<int> vertices;
		for (int v = 0; v < V; ++v)
			if (!adj[v].empty() && find(vertices.begin(), vertices.end(), v) == vertices.end())
				vertices.push_back(v);

		return vertices;
	}

	vector<pair<int, int>> getEdgeList(const set<int>& vertices)
	{
		vector<pair<int, int>> edges;
		for (int v : vertices)
		{
			for (int neighbor : adj[v])
			{
				// Ensure that the edge is not duplicated and both vertices are in the given set
				if (v < neighbor && vertices.count(neighbor) > 0) {
					edges.emplace_back(v, neighbor);
				}
			}
		}
		return edges;
	}

private:
	int V;  // Number of vertices
	vector<set<int>> adj;

	void dfs(int v, vector<bool>& visited, set<int>& component) {
		visited[v] = true;
		component.insert(v);

		for (int neighbor : adj[v]) {
			if (!visited[neighbor]) {
				dfs(neighbor, visited, component);
			}
		}
	}

	// Function to check if a given set is a dominating set
	bool isDominatingSet(const set<int>& subset) {
		vector<bool> dominated(V, false);  // Tracks which vertices are dominated

		// Mark all vertices in the subset and their neighbors as dominated
		for (int u : subset) {
			dominated[u] = true;  // The vertex itself is dominated
			for (int v : adj[u]) {  // Mark all neighbors of u as dominated
				dominated[v] = true;
			}
		}

		// Check if all vertices in the graph are dominated
		for (int v = 0; v < V; ++v) {
			if (!dominated[v]) {  // If any vertex is not dominated, return false
				return false;
			}
		}

		return true;  // All vertices are dominated
	}

};

void GetData(ifstream& file, ifstream& file2, ifstream& file3, string line,  int* d, int* f, int** cij, int** adj)
{
	//for orlib
	if (file.is_open())
	{
		while (getline(file, line))
		{
			istringstream Stream(line); //ID pj rj dj Ymax
			int i, j, val;
			Stream >> i >> j >> val;
			//4.4.24-control!
			if (i == j)
				continue;
			cij[i][j] = val;
			cij[j][i] = cij[i][j];
			adj[i][j] = 1;
			adj[j][i] = 1;
			//for dijkstra
			//add_edge(i, j, val, gd); //for weighted graph
//			add_edge(i, j, 1, gd); //for unit graph
		}

		file.close();
	}



	//for orlib
	if (file2.is_open())
	{
		int i = 0;
		while (getline(file2, line))
		{
			int val;
			istringstream Stream(line); //ID pj rj dj Ymax
			Stream >> val;
			d[i] = val;
			++i;
		}

		file2.close();
	}




	//for orlib
	if (file3.is_open())
	{
		int i = 0;
		while (getline(file3, line))
		{
			int val;
			istringstream Stream(line); //ID pj rj dj Ymax
			Stream >> val;
			f[i] = val;
			++i;
		}

		file3.close();
	}
}

void SolveMTZ(int** cij, int** conij, int** adj, int* d, int* f, int N, int M, int k, const vector<pair<int, int>>& edges, const vector<int>& vertices, vector<double>& output)
{
	IloEnv env;
	IloModel model(env);
	BoolVarArray2 x = CreateBoolVarArray2(env, N, N, "x");
	BoolVarArray2 Z = CreateBoolVarArray2(env, N, N, "Z");
	IloBoolVarArray y = CreateBoolVarArray(env, N, "y");
	IloNumVarArray u = CreateNumVarArray(env, N, "u", 0, IloInfinity);

	IloExpr objExp(env);

	int** adj2;
	adj2 = new int*[N];

	for (int i = 0; i < N; i++)
		adj2[i] = new int[N];

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
		{
			if (i == j)
				adj2[i][j] = 0;
			else
				adj2[i][j] = 100000;
		}

	for (auto edge : edges)
	{
		int i = edge.first;
		int j = edge.second;
		adj2[i][j] = 1;
		adj2[j][i] = 1;
	}

	for (int i = 0; i < N; ++i)
		for (int j : vertices)
			objExp += d[i] * cij[i][j] * x[i][j];

	//fixed opening cost
	for (int i : vertices)
		objExp += f[i] * y[i];

	for (int i : vertices)
		for (int j : vertices)
			if (adj2[i][j] == 1)
				objExp += M *conij[i][j] * Z[i][j];

	IloObjective obj = IloMinimize(env, objExp, "obj");
	model.add(obj);
	objExp.end();

	//Baseline
	// Coverage constraint		
	for (int i = 0; i < N; ++i)
	{
		IloExpr costcover(env);
		for (int j : vertices)
			costcover += x[i][j];
		//30.1.23
		model.add(costcover == 1);
		costcover.end();
	}

	for (int i = 0; i < N; ++i)
		for (int j : vertices)
			model.add(x[i][j] <= y[j]);

	//select k facility
	IloExpr facilsum(env);
	for (int i : vertices) //foreach
		facilsum += y[i];
	model.add(facilsum == k);
	facilsum.end();

	//new-2-min span tree size
	//2-37
	IloExpr z_ijextendedtree(env);
	for (int i : vertices)
		for (int j : vertices)
			if (adj2[i][j] == 1)
				z_ijextendedtree += Z[i][j];

	model.add(z_ijextendedtree == k - 1);
	z_ijextendedtree.end();

	//at most one incoming
	for (int i : vertices)
	{
		IloExpr z_ijincoming(env);
		for (int j : vertices)
			if (adj2[j][i] == 1)
				z_ijincoming += Z[j][i];

		model.add(z_ijincoming <= y[i]);
		z_ijincoming.end();
	}

	//validineq1
	for (int i = 0; i < N; ++i)
		for (int j = i + 1; j < N; ++j)
			if (adj2[i][j] == 1)
				model.add(Z[i][j] + Z[j][i] <= y[i]);


	//validineq2
	for (int i = 0; i < N; ++i)
		for (int j = i + 1; j < N; ++j)
			if (adj2[i][j] == 1)
				model.add(Z[i][j] + Z[j][i] <= y[j]);


	//adjusted version wrt bileteral flow
	for (int i : vertices)
	{
		IloExpr Zijneighborsum(env);
		for (int j : vertices)
			if (adj2[i][j] == 1)
				Zijneighborsum += Z[i][j] + Z[j][i];

		model.add(y[i] <= Zijneighborsum);
		Zijneighborsum.end();
	}

	//valid inequality
	for (int i : vertices)
	{
		IloExpr yneighborsum(env);
		for (int j : vertices)
			if (adj2[i][j] == 1)
				yneighborsum += y[j];

		model.add(y[i] <= yneighborsum);
		yneighborsum.end();
	}

	for (int i : vertices)
		for (int j : vertices)
			if (adj2[i][j] == 1)
				model.add(u[i] - u[j] + k*Z[i][j] <= k - 1);

	IloCplex cplex(model);
	double cplextotal = 60 * 60 * 2;
	cplex.setParam(IloCplex::TiLim, cplextotal);

	IloBool success = cplex.solve();
	double lowerbound = pow(10, -5);

	if (success && cplex.isPrimalFeasible())
	{
		//transportation cost
		double trans = 0;
		double con = 0;
		double op = 0;

		//transport cost
		for (int i = 0; i < N; ++i)
			for (int j : vertices)
				trans += d[i] * cij[i][j] * cplex.getValue(x[i][j]);

		//opening cost
		for (int i : vertices)
			op += f[i] * cplex.getValue(y[i]);

		//connection cost
		for (int i = 0; i < N; ++i)
			for (int j = 0; j < N; ++j)
				if (adj2[i][j] == 1)
					con += M *conij[i][j] * cplex.getValue(Z[i][j]);

		output.push_back(trans);
		output.push_back(con);
		output.push_back(op);
		output.push_back(cplex.getObjValue());
	}

	else
	{
		output.push_back(-1.0);
		output.push_back(-1.0);
	}

	cplex.end();
	model.end();
	env.end();


}

void SolveFlow(int** cij, int** conij, int** adj, int* d, int* f, int N, int M, int k, const vector<pair<int, int>>& edges, const vector<int>& vertices, vector<double>& output)
{
	IloEnv env;
	IloModel model(env);
	NumVarArray2 a = CreateNumVarArray2(env, N, N, "a", 0, IloInfinity);
	IloBoolVarArray s = CreateBoolVarArray(env, N, "s"); //flow source
	BoolVarArray2 x = CreateBoolVarArray2(env, N, N, "x");
	BoolVarArray2 Z = CreateBoolVarArray2(env, N, N, "Z");
	IloBoolVarArray y = CreateBoolVarArray(env, N, "y");

	IloExpr objExp(env);

	int** adj2;
	adj2 = new int*[N];

	for (int i = 0; i < N; i++)
		adj2[i] = new int[N];

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
		{
			if (i == j)
				adj2[i][j] = 0;
			else
				adj2[i][j] = 100000;
		}

	for (auto edge : edges)
	{
		int i = edge.first;
		int j = edge.second;
		adj2[i][j] = 1;
		adj2[j][i] = 1;
	}

	for (int i = 0; i < N; ++i)
		for (int j : vertices)
			objExp += d[i] * cij[i][j] * x[i][j];

	//fixed opening cost
	for (int i : vertices)
		objExp += f[i] * y[i];

	for (int i : vertices)
		for (int j : vertices)
			if (adj2[i][j] == 1)
				objExp += M *conij[i][j] * Z[i][j];

	IloObjective obj = IloMinimize(env, objExp, "obj");
	model.add(obj);
	objExp.end();

	//Baseline
	// Coverage constraint		
	for (int i = 0; i < N; ++i)
	{
		IloExpr costcover(env);
		for (int j : vertices)
			costcover += x[i][j];
		//30.1.23
		model.add(costcover == 1);
		costcover.end();
	}

	for (int i = 0; i < N; ++i)
		for (int j : vertices)
			model.add(x[i][j] <= y[j]);

	//select k facility
	IloExpr facilsum(env);
	for (int i : vertices) //foreach
		facilsum += y[i];
	model.add(facilsum == k);
	facilsum.end();

	//new-2-min span tree size
	//2-37
	IloExpr z_ijextendedtree(env);
	for (int i : vertices)
		for (int j : vertices)
			if (adj2[i][j] == 1)
				z_ijextendedtree += Z[i][j];

	model.add(z_ijextendedtree == k - 1);
	z_ijextendedtree.end();

	//at most one incoming
	for (int i : vertices)
	{
		IloExpr z_ijincoming(env);
		for (int j : vertices)
			if (adj2[j][i] == 1)
				z_ijincoming += Z[j][i];

		model.add(z_ijincoming <= y[i]);
		z_ijincoming.end();
	}

	//validineq1
	for (int i = 0; i < N; ++i)
		for (int j = i + 1; j < N; ++j)
			if (adj2[i][j] == 1)
				model.add(Z[i][j] + Z[j][i] <= y[i]);


	//validineq2
	for (int i = 0; i < N; ++i)
		for (int j = i + 1; j < N; ++j)
			if (adj2[i][j] == 1)
				model.add(Z[i][j] + Z[j][i] <= y[j]);


	//adjusted version wrt bileteral flow
	for (int i : vertices)
	{
		IloExpr Zijneighborsum(env);
		for (int j : vertices)
			if (adj2[i][j] == 1)
				Zijneighborsum += Z[i][j] + Z[j][i];

		model.add(y[i] <= Zijneighborsum);
		Zijneighborsum.end();
	}

	//valid inequality
	for (int i : vertices)
	{
		IloExpr yneighborsum(env);
		for (int j : vertices)
			if (adj2[i][j] == 1)
				yneighborsum += y[j];

		model.add(y[i] <= yneighborsum);
		yneighborsum.end();
	}


	//Flow constraints
	//17a
	for (int j : vertices)
		model.add(s[j] <= y[j]);

	//17b
	IloExpr dummyexpr2(env);
	for (int j : vertices)
		dummyexpr2 += s[j];
	model.add(dummyexpr2 == 1);
	dummyexpr2.end();

	////11.10.23
	for (int i : vertices)
	{
		IloExpr inflowexpr(env);
		for (int j : vertices)
			if (adj2[j][i] == 1)
				inflowexpr += a[j][i];

		IloExpr outflowexpr(env);
		for (int j : vertices)
			if (adj2[i][j] == 1)
				outflowexpr += a[i][j];

		model.add(inflowexpr - outflowexpr == y[i] - k*s[i]);
		inflowexpr.end();
		outflowexpr.end();
	}

	for (int i : vertices)
		for (int j : vertices)
			if (adj2[i][j] == 1)
				model.add(a[i][j] <= (k - 1)*Z[i][j]);

	IloCplex cplex(model);
	double cplextotal = 60 * 60 * 2;
	cplex.setParam(IloCplex::TiLim, cplextotal);

	IloBool success = cplex.solve();

	if (success && cplex.isPrimalFeasible())
	{
		//transportation cost
		double trans = 0;
		double con = 0;
		double op = 0;

		//transport cost
		for (int i = 0; i < N; ++i)
			for (int j : vertices)
				trans += d[i] * cij[i][j] * cplex.getValue(x[i][j]);

		//opening cost
		for (int i : vertices)
			op += f[i] * cplex.getValue(y[i]);

		//connection cost
		for (int i = 0; i < N; ++i)
			for (int j = 0; j < N; ++j)
				if (adj2[i][j] == 1)
					con += M *conij[i][j] * cplex.getValue(Z[i][j]);

		output.push_back(trans);
		output.push_back(con);
		output.push_back(op);
		output.push_back(cplex.getObjValue());
	}

	else
	{
		output.push_back(-1.0);
		output.push_back(-1.0);
	}

	cplex.end();
	model.end();
	env.end();

}

double LPMTZSolve(int** cij, int** conij, int** adj, int* d, int* f, vector<std::tuple<int, int, double>>& Zhead, vector<std::tuple<int, int, double>>& Xhead, vector<pair<int, double>>& Yhead, int N, int M, int k)
{
	IloEnv env;
	IloModel model(env);
	BoolVarArray2 x = CreateBoolVarArray2(env, N, N, "x");
	BoolVarArray2 Z = CreateBoolVarArray2(env, N, N, "Z");
	IloBoolVarArray y = CreateBoolVarArray(env, N, "y");
	IloNumVarArray u = CreateNumVarArray(env, N, "u", 0, IloInfinity);

	IloConversion yconversion(env, y, ILOFLOAT);
	model.add(yconversion);

	for (int i = 0; i < N; ++i)
	{
		IloConversion Zconversion(env, Z[i], ILOFLOAT);
		model.add(Zconversion);
		IloConversion xconversion(env, x[i], ILOFLOAT);
		model.add(xconversion);
	}

	IloExpr objExp(env);

	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j)
			objExp += d[i] * cij[i][j] * x[i][j];

	//fixed opening cost
	for (int i = 0; i < N; ++i)
		objExp += f[i] * y[i];

	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j)
			if (adj[i][j] == 1)
				objExp += M *conij[i][j] * Z[i][j];

	IloObjective obj = IloMinimize(env, objExp, "obj");
	model.add(obj);
	objExp.end();

	//Baseline
	// Coverage constraint		
	for (int i = 0; i < N; ++i)
	{
		IloExpr costcover(env);
		for (int j = 0; j < N; ++j)
			costcover += x[i][j];
		//30.1.23
		model.add(costcover == 1);
		costcover.end();
	}

	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
			model.add(x[i][j] <= y[j]);
	}

	//select k facility
	IloExpr facilsum(env);
	for (int i = 0; i < N; ++i) //foreach
		facilsum += y[i];
	model.add(facilsum == k);
	facilsum.end();

	//new-2-min span tree size
	//2-37
	IloExpr z_ijextendedtree(env);
	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j)
			if (adj[i][j] == 1)
				z_ijextendedtree += Z[i][j];

	model.add(z_ijextendedtree == k - 1);
	z_ijextendedtree.end();

	//at most one incoming
	for (int i = 0; i < N; ++i)
	{
		IloExpr z_ijincoming(env);
		for (int j = 0; j < N; ++j)
			if (adj[j][i] == 1)
				z_ijincoming += Z[j][i];

		model.add(z_ijincoming <= y[i]);
		z_ijincoming.end();
	}

	//validineq1
	for (int i = 0; i < N; ++i)
		for (int j = i + 1; j < N; ++j)
			if (adj[i][j] == 1)
				model.add(Z[i][j] + Z[j][i] <= y[i]);


	//validineq2
	for (int i = 0; i < N; ++i)
		for (int j = i + 1; j < N; ++j)
			if (adj[i][j] == 1)
				model.add(Z[i][j] + Z[j][i] <= y[j]);


	//adjusted version wrt bileteral flow
	for (int i = 0; i < N; ++i)
	{
		IloExpr Zijneighborsum(env);
		for (int j = 0; j < N; ++j)
			if (adj[i][j] == 1)
				Zijneighborsum += Z[i][j] + Z[j][i];

		model.add(y[i] <= Zijneighborsum);
		Zijneighborsum.end();
	}

	//valid inequality
	for (int i = 0; i < N; ++i)
	{
		IloExpr yneighborsum(env);
		for (int j = 0; j < N; ++j)
			if (adj[i][j] == 1)
				yneighborsum += y[j];

		model.add(y[i] <= yneighborsum);
		yneighborsum.end();
	}

	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j)
			if (adj[i][j] == 1)
				model.add(u[i] - u[j] + k*Z[i][j] <= k - 1);

	IloCplex cplex(model);
	double cplextotal = 7813;
	cplex.setParam(IloCplex::TiLim, cplextotal);

	IloBool success = cplex.solve();
	double zlp;

	if (success && cplex.isPrimalFeasible())
	{
		for (int i = 0; i < N; ++i)
		  for (int j = 0; j < N; ++j)
			  Xhead.emplace_back(i, j, cplex.getValue(x[i][j]) );

		for (int i = 0; i < N; ++i)
			Yhead.push_back({ i, cplex.getValue(y[i]) });

		for (int i = 0; i < N; ++i)
			for (int j = 0; j < N; ++j)
				if (adj[i][j] == 1)
					Zhead.emplace_back(i, j, cplex.getValue(Z[i][j]));;

		zlp = cplex.getObjValue(); //1
		
	}

	cplex.end();
	model.end();
	env.end();

	return zlp;
}

double LPFlowSolve(int** cij, int** conij, int** adj, int* d, int* f, vector<std::tuple<int, int, double>>& Zhead, vector<std::tuple<int, int, double>>& Xhead, vector<pair<int, double>>& Yhead, int N, int M, int k)
{
	IloEnv env;
	IloModel model(env);
	NumVarArray2 a = CreateNumVarArray2(env, N, N, "a", 0, IloInfinity);
	IloBoolVarArray s = CreateBoolVarArray(env, N, "s"); //flow source
	BoolVarArray2 x = CreateBoolVarArray2(env, N, N, "x");
	BoolVarArray2 Z = CreateBoolVarArray2(env, N, N, "Z");
	IloBoolVarArray y = CreateBoolVarArray(env, N, "y");

	IloConversion yconversion(env, y, ILOFLOAT);
	model.add(yconversion);

	IloConversion sconversion(env, s, ILOFLOAT);
	model.add(sconversion);

	for (int i = 0; i < N; ++i)
	{
		IloConversion Zconversion(env, Z[i], ILOFLOAT);
		model.add(Zconversion);
		IloConversion xconversion(env, x[i], ILOFLOAT);
		model.add(xconversion);
	}

	IloExpr objExp(env);

	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j)
			objExp += d[i] * cij[i][j] * x[i][j];

	//fixed opening cost
	for (int i = 0; i < N; ++i)
		objExp += f[i] * y[i];

	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j)
			if (adj[i][j] == 1)
				objExp += M *conij[i][j] * Z[i][j];

	IloObjective obj = IloMinimize(env, objExp, "obj");
	model.add(obj);
	objExp.end();

	//Baseline
	// Coverage constraint		
	for (int i = 0; i < N; ++i)
	{
		IloExpr costcover(env);
		for (int j = 0; j < N; ++j)
			costcover += x[i][j];
		//30.1.23
		model.add(costcover == 1);
		costcover.end();
	}

	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
			model.add(x[i][j] <= y[j]);
	}

	//select k facility
	IloExpr facilsum(env);
	for (int i = 0; i < N; ++i) //foreach
		facilsum += y[i];
	model.add(facilsum == k);
	facilsum.end();

	//new-2-min span tree size
	//2-37
	IloExpr z_ijextendedtree(env);
	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j)
			if (adj[i][j] == 1)
				z_ijextendedtree += Z[i][j];

	model.add(z_ijextendedtree == k - 1);
	z_ijextendedtree.end();

	//at most one incoming
	for (int i = 0; i < N; ++i)
	{
		IloExpr z_ijincoming(env);
		for (int j = 0; j < N; ++j)
			if (adj[j][i] == 1)
				z_ijincoming += Z[j][i];

		model.add(z_ijincoming <= y[i]);
		z_ijincoming.end();
	}

	//validineq1
	for (int i = 0; i < N; ++i)
		for (int j = i + 1; j < N; ++j)
			if (adj[i][j] == 1)
				model.add(Z[i][j] + Z[j][i] <= y[i]);


	//validineq2
	for (int i = 0; i < N; ++i)
		for (int j = i + 1; j < N; ++j)
			if (adj[i][j] == 1)
				model.add(Z[i][j] + Z[j][i] <= y[j]);


	//adjusted version wrt bileteral flow
	for (int i = 0; i < N; ++i)
	{
		IloExpr Zijneighborsum(env);
		for (int j = 0; j < N; ++j)
			if (adj[i][j] == 1)
				Zijneighborsum += Z[i][j] + Z[j][i];

		model.add(y[i] <= Zijneighborsum);
		Zijneighborsum.end();
	}

	//valid inequality
	for (int i = 0; i < N; ++i)
	{
		IloExpr yneighborsum(env);
		for (int j = 0; j < N; ++j)
			if (adj[i][j] == 1)
				yneighborsum += y[j];

		model.add(y[i] <= yneighborsum);
		yneighborsum.end();
	}


	//Flow constraints
	//17a
	for (int j = 0; j < N; ++j)
		model.add(s[j] <= y[j]);

	//17b
	IloExpr dummyexpr2(env);
	for (int j = 0; j < N; ++j)
		dummyexpr2 += s[j];
	model.add(dummyexpr2 == 1);
	dummyexpr2.end();

	////11.10.23
	for (int i = 0; i < N; ++i)
	{
		IloExpr inflowexpr(env);
		for (int j = 0; j < N; ++j)
			if (adj[j][i] == 1)
				inflowexpr += a[j][i];

		IloExpr outflowexpr(env);
		for (int j = 0; j < N; ++j)
			if (adj[i][j] == 1)
				outflowexpr += a[i][j];

		model.add(inflowexpr - outflowexpr == y[i] - k*s[i]);
		inflowexpr.end();
		outflowexpr.end();
	}

	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j)
			if (adj[i][j] == 1)
				model.add(a[i][j] <= (k - 1)*Z[i][j]);

	IloCplex cplex(model);
	double cplextotal = 7813;
	cplex.setParam(IloCplex::TiLim, cplextotal);
	double zlp;

	IloBool success = cplex.solve();

	if (success && cplex.isPrimalFeasible())
	{
		for (int i = 0; i < N; ++i)
			for (int j = 0; j < N; ++j)
				Xhead.emplace_back(i, j, cplex.getValue(x[i][j]));

		for (int i = 0; i < N; ++i)
			Yhead.push_back({ i, cplex.getValue(y[i]) });

		for (int i = 0; i < N; ++i)
			for (int j = 0; j < N; ++j)
				if (adj[i][j] == 1)
					Zhead.emplace_back(i, j, cplex.getValue(Z[i][j]));;

		zlp = cplex.getObjValue(); //1
	}

	cplex.end();
	model.end();
	env.end();

	return zlp;

}


int main()
{
	ofstream errorfile("eror.txt");

	try
	{
		int M = 100;
		string computer = "burakpc";
		string url;

		if (computer == "hp")
			url = "C:\\Users\\murat\\Desktop\\sensornetwork\\networkoptim\\input\\";
		else
			url = "C:\\Users\\bkocuk\\Desktop\\networkoptim\\input\\";

		string pathheur = url + "heursettings.txt";
		ifstream heurloop(pathheur);

		double epsilon = 1.0;
		double tau = 0.01;
		double gamma = pow(10, -5);

		ofstream summaryfile("summary.txt");

		if (heurloop.is_open())
		{
			string row;

			while (getline(heurloop, row))
			{
				//changable
				string edgefilt;
				string reduct;
				string verfilt;
				string model;

				double zsum = 0;
				double timesum = 0;
				double resultsum = 0;

				istringstream Heurloop(row);
				Heurloop >> model >> verfilt >> edgefilt >> reduct;
				
				int faclb = 1;
				int facub = 5;

				for (int fac = faclb; fac < facub; ++fac) // 10-20-30-40
				{
					///////////////////////////////////////////////777
					int k = 10 * fac;
					summaryfile << model << "_" << verfilt << "_" << edgefilt << "_" << reduct << "_" << k;
					//begin of Node limit loop
					string iplist2 = "trials-" + to_string(k);
					string iplist = "trials";
					string pathloop = url + iplist2 + ".txt";

					string line;

					ofstream resultfile(model + "_" +
						verfilt + "_" + edgefilt + "_" + reduct + "_" + to_string(k) + "_" + computer + "_.txt");
					ifstream fileloop(pathloop);
					resultfile << "fileID\t" << "preprocess\t" << "Tottime\t" << "Trans\t" << "Con\t" << "Open\t" <<
						"Z_LP\t" << "Z\t" << "gaplb%\t" << "Ec\t" << "y" << endl;

					///////////////////////////////////////////////777
					int count = 0;
					if (fileloop.is_open())
					{
						while (getline(fileloop, line))
						{		
							++count;
							vector<int> V; //7.3.22
							string filename;
							int N;
							istringstream Streamloop(line);
							Streamloop >> filename >> N;

							string path = url + iplist + "\\" + filename + "\\" + filename + ".txt";
							string path2 = url + iplist + "\\" + filename + "\\d.txt";
							string path3 = url + iplist + "\\" + filename + "\\f.txt";

							ifstream file(path);
							string line;

							//resultfile << "\n file: " + filename << endl;

							int** conij;
							int** cij;
							int** adj;
							int** adjcopy;

							int routecost = 100000; //this determines the general objective!

							cij = new int*[N];
							conij = new int*[N];
							adj = new int*[N];
							adjcopy = new int*[N];

							for (int i = 0; i < N; i++)
							{
								cij[i] = new int[N];
								adj[i] = new int[N];
								conij[i] = new int[N];
								adjcopy[i] = new int[N];
							}

							for (int i = 0; i < N; i++)
								for (int j = 0; j < N; j++)
								{
									if (i == j)
									{
										cij[i][j] = 0;
										adj[i][j] = 0;
										adjcopy[i][j] = 0;
										conij[i][j] = 0;
									}

									else
									{
										cij[i][j] = routecost; //arbitrary large number
										adj[i][j] = routecost;
										conij[i][j] = routecost;
										adjcopy[i][j] = routecost;
									}

								}

							ifstream file3(path3);
							ifstream file2(path2);

							int* d = new int[N];
							int* f = new int[N];

							GetData(file, file2, file3, line, d, f, cij, adj);

							//con=cij
							for (int i = 0; i < N; i++)
								for (int j = 0; j < N; j++)
									conij[i][j] = cij[i][j];

							allpairshort(cij, N); //buradaki hatalar c_ij diagonalleri d���nmedi�inden!

							auto startmodel = chrono::steady_clock::now();
							double result;
							//construct S_j
							vector<std::tuple<int, int, double>> Zhead;
							vector<std::tuple<int, int, double>> Xhead;
							vector<std::pair<int, double>> Yhead;

							//Phase 1
							if (model == "MTZ")
								result = LPMTZSolve(cij, conij, adj, d, f, Zhead, Xhead, Yhead, N, M, k);
							else
								result = LPFlowSolve(cij, conij, adj, d, f, Zhead, Xhead, Yhead, N, M, k);
							
							//output: x_ij, y_i, z_ij fractionals
							//////

							double* Chead = new double[N];

							for (int i = 0; i < N; ++i)
								Chead[i] = 0;

							for (const auto& tuple : Xhead)
							{
								int i = get<0>(tuple);
								int j = get<1>(tuple);
								double val = get<2>(tuple);

								//cout << i << " " << j << " " << val << endl;
								Chead[i] += cij[i][j] * val;
							}

							/*for (int i = 0; i < N;++i)
							cout << Chead[i]<< endl;*/

							vector<vector<int>> S;
							S.resize(N);
							vector<int> Fc;
							vector<pair<int, int>> Ec;
							vector<set<int>> components;
							Graph g(N);

							if (verfilt == "LV" || verfilt == "LV-g")
							{
								if (verfilt == "LV-g")
									epsilon = 1 / epsilon;

								for (const auto& pair : Yhead)
								{
									int j = pair.first;
									double val = pair.second;

									if (val <= gamma)
										continue;

									for (int i = 0; i < N; ++i)
										if (cij[i][j] <= (1 + epsilon)*Chead[i])
											S[j].push_back(i);
								}


								if (verfilt == "LV-g")
								{
									vector<vector<int>> deltaS;
									deltaS.resize(N);

									for (int i = 0; i < N; ++i)
										for (int j = i + 1; j < N; ++j)
										{
											bool intersect = FindIntersection(S[i], S[j]);
											if (intersect)
											{
												deltaS[i].push_back(j);
												deltaS[j].push_back(i);
											}
										}

									for (int i = 0; i < N; ++i)
										for (int el : deltaS[i])
											S[i].push_back(el);
								}

								vector<pair<int, double>> indexValuePairs;

								for (int i = 0; i < S.size(); ++i)
									indexValuePairs.push_back({ i, Chead[i] });

								// Sort the vector of pairs based on values in increasing order
								sort(indexValuePairs.begin(), indexValuePairs.end(), [](const pair<int, double>& a, const pair<int, double>& b)
								{
									return a.second < b.second; // Sort in ascending order
								});

								//Step5
								Fc = FindCover(S, indexValuePairs);
							}

							if (verfilt == "Ch")
							{
								vector<pair<int, double>> Cheadvec;
								vector<std::tuple<int, int, int, int, int>> Realdemandpos; // i, d_i, s_i, c_isi, fi
								vector<pair<int, int>> Demandpos;

								for (int i = 0; i < N; ++i)
								{
									Cheadvec.push_back({ i, Chead[i] });
									Demandpos.emplace_back(i, d[i]);
								}

								sort(Cheadvec.begin(), Cheadvec.end(), [](const pair<int, double>& a, const pair<int, double>& b)
								{
									return a.second < b.second; // Sort in ascending order
								});

								//1.11.24
								for (int i = 0; i < N; ++i)
									for (int j = 0; j < i; ++j)
									{
										if (Demandpos[Cheadvec[j].first].second>0 && cij[Cheadvec[i].first][Cheadvec[j].first] <= 4 * Cheadvec[j].second)
										{
											Demandpos[Cheadvec[j].first].second = Demandpos[Cheadvec[i].first].second + Demandpos[Cheadvec[j].first].second;
											Demandpos[Cheadvec[i].first].second = 0;
										}

									}

								//vector<i> Realdemandpos;
								vector<pair<int, int>> Demandposcopy = Demandpos;
								Demandpos.clear();
								for (pair<int, int> mypair : Demandposcopy)
									if (mypair.second > 0)
										Demandpos.emplace_back(mypair.first, mypair.second);

								//vector<std::tuple<int,int,int, int>> Realdemandpos;
								for (pair<int, int> pairi : Demandpos)
								{
									int min = INT_MAX;
									int nearneighbor = -1;
									int first = pairi.first;
									for (pair<int, int> pairj : Demandpos)
									{
										int second = pairj.first;
										if (first == second)
											continue;
										if (cij[first][second] < min)
										{
											nearneighbor = second;
											min = cij[first][second];
										}

									}

									Realdemandpos.emplace_back(first, pairi.second, nearneighbor, min, f[first]); // i, d_i, s_i, c_isi, fi
								}

								int realposize = Realdemandpos.size();
								//vector<int> U;

								if (realposize > k)
								{
									sort(Realdemandpos.begin(), Realdemandpos.end(), [](const std::tuple<int, int, int, int, int>& a, const std::tuple<int, int, int, int, int>& b) {
										int productA = get<1>(a)*get<3>(a)-get<4>(a); //first one is consolidation, second one is to minimize open cost
										int productB = get<1>(b)*get<3>(b)-get<4>(b);
										return productA > productB; // Sort in descending order
									});

									set<int> yone;
									Graph g(realposize);
									//Rounding firststage-solis oba 5
									vector<pair<int, double>> ypos; //index, yhead
									int count = 1;

									//y_=1 first 2k-Demandpos
									//y-=1/2 2(Demandpost-k) locs in Demandpos
									//rest is zero 

									// Define a map to store key-value pairs
									std::map<int, int> myMap;
									map<int, int> reversemap;

									int index = 0;

									for (const auto& tuple : Realdemandpos)
									{
										myMap[get<0>(tuple)] = index;
										reversemap[index] = get<0>(tuple);
										if (count <= 2 * k - realposize)
										{
											ypos.push_back({ get<0>(tuple), 1.0 });
											yone.insert(index);
											++count;
											++index;
											continue;
										}
										ypos.push_back({ get<0>(tuple), 0.5 });
										//g.addEdge(myMap[get<0>(tuple)], myMap[get<2>(tuple)]);
										++count;
										++index;
									}

									for (int i = 0; i < realposize; ++i)
									{
										if (i < 2 * k - realposize)
											continue;

										auto tuple = Realdemandpos[i];
										//int first = get<0>(tuple);
										int si = get<2>(tuple);
										g.addEdge(i, myMap[si]);
									}

									set<int> dominatingSet = g.findDominatingSet(yone); //u is candidate vertex set

									for (const auto& vertex : dominatingSet)
										Fc.push_back(reversemap[vertex]);
								}

								else //N<k
								{
									for (int i = 0; i < realposize; ++i)
									{
										auto tuple = Realdemandpos[i];
										//int first = get<0>(tuple);
										int vertex = get<0>(tuple);
										Fc.push_back(vertex);
									}
								}
							}

							if (edgefilt == "True")
							{
								for (const auto& tuple : Zhead) //ok
								{
									int i = get<0>(tuple);
									int j = get<1>(tuple);
									double val = get<2>(tuple);

									if (val >= tau)
										Ec.push_back({ i, j });
								}

								if (verfilt != "Ch" && verfilt != "LV" && verfilt != "LV-g") //verfilt false
								{
									for (const auto& pair : Ec)
									{
										int u = pair.first;
										int v = pair.second;
										g.addEdge(u, v);
									} //no filter end

									components = g.findConnectedComponents();
									Fc = g.getAllVertices();
								}

								else //verfilt true
								{
									vector<int> Fccopy = Fc;

									//no filter start
									sort(Fc.begin(), Fc.end());
									for (const auto& pair : Ec)
									{
										int u = pair.first;
										int v = pair.second;

										if (binary_search(Fc.begin(), Fc.end(), u))
											if (binary_search(Fc.begin(), Fc.end(), v))
											{
												g.addEdge(u, v);
												Fccopy.erase(std::remove(Fccopy.begin(), Fccopy.end(), u), Fccopy.end());
												Fccopy.erase(std::remove(Fccopy.begin(), Fccopy.end(), v), Fccopy.end());
											}
									} //no filter end

									components = g.findConnectedComponents();
									Ec = g.getAllEdges();

									//add isolated ones!
									for (int node : Fccopy)
									{
										set<int> nodes;
										nodes.insert(node);
										components.push_back(nodes);
									}
								}
							}

							else //edgefilt false
							{
								//add isolated ones!
								for (int node : Fc)
								{
									set<int> nodes;
									nodes.insert(node);
									components.push_back(nodes);
								}
							}

							//Phase 3: Obtain a pseudo feasible sol
							for (int i = 0; i < N; ++i)
								for (int j = 0; j < N; ++j)
									adjcopy[i][j] = adj[i][j];

							allpairshort(adj, N);

							while (components.size() > 1)
							{
								//using a_ij, create upper triangular smallest distance between components
								int totdist = 0;

								//vector<DijkstraResults> paths;
								vector<Path> endpoints;
								vector<vector<int>> matrix(components.size(), vector<int>(components.size()));

								for (int i = 0; i < components.size(); ++i)
									for (int j = i + 1; j < components.size(); ++j)
									{
										int incumbent = INT_MAX;
										int u, v;
										for (int a : components[i])
											for (int b : components[j])
											{
												int shortdist = adj[a][b];
												if (shortdist < incumbent)
												{
													incumbent = shortdist;
													u = a;
													v = b;
												}

											}

										matrix[i][j] = incumbent; //unit distance between i-j
										endpoints.emplace_back(u, v, incumbent, i, j);
										totdist += incumbent;
									}


								//spanning tree between components
								int** inducedgraph = new int*[components.size()];
								for (int i = 0; i < components.size(); ++i)
									inducedgraph[i] = new int[components.size()];

								for (int i = 0; i < components.size(); ++i)
									for (int j = i + 1; j < components.size(); ++j)
										if (matrix[i][j] > 0)
										{
											inducedgraph[i][j] = matrix[i][j];
											inducedgraph[j][i] = matrix[i][j];
										}
										else
										{
											inducedgraph[i][j] = 0;
											inducedgraph[j][i] = 0;
										}

								vector<vector<int>> G = MinSpanTree(inducedgraph, components.size());

								GraphDijkstra gd(N);
								for (int i = 0; i < N; ++i)
									for (int j = 0; j < N; ++j)
										if (adjcopy[i][j] == 1)
											add_edge(i, j, 1, gd); //1 for unit graph, O.W use c_ij

								for (vector<int> edges : G)
									for (Path p : endpoints)
									{
										if (edges.size() == 0)
											continue;
										int compi = edges[0];
										int compj = edges[1];

										if (compi == p.compi && compj == p.compj)
										{
											if (p.val > 1)
											{
												DijkstraResults route = dijkstra_shortest_path(gd, p.i, p.j); //unit shortest distance over real graph

												for (int i = 0; i < route.path.size() - 1; ++i)
												{
													int u = route.path[i];
													int v = route.path[i + 1];
													if (g.DoesExist(u, v) == false)
														g.addEdge(u, v);
												}

											}
											else
												g.addEdge(p.i, p.j);
											//resultdetail << p.i << " " << p.j << endl;
										}

									}

								components = g.findConnectedComponents();
							}


							Ec = g.getAllEdges();
							Fc = g.getAllVertices();

							if (reduct == "True") //reduction is True
							{
								
								//if s_j is empty, remove from G_c 18.7.24
								//vector<int> vertexcopy = Fc;
								for (int vertex : Fc)
								{
									if (S[vertex].size() == 0)
										S[vertex].push_back(vertex);

									//REFLECT in overleaf

									//Fc.erase(std::remove(Fc.begin(), Fc.end(), vertex), Fc.end());
									////remove edges rooted from v
									//vector<pair<int, int>> edgescopy = Ec;
									//for (int j = 0; j < edgescopy.size(); ++j)
									//{
									//	auto edge = edgescopy[j];
									//	if (edge.first == vertex)
									//	{
									//		std::pair<int, int> pairToRemove = { vertex, edge.second };
									//		Ec.erase(std::remove(Ec.begin(), Ec.end(), pairToRemove), Ec.end());
									//	}
									//}
								}

								//REFLECT in overleaf
								FillNoncoveredclients(S, Fc, cij);
								ConnectedSetCover(Ec, Fc, S); //Ren and Zhao's CSC
							}

							//Phase 4
							int bfredges = Ec.size();
							int bfrvertices = Fc.size();

							if (Fc.size() <= k)
							{
								if (Fc.size() < k)
									Repair(Ec, N, Fc, cij, k);

								if (Fc.size() == k)
								{
									//Find MST

									vector<vector<int>> adjacency_matrix;
									adjacency_matrix.resize(k);

									for (int i = 0; i < k; ++i)
										for (int j = 0; j < k; ++j)
											adjacency_matrix[i].push_back(0);


									std::map<int, int> myMap;
									map<int, int> reversemap;

									//myMap[get<0>(tuple)] = index;
									for (int i = 0; i < Fc.size(); ++i)//(pair<int, int> edge : edges)
									{
										myMap[Fc[i]] = i;
										reversemap[i] = Fc[i];
									}

									for (pair<int, int> edge : Ec)
									{
										int i = myMap[edge.first];
										int j = myMap[edge.second];
										adjacency_matrix[i][j] = conij[edge.first][edge.second];
									}

									Ec.clear();

									// Define the Boost graph type (undirected weighted graph)
									using Graph = boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS,
										boost::no_property, boost::property < boost::edge_weight_t, int >> ;

									Graph g(k);

									// Add edges from the adjacency matrix to the graph
									for (int i = 0; i < k; ++i) {
										for (int j = i + 1; j < k; ++j) {
											if (adjacency_matrix[i][j] != 0) {
												boost::add_edge(i, j, adjacency_matrix[i][j], g);
											}
										}
									}

									// Vector to store the MST parent information
									std::vector<boost::graph_traits<Graph>::vertex_descriptor> parent(k);

									// Run Prim's algorithm
									boost::prim_minimum_spanning_tree(g, &parent[0]);

									for (int i = 0; i < k; ++i)
									{
										if (parent[i] != i)
										{ // Exclude the starting node itself
											// std::cout << reversemap[parent[i]] << " - " << reversemap[i] << std::endl;
											Ec.push_back({ reversemap[parent[i]], reversemap[i] });
										}
									}
								}

								auto end = chrono::steady_clock::now();
								auto totdif = end - startmodel;
								double tottime = chrono::duration_cast<chrono::milliseconds>(totdif).count() / 1000;

								//transportation
								int trans = 0;
								for (int i = 0; i < N; ++i)
								{
									int min = INT_MAX;
									for (int j : Fc)
										if (cij[i][j] < min)
											min = cij[i][j];

									trans += d[i] * min;
								}

								int op = 0;
								for (int i : Fc)
									op += f[i];

								int con = 0;

								//connection cost
								for (const auto& pair : Ec)
									con += M*conij[pair.first][pair.second];

								int Z = trans + con + op;
								double gaplp = ((Z - result) / Z) * 100;

								resultfile << filename << "\t" << tottime << "\t" <<
									tottime << "\t" << trans << "\t" << con << "\t"
									<< op << "\t" << result << "\t" << Z << "\t" <<
									gaplp << "\t" << bfredges << "\t" << bfrvertices << endl;

								zsum += Z;
								timesum += tottime;
								resultsum += result;
							}

							else
							{
								vector<double> output;

								auto endprep = chrono::steady_clock::now();
								auto totprepdif = endprep - startmodel;
								double totprep = chrono::duration_cast<chrono::milliseconds>(totprepdif).count() / 1000;

								if (model == "MTZ")
									SolveMTZ(cij, conij, adj, d, f, N, M, k, Ec, Fc, output);
								else
									SolveFlow(cij, conij, adj, d, f, N, M, k, Ec, Fc, output);

								auto end = chrono::steady_clock::now();
								auto totdif = end - startmodel;
								double tottime = chrono::duration_cast<chrono::milliseconds>(totdif).count() / 1000;
								double trans = output[0];
								double con = output[1];
								double op = output[2];
								double Z = output[3];
								double gaplp = ((Z - result) / Z) * 100;

								resultfile << filename << "\t" << totprep << "\t" <<
									tottime << "\t" << trans << "\t" << con << "\t"
									<< op << "\t" << result << "\t" << Z << "\t" <<
									gaplp << "\t" << bfredges << "\t" << bfrvertices << endl;

								zsum += Z;
								timesum += tottime;
								resultsum += result;
							}

							//output.push_back(cplex.getObjValue());
							//garbage collection
							delete[] d;
							delete[] f;
							delete[] Chead;
							for (int i = 0; i < N; ++i)
							{
								delete[] cij[i];
								delete[] conij[i];
								delete[] adj[i];
								delete[] adjcopy[i];
							}

						}

					} //endtrialfileloop

					summaryfile << "\t"  << resultsum / count << "\t" << zsum / count << "\t" << timesum / count << endl;

				} //endforloop

			}
		}

	} //endtry

	catch (IloException& exception)
	{
		cout << exception.getMessage() << endl;
	}

	return 0;
}