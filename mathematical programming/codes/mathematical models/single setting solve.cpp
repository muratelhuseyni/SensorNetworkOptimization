#include "Common.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <stack>
#include <chrono>
#include <unordered_set>
#include<array> 

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

int pairshort(int** a, int i, int j, int n)
{
	//as k increases, we add the (0,k-1) ones to permenant
	//at each step, if shortest up to i-k and k-j is less than i-j incumbent, make the new value incumbent

	//subset of allpairs shortestpath
	for (k = 0; k < n; k++)
	{
		if (a[i][k] + a[k][j] < a[i][j])
			a[i][j] = a[i][k] + a[k][j];	
	}

	return a[i][j];
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

ILOLAZYCONSTRAINTCALLBACK5(BendersLazyCallback, BoolVarArray2&, Z, IloBoolVarArray&, y,  int**, adj, int, N, int, k)
{
	vector<int> vertices;
	vector<vector<int>> edges;

	vector<vector<vector<int>>> cycles;

	for (int i = 0; i<N; ++i)
		for (int j = 0; j<N; ++j )
			if (adj[i][j] == 1 && getValue(Z[i][j]) > 0.5)
				edges.push_back({ i, j });

	for (int i = 0; i < N; ++i)
	{
		if (getValue(y[i])>0.5)
			vertices.push_back(i);
	}

	bool found = false;
	for (int w = 3; w < k; ++w)
	{
		vector<vector<int>> subsets = get_k_element_subsets(vertices, w);

		for (int i = 0; i < subsets.size(); ++i)
		{
			vector<int> subset = subsets[i];
			//travel in edges
			int* countdegree = new int[N];
			for (int i = 0; i < N; ++i)
				countdegree[i] = 0;

			for (vector<int> edge : edges)
				for (int vertex1 : subset)
					for (int vertex2 : subset)
						if (vertex1 == edge[0] && vertex2 == edge[1])
						{
							++countdegree[vertex1];
							++countdegree[vertex2];
						}


			//count degrees in the subgraph
			int memberwith2 = 0;
			for (int subsetvertex : subset)
				if (countdegree[subsetvertex] == 2)
					++memberwith2;

			delete countdegree;

			vector<vector<int>> cycle;

			if (memberwith2 == w) //subset size, then cycle found
			{
				for (int i = 0; i < w; ++i)
				{
					int first = subset[i];
					for (int j = 0; j < w; ++j)
					{
						int second = subset[j];
						int cand[2] = { i, j };
						for (vector<int> edge : edges)
							if (first == edge[0] && second == edge[1])
								cycle.push_back(edge);
					}
				}

				cycles.push_back(cycle);
			}
			else
				cycle.clear();
		}

	}

	if (!cycles.empty())
	{
		for (vector<vector<int>> cycle : cycles)
		{
			IloExpr z_ijsol(getEnv());
			//modeladd
			for (vector<int> edge : cycle)
				z_ijsol += Z[edge[0]][edge[1]];

			int cyclesize = cycle.size();
			add(z_ijsol <= cyclesize - 1); //to add the cut on the fly!
			z_ijsol.end();
		}

	}
}

void SolveLazy(int** cij, int** conij, int** adj, int* d, int* f, int N, int M, int k, vector<double>& output)
{
	IloEnv env;
	IloModel model(env);
	BoolVarArray2 x = CreateBoolVarArray2(env, N, N, "x");
	BoolVarArray2 Z = CreateBoolVarArray2(env, N, N, "Z");
	IloBoolVarArray y = CreateBoolVarArray(env, N, "y");
	IloNumVarArray u = CreateNumVarArray(env, N, "u", 0, IloInfinity);

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

	IloCplex cplex(model);
	double cplextotal = 60 * 60 * 2;
	cplex.setParam(IloCplex::TiLim, cplextotal);
	cplex.setParam(IloCplex::Threads, 32);

	cplex.use(BendersLazyCallback(env, Z, y, adj, N, k));

	IloBool success = cplex.solve();

	if (success && cplex.isPrimalFeasible())
	{
		//transportation cost
		double trans = 0;
		double con = 0;
		double op = 0;

		//transport cost
		for (int i = 0; i < N; ++i)
			for (int j = 0; j < N; ++j)
				trans += d[i] * cij[i][j] * cplex.getValue(x[i][j]);

		//opening cost
		for (int i = 0; i < N; ++i)
			op += f[i] * cplex.getValue(y[i]);

		//connection cost
		for (int i = 0; i < N; ++i)
			for (int j = 0; j < N; ++j)
				if (adj[i][j] == 1)
					con += M *conij[i][j] * cplex.getValue(Z[i][j]);

		output.push_back(trans);
		output.push_back(con);
		output.push_back(op);
		output.push_back(cplex.getObjValue());
		output.push_back(cplex.getMIPRelativeGap());

		//resultfile << "fileID\t" << "Trans\t" << "Con\t" << "Open\t" << "Z\t" << "Gap%\t" << "Time" << endl;
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


void SolveMTZ(int** cij, int** conij, int** adj, int* d, int* f, int N, int M, int k, vector<double>& output)
{
	IloEnv env;
	IloModel model(env);
	BoolVarArray2 x = CreateBoolVarArray2(env, N, N, "x");
	BoolVarArray2 Z = CreateBoolVarArray2(env, N, N, "Z");
	IloBoolVarArray y = CreateBoolVarArray(env, N, "y");
	IloNumVarArray u = CreateNumVarArray(env, N, "u", 0, IloInfinity);

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
			for (int j = 0; j < N; ++j)
				trans += d[i] * cij[i][j] * cplex.getValue(x[i][j]);

		//opening cost
		for (int i = 0; i < N; ++i)
			op += f[i] * cplex.getValue(y[i]);

		//connection cost
		for (int i = 0; i < N; ++i)
			for (int j = 0; j < N; ++j)
				if (adj[i][j] == 1)
					con += M *conij[i][j] * cplex.getValue(Z[i][j]);

		output.push_back(trans);
		output.push_back(con);
		output.push_back(op);
		output.push_back(cplex.getObjValue());
		output.push_back(cplex.getMIPRelativeGap());

		//resultfile << "fileID\t" << "Trans\t" << "Con\t" << "Open\t" << "Z\t" << "Gap%\t" << "Time" << endl;
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

void SolveFlow(int** cij, int** conij, int** adj, int* d, int* f, int N, int M, int k, vector<double>& output)
{
	IloEnv env;
	IloModel model(env);
	NumVarArray2 a = CreateNumVarArray2(env, N, N, "a", 0, IloInfinity);
	IloBoolVarArray s = CreateBoolVarArray(env, N, "s"); //flow source
	BoolVarArray2 x = CreateBoolVarArray2(env, N, N, "x");
	BoolVarArray2 Z = CreateBoolVarArray2(env, N, N, "Z");
	IloBoolVarArray y = CreateBoolVarArray(env, N, "y");

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
	double cplextotal = 60*60*2;
	cplex.setParam(IloCplex::TiLim, cplextotal);

	IloBool success = cplex.solve();

	//resultfile << "fileID\t" << "gap%\t" << "time\t" << "Z_feas\t" << endl;

	if (success && cplex.isPrimalFeasible())
	{
		//transportation cost
		double trans = 0;
		double con = 0;
		double op = 0;

		//transport cost
		for (int i = 0; i < N; ++i)
			for (int j = 0; j < N; ++j)
				trans += d[i] * cij[i][j] * cplex.getValue(x[i][j]);

		//opening cost
		for (int i = 0; i < N; ++i)
			op += f[i] * cplex.getValue(y[i]);

		//connection cost
		for (int i = 0; i < N; ++i)
			for (int j = 0; j < N; ++j)
				if (adj[i][j] == 1)
					con += M *conij[i][j] * cplex.getValue(Z[i][j]);

		output.push_back(trans);
		output.push_back(con);
		output.push_back(op);
		output.push_back(cplex.getObjValue());
		output.push_back(cplex.getMIPRelativeGap());

		//resultfile << "fileID\t" << "Trans\t" << "Con\t" << "Open\t" << "Z\t" << "Gap%\t" << "Time" << endl;
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

		//select an experiment path
		int faclb = 3;
		int facub = 5;

		int metlb = 0;
		int metub = 2;

		for (int fac = faclb; fac < facub; ++fac) // 10-20-30-40
		for (int met = metlb; met < metub; ++met) //3 setting: 0-mtz, 1-flow, 2-lazy
		{

			///////////////////////////////////////////////777
			int k = 10 * fac;
			//begin of Node limit loop
			string iplist2 = "trials-" + to_string(k);
			string iplist = "trials";

			string pathloop = url + iplist2 + ".txt";
			 
			//string pathloop = url + iplist + ".txt";
			//string modeltype = "mtz";
			//string root = "Rootnode";
			string method;

			string line;
			//int k = 10;
			//int k = 10;
			int N;

			switch (met)
			{
				case 0:
					method = "MTZ";
					break;
				case 1:
					method = "flow";
					break;
				case 2:
					method = "lazy";
					break;
				default:
					method = "Unknown method for ip";
					break;
			}

			ofstream resultfile(method + "_" + to_string(k) + "_" + computer + "_.txt");
			ifstream fileloop(pathloop);

			resultfile << "fileID\t" << "Trans\t" << "Con\t" << "Open\t" << "Z\t" << "Gap%\t" << "Time" << endl;

			if (fileloop.is_open())
			{
				while (getline(fileloop, line))
				{
					//tsp225
					vector<int> V; //7.3.22
					string filename;
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

					int routecost = 100000; //this determines the general objective!
					
					cij = new int*[N];
					conij = new int*[N];
					adj = new int*[N];

					for (int i = 0; i < N; i++)
					{
						cij[i] = new int[N];
						adj[i] = new int[N];
						conij[i] = new int[N];
					}

					for (int i = 0; i < N; i++)
						for (int j = 0; j < N; j++)
						{
							if (i == j)
							{
								cij[i][j] = 0;
								adj[i][j] = 0;
								conij[i][j] = 0;
							}

							else
							{
								cij[i][j] = routecost; //arbitrary large number
								adj[i][j] = routecost;
								conij[i][j] = routecost;
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

					allpairshort(cij, N); //buradaki hatalar c_ij diagonalleri düţünmediđinden!

					auto startmodel = chrono::steady_clock::now();
					double result;
					vector<double> output;

					switch (met)
					{
						case 0:
							SolveMTZ(cij, conij, adj, d, f, N, M, k, output);
							break;
						case 1:
							SolveFlow(cij, conij, adj, d, f, N, M, k, output);
							break;
						case 2:
							SolveLazy(cij, conij, adj, d, f, N, M, k, output);
							break;
						default:
							method = "Unknown method for ip";
							break;
					}

					auto end2 = chrono::steady_clock::now();
					auto diff2 = end2 - startmodel;

					resultfile << filename << "\t" << output[0] << "\t" << output[1] << "\t" << output[2] << "\t" << output[3] << "\t" <<
						output[4] << "\t" << chrono::duration_cast<chrono::milliseconds>(diff2).count() / 1000 << "\t"
						<< endl;

					//garbage collection
					delete[] d;
					delete[] f;
					for (int i = 0; i < N; ++i)
					{
						delete[] cij[i];
						delete[] conij[i];
						delete[] adj[i];
					}

				}

			   } //endtrialfileloop
	   } //endmethodfor
	} //endtry

	catch (IloException& exception)
	{
		cout << exception.getMessage() << endl;
	}

	return 0;
}