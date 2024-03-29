import networkx as nx
#from random import randint

# Parameters for LFR_benchmark_graph
n = 400
tau1 = 3
tau2 = 2
mu = 0.4
avgdeg = 35 
mincom = 0.20*n
#avgdeg=15, n=100, mincom=0.1n,sed=1
#avgdeg=15, n=200, mincom=0.15n,sed=2
#avgdeg=25, n=300, mincom=0.20n,sed=2
#avgdeg=35, n=400, mincom=0.20n,sed=2

sed = 12

file_path = "benchmark_out.txt"
my_string = "file Node seed edge dens tau1 tau2 mu avgdeg mincom \n"

# Generate the graph
G = nx.generators.community.LFR_benchmark_graph(n, tau1, tau2, mu, average_degree=avgdeg, min_community= mincom, seed=sed)
edge_list = G.edges()

size = len(edge_list)
density = size / (n * n) * 2

result = "bench-" + str(n) + "-"+str(sed)
my_string += result + " " + str(n) + " " + str(sed) + " " + str(size) + " " + str(density) + " " + str(tau1) + " " + str(tau2) + " " + str(mu) + " " + str(avgdeg) + " " + str(mincom) + "\n"

result +=".txt"
# Write the edge list to a file in edgelist format
with open(result, "w") as f:
   for edge in edge_list:
       f.write(str(edge[0]) + " " + str(edge[1]) + "\n")

comsize = len({frozenset(G.nodes[v]["community"]) for v in G})

# Write the string to the text file
with open(file_path, "w") as f:
    f.write(my_string)