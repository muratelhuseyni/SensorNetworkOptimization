import networkx as nx
#import random

#n: the number of nodes 
#m: the number of attachment for the each incoming node to the existing nodes

m=15

file_path = "output.txt"
my_string = "file Node seed edge dens n m\n"

for i in range(1, 5):
    for j in range(1, 6):
#        random.seed(j)
        n = 100 * i
        
        G= nx.barabasi_albert_graph(n,m,seed=j)
        edge_list = G.edges()
        
        size = len(edge_list)
        density = size / (n * n) * 2
        
        result = "ab-n-" + str(n) + "-"+str(j)
        my_string += result + " " + str(n) + " " + str(j) + " " + str(size) + " " + str (density) + " " + str(n) + " " + str (m) + "\n"
        m=m+2
        result +=".txt"
        # Write the edge list to a file in edgelist format
        with open(result, "w") as f:
           for edge in edge_list:
               f.write(str(edge[0]) + " " + str(edge[1]) + "\n")

# Write the string to the text file
with open(file_path, "w") as f:
    f.write(my_string)