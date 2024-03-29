#source: leskovec et al- Graphs over Time: Densification Laws, Shrinking
#Diameters and Possible Explanations

#two parameters
#forward burning probability p
#backward burning probability pb:     p * bw.factor(backward burning ratio)

#forest fire R: https://search.r-project.org/CRAN/refmans/igraph/html/sample_forestfire.html
#package install: conda->environment->search igraph->select and install
#package:igraph : select from the R library to run

#leskovec et al. at Fig 5, see dense param: p=0.38, pb = 0.35

fw.prob = 0.50
pb = 0.40 
bw.factor = pb / fw.prob
#myseed=5
# Specify the file path
file_path <- "output.txt"

my_string <- paste("file" , "Node ","seed ","edge ","dens ","fwprob ","pb\n")

ubi=4
ubj=5

for (i in 1:ubi)
  for (j in 1:ubj) 
  {
    set.seed(j)
    n=100*i
    g <- sample_forestfire(n, fw.prob, bw.factor) #bwfactor = p_b/p
    edge_list <- as_edgelist(g)
    size <- length(edge_list)/2
    density <- size/(n*n)*2 
    # Specify the file path
    result <- paste("ff-n-", n, "-", j)
    
    my_string <- paste(my_string,paste(result, " ", n, " ", j, " ", size, " ", density, " ", fw.prob, " ", pb, "\n"))
    
    result <- paste(result,".txt")
    # Write the graph to a file in edgelist format
    write_graph(g, file = result, format = "edgelist")
  }

# Write the string to the text file
writeLines(my_string, file_path)