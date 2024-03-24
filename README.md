# SensorNetworkOptimization
Relevant Files for Connected Gateway Deployment Problem in Sensor Networks

Folder explanations are as follows:

ARRS_Simulation_env: Simulation files for ARRS

codes:

1- Unit graph generation codes
2-"unit graphs to input.cpp"
3- In all files in "parameters" folder, there is a pattern given below:

  3.1. Up to last column: Unit graph generation parameters
  3.2. Last column: Random seed for "unit graphs to input.cpp"

4- Erdos rhenyi

unit graphs: Unit graphs obtained by with the help of unit graph generation codes. Graphs are as follows:

1-albert barabasi
2-benchmark
3-forest fire

graphs: Input folder for our optimization experiments. 

Each folder is comprised of three elements:

1. d-demand for each node
2. f-Deployment cost of each gateway for each node
3. Weighted graph:
 First two column: unit graphs
 Third column: Routing cost of corresponding edge

Algorithm: Turns unit graphs into graphs through "unit graphs to input.cpp"

