clear all; clc;

%% Grid example
seed          = 420;
gridSize      = 5;
tFIR          = 10;
connectThresh = 0.65;
actDens       = 1.0;
Ts            = 0.2;

paramsGrid       = MPCParams();
paramsGrid.tFIR_ = tFIR;

numNodes      = gridSize * gridSize; 
numActs       = round(actDens*numNodes);
actuatedNodes = randsample(numNodes, numActs);
[adjMtx, nodeCoords, susceptMtx, inertiasInv, dampings] = generate_grid_topology(gridSize, connectThresh, seed);
sysGrid = generate_grid_plant(actuatedNodes, adjMtx, susceptMtx, inertiasInv, dampings, Ts);

adjustLocality = true;
localityGrid   = get_ideal_locality(sysGrid, paramsGrid, adjustLocality)

%plot_graph(adjMtx, nodeCoords, 'k')