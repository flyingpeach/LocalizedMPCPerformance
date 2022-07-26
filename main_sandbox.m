clear all; clc;

%% Chain example
rho        = 2.0;
actDens    = 0.5;
Nx         = 15;
tFIR       = 8;

paramsChain       = MPCParams();
paramsChain.tFIR_ = tFIR;
    
sysChain = LTISystem(); sysChain.Nx = Nx; 
generate_rand_chain(sysChain, rho, actDens);

localityChain = get_ideal_locality(sysChain, paramsChain)

%% Grid example
seed          = 420;
gridSize      = 4;
tFIR          = 5;
connectThresh = 0.65;
actDens       = 0.6;
Ts            = 0.2;

paramsGrid       = MPCParams();
paramsGrid.tFIR_ = tFIR;

numNodes      = gridSize * gridSize; 
numActs       = round(actDens*numNodes);
actuatedNodes = randsample(numNodes, numActs);
[adjMtx, nodeCoords, susceptMtx, inertiasInv, dampings] = generate_grid_topology(gridSize, connectThresh, seed);
sysGrid = generate_grid_plant(actuatedNodes, adjMtx, susceptMtx, inertiasInv, dampings, Ts);

localityGrid = get_ideal_locality(sysGrid, paramsGrid)

plot_graph(adjMtx, nodeCoords, 'k')