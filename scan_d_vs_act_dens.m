clear all; clc;

%% User-specified parameters
numSimsPerPt  = 5;
tFIR          = 10;
gridSize      = 4;
actDensities  = [0.9, 0.8, 0.7, 0.5];
connectThresh = 0.65; % To make sure grid is connected 
Ts            = 0.2;

seeds = [726, 730, 731, 732, 733];
% 727-729 are ommitted as they produce a non-connected topology

plotTopology = false;

%% Generate and visualize plants
numActDens  = length(actDensities);
numNodes    = gridSize*gridSize;
systems     = cell(numActDens, numSimsPerPt);

for j=1:numSimsPerPt
    seed = seeds(j);
    rng(seed);

    [adjMtx, nodeCoords, susceptMtx, inertiasInv, dampings] = generate_grid_topology(gridSize, connectThresh, seed);
    
    if plotTopology
        figure(j);
        plot_graph(adjMtx, nodeCoords, 'k');
    end

    for i=1:numActDens
        numActs       = round(actDensities(i)*numNodes);
        actuatedNodes = randsample(numNodes, numActs);    
        systems{i,j} = generate_grid_plant(actuatedNodes, adjMtx, susceptMtx, inertiasInv, dampings, Ts);
    end
end

%% Simulations
params       = MPCParams();
params.tFIR_ = tFIR;
locSizes     = zeros(numActDens, numSimsPerPt);

for i=1:numActDens
    fprintf('Simulating actuation density size %d of %d\n', i, numActDens);
    for j=1:numSimsPerPt
        fprintf('\tSim %d of %d\n', j, numSimsPerPt);
        locSizes(i,j) = get_ideal_locality(systems{i,j}, params);
    end
end

%% Plots
% Remember that our "locality" is actually different from paper
% paper: d=0 means only self communication
% us   : d=1 means only self communication

figure();
plot(actDensities, mean(locSizes,2) - 1);
xlabel('Actuation density (freq only)');
ylabel('Minimum locality size');