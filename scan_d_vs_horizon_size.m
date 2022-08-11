clear all; clc;

%% User-specified parameters
numSimsPerPt  = 5;
tFIRs         = [4, 8, 12, 16];
gridSize      = 4;
actDens       = 1.0;  % Translates to 50% actuation on nodes, 
                      % since only freq is actuated
connectThresh  = 0.65; % To make sure grid is connected 
Ts             = 0.2;
adjustLocality = true;

seeds = [726, 730, 731, 732, 733];
% 727-729 are ommitted as they produce a non-connected topology

plotTopology = false;

%% Generate and visualize plants
numNodes      = gridSize*gridSize;
numActs       = round(actDens*numNodes);
systems       = cell(numSimsPerPt, 1);

for j=1:numSimsPerPt
    seed = seeds(j);
    rng(seed);

    [adjMtx, nodeCoords, susceptMtx, inertiasInv, dampings] = generate_grid_topology(gridSize, connectThresh, seed);
    
    if plotTopology
        figure(j);
        plot_graph(adjMtx, nodeCoords, 'k');
    end

    actuatedNodes = randsample(numNodes, numActs);
    systems{j}    = generate_grid_plant(actuatedNodes, adjMtx, susceptMtx, inertiasInv, dampings, Ts);
end

%% Simulations
params          = MPCParams();
numHorizonSizes = length(tFIRs);
locSizes        = zeros(numHorizonSizes, numSimsPerPt);

for i=1:numHorizonSizes
    params.tFIR_ = tFIRs(i);
    fprintf('Simulating horizon size %d of %d\n', i, numHorizonSizes);
    for j=1:numSimsPerPt
        fprintf('\tSim %d of %d\n', j, numSimsPerPt)
        locSizes(i,j) = get_ideal_locality(systems{j}, params, adjustLocality);
    end
end

%% Plots
% Remember that our "locality" is actually different from paper
% paper: d=0 means only self communication
% us   : d=1 means only self communication

figure();
plot(tFIRs, mean(locSizes,2) - 1);
xlabel('Horizon lengths');
ylabel('Minimum locality size');