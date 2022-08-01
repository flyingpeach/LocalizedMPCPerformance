clear all; clc;

%% User-specified parameters
numSimsPerPt  = 5;
tFIR          = 10;
gridSize      = 4;
specRads      = [1.0, 1.5, 2.0, 4.0];
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
params       = MPCParams();
params.tFIR_ = tFIR;
numSpecRads  = length(specRads);
locSizes     = zeros(numSpecRads, numSimsPerPt);

for i=1:numSpecRads
    fprintf('Simulating spectral radius size %d of %d\n', i, numSpecRads);
    for j=1:numSimsPerPt
        fprintf('\tSim %d of %d\n', j, numSimsPerPt)
        
        sys             = systems{j}.copy(); % deepcopy
        specRadOriginal = max(abs(eig(sys.A)));
        sys.A           = sys.A / specRadOriginal * specRads(i);
        
        locSizes(i,j) = get_ideal_locality(sys, params, adjustLocality);
    end
end