clear all; clc;
warning off;

%% User-specified parameters
numSimsPerPt  = 5;
tFIR          = 15;
gridSize      = 5;
specRads      = [0.5, 1.0, 1.5, 2.0, 2.5];
actDens       = 1.0; % and 0.8
connectThresh  = 0.65;
Ts             = 0.2;

% Need to pick seeds such that grid is fully connected
seeds = [700, 703, 704, 705, 706]; % for gridSize = 5
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
        
        % Use custom communication structure for grid
        sys.AComm = adjust_grid_sys_locality(sys.A);

        locSizes(i,j) = get_ideal_locality(sys, params);
    end
end

%% Plots
% Remember that our "locality" is actually different from paper
% paper: d=0 means only self communication
% us   : d=1 means only self communication

figure();
plot(specRads, mean(locSizes,2) - 1);
xlabel('Spectral radius');
ylabel('Minimum locality size');