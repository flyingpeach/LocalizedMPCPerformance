clear all; clc;
warning off;

%% User-specified parameters
numSimsPerPt   = 5;
tFIRs          = [5, 10, 15, 20, 25, 30];
gridSize       = 5;
actDens        = 1.0; % and 0.8
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
            
    % Use custom communication structure for grid
    systems{j}.AComm = adjust_grid_sys_locality(systems{j}.A);
end

%% Simulations
params          = MPCParams();
numHorizonSizes = length(tFIRs);
locSizes        = zeros(numHorizonSizes, numSimsPerPt);

for i=1:numHorizonSizes
    params.tFIR_ = tFIRs(i);
    fprintf('Simulating horizon size %d of %d\n', i, numHorizonSizes);
    for j=1:numSimsPerPt
        fprintf('\tSim %d of %d\n', j, numSimsPerPt);     
        locSizes(i,j) = get_ideal_locality(systems{j}, params);
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