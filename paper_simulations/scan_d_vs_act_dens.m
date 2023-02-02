clear all; clc;
warning off;

%% User-specified parameters
numSimsPerPt   = 5;
T              = 15;
gridSize       = 5;
actDensities   = [1.0, 0.8, 0.6, 0.4, 0.2];
connectThresh  = 0.65;
Ts             = 0.2;

% Need to pick seeds such that grid is fully connected
seeds = [700, 703, 704, 705, 706]; % for gridSize = 5

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
        
        % Use custom communication structure for grid
        systems{i,j}.AComm = adjust_grid_sys_locality(systems{i,j}.A);
    end
end

%% Simulations
params       = MPCParams();
params.tFIR_ = T+1; % Code and paper use different conventions
locSizes     = zeros(numActDens, numSimsPerPt);

for i=1:numActDens
    fprintf('Simulating actuation density size %d of %d\n', i, numActDens);
    for j=1:numSimsPerPt
        fprintf('\tSim %d of %d\n', j, numSimsPerPt);
        locSizes(i,j) = get_ideal_locality(systems{i,j}, params);
    end
end

save('data/scan_d_vs_act_dens.mat');

%% Plots
% Remember that our "locality" is actually different from paper
% paper: d=0 means only self communication
% us   : d=1 means only self communication

load('data/scan_d_vs_act_dens.mat');
figure();
plot(actDensities, mean(locSizes,2) - 1);
xlabel('Actuation density (freq only)');
ylabel('Minimum locality size');