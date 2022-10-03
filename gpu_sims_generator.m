clear all; clc;
warning off;

%% User-specified parameters
numSims        = 20;
tFIR           = 15;
actDensity     = 1.0;
gridSize       = 11;
connectThresh  = 0.6;
Ts             = 0.2;

% Need to pick seeds such that grid is fully connected
% These work for gridSize = 11, connectThresh=0.6
seeds = [700, 702, 707, 710, 712, ...
         715, 717, 718, 719, 720, ...
	     727, 730, 731, 733, 734, ...
	     735, 738, 740, 741, 742];

plotTopology = false;

% Generate and visualize plants
numNodes    = gridSize*gridSize;
systems     = cell(numSims, 1);

for i=1:numSims
    seed = seeds(i);
    rng(seed);

    [adjMtx, nodeCoords, susceptMtx, inertiasInv, dampings] = generate_grid_topology(gridSize, connectThresh, seed);
    
    if plotTopology
        figure(i);
        plot_graph(adjMtx, nodeCoords, 'k');
    end

    numActs       = round(actDensity*numNodes);
    actuatedNodes = randsample(numNodes, numActs);    
    systems{i}    = generate_grid_plant(actuatedNodes, adjMtx, susceptMtx, inertiasInv, dampings, Ts);
        
    % Use custom communication structure for grid
    systems{i}.AComm = adjust_grid_sys_locality(systems{i}.A);
end

%% Output into python-friendly format
Nx     = 2 * numNodes;
Nu     = numNodes;
As     = cell(numSims, 1);
B2s    = cell(numSims, 1);
AComms = cell(numSims, 1);

for i=1:numSims
    As{i}     = systems{i}.A;
    B2s{i}    = systems{i}.B2;
    AComms{i} = systems{i}.AComm;
end

save('gpu_sims_systems.mat', 'Nx', 'Nu', 'As', 'B2s', 'AComms');

%% Check locality sizes (sanity check)
params       = MPCParams();
params.tFIR_ = tFIR;
locSizes     = zeros(1, numSims);

for i=1:numSims
    fprintf('Checking system %d of %d\n', i, numSims);
    locSizes(i) = get_ideal_locality(systems{i}, params);
end

% Remember that our "locality" is actually different from paper
% paper: d=0 means only self communication
% us   : d=1 means only self communication

locSizes - 1