clear all; clc;
warning off;

%% User-specified parameters
numSimsPerPt   = 5;
T              = 15;
gridSizes      = [4, 5, 6, 8, 11];
actDens        = 1.0;
connectThresh  = 0.65;
Ts             = 0.2;

% Need to pick seeds such that grid is fully connected
seeds = {[700, 701, 702, 703, 704], % gridSize = 4
         [700, 703, 704, 705, 706], % gridSize = 5
         [700, 702, 706, 707, 709], % gridSize = 6
         [700, 705, 707, 709, 713], % gridSize = 8
         [700, 702, 707, 710, 712]  % gridSize = 11, connectThresh=0.6
        };

plotTopology = false;

%% Generate and visualize plants
numGridSizes = length(gridSizes);
systems      = cell(numGridSizes, numSimsPerPt);

% If gridSize >= 11, use connectThresh 0.6
% at that size it's too likely to make an unconnected graph
for i=1:numGridSizes
    gridSize = gridSizes(i);
    numNodes = gridSize*gridSize;
    numActs  = round(actDens*numNodes);

    ct = connectThresh;
    if gridSize >= 11
        ct = 0.6;
    end
    for j=1:numSimsPerPt
        seed = seeds{i}(j);
        rng(seed);
        [adjMtx, nodeCoords, susceptMtx, inertiasInv, dampings] = generate_grid_topology(gridSize, ct, seed);

        if plotTopology
            figure((i-1)*numSimsPerPt+j);
            plot_graph(adjMtx, nodeCoords, 'k');
        end

        actuatedNodes = randsample(numNodes, numActs);
        systems{i,j}  = generate_grid_plant(actuatedNodes, adjMtx, susceptMtx, inertiasInv, dampings, Ts);
        
        % Use custom communication structure for grid
        systems{i,j}.AComm = adjust_grid_sys_locality(systems{i,j}.A);
    end
end

%% Simulations
params       = MPCParams();
params.tFIR_ = T+1; % Code and paper use different conventions
locSizes     = zeros(numGridSizes, numSimsPerPt);
parTimes     = zeros(numGridSizes, numSimsPerPt);
rankTimes    = zeros(numGridSizes, numSimsPerPt);

for i=1:numGridSizes
    fprintf('Simulating network size %d of %d\n', i, numGridSizes);
    for j=1:numSimsPerPt
        fprintf('\tSim %d of %d\n', j, numSimsPerPt)
        [locSizes(i,j), parTimes(i,j), rankTimes(i,j)] = get_optimal_locality(systems{i,j}, params);
    end
end

save('data/runtime_vs_network_size.mat');

%% Plots
% Remember that our "locality" is actually different from paper
% paper: d=0 means only self communication
% us   : d=1 means only self communication

load('data/runtime_vs_network_size.mat');
figure(); hold on;
plot(gridSizes.^2, log(mean(parTimes,2)));
plot(gridSizes.^2, log(mean(rankTimes,2)));
legend('Matrix construction (parallelized)', 'Rank determination (non-parallelized)');
xlabel('# Subsystems');
ylabel('Time (s)');