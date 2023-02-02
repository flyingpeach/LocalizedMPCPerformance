clear all; clc;
warning off;

%% User-specified parameters
numSimsPerPt   = 5;
Ts             = [5, 10, 15, 20, 25, 30];
gridSize       = 5;
actDens        = 1.0;
connectThresh  = 0.65;
TsSamp         = 0.2; % Sampling time

% Need to pick seeds such that grid is fully connected
seeds = [700, 703, 704, 705, 706]; % for gridSize = 5
plotTopology = false;

%% Generate and visualize plants
numNodes = gridSize*gridSize;
numActs  = round(actDens*numNodes);
systems  = cell(numSimsPerPt, 1);

for i=1:numSimsPerPt
    seed = seeds(i);
    rng(seed);
    [adjMtx, nodeCoords, susceptMtx, inertiasInv, dampings] = generate_grid_topology(gridSize, connectThresh, seed);

    if plotTopology
        figure(i);
        plot_graph(adjMtx, nodeCoords, 'k');
    end

    actuatedNodes = randsample(numNodes, numActs);
    systems{i}  = generate_grid_plant(actuatedNodes, adjMtx, susceptMtx, inertiasInv, dampings, TsSamp);
        
    % Use custom communication structure for grid
    systems{i}.AComm = adjust_grid_sys_locality(systems{i}.A);
end


%% Simulations
params          = MPCParams();
numHorizonSizes = length(Ts);

locSizes     = zeros(numHorizonSizes, numSimsPerPt);
parTimes     = zeros(numHorizonSizes, numSimsPerPt);
rankTimes    = zeros(numHorizonSizes, numSimsPerPt);

for i=1:numHorizonSizes
    params.tFIR_ = Ts(i) + 1; % Code and paper use different conventions
    fprintf('Simulating horizon size %d of %d\n', i, numHorizonSizes);
    for j=1:numSimsPerPt
        fprintf('\tSim %d of %d\n', j, numSimsPerPt)
        [locSizes(i,j), parTimes(i,j), rankTimes(i,j)] = get_optimal_locality(systems{j}, params);
    end
end

save('data/runtime_vs_horizon_size.mat');

%% Plots
% Remember that our "locality" is actually different from paper
% paper: d=0 means only self communication
% us   : d=1 means only self communication

load('data/runtime_vs_horizon_size.mat');
figure(); hold on;
plot(Ts, log(mean(parTimes,2)));
plot(Ts, log(mean(rankTimes,2)));
legend('Matrix construction (parallelized)', 'Rank determination (non-parallelized)');
xlabel('Horizon length');
ylabel('Time (s)');