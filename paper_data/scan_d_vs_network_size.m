clear all; clc;
warning off;

%% User-specified parameters
numSimsPerPt   = 5;
tFIR           = 15;
gridSizes      = [4, 5, 6, 8, 11];
actDensities   = [1.0, 0.8, 0.6];
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
numActDens   = length(actDensities);
systems      = cell(numActDens, numGridSizes, numSimsPerPt);

% If gridSize >= 11, use connectThresh 0.6
% at that size it's too likely to make an unconnected graph
for i=1:numGridSizes
    gridSize = gridSizes(i);
    numNodes = gridSize*gridSize;

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

        for actIdx = 1:numActDens
            actDens = actDensities(actIdx);
            numActs = round(actDens*numNodes);
            actuatedNodes = randsample(numNodes, numActs);
            systems{actIdx,i,j}  = generate_grid_plant(actuatedNodes, adjMtx, susceptMtx, inertiasInv, dampings, Ts);

            % Use custom communication structure for grid
            systems{actIdx,i,j}.AComm = adjust_grid_sys_locality(systems{actIdx,i,j}.A);
        end
    end
end

%% Simulations
params       = MPCParams();
params.tFIR_ = tFIR;
locSizes     = cell(numActDens, 1);

for actIdx=1:numActDens
    locSizes{actIdx} = zeros(numGridSizes, numSimsPerPt);
    actDens = actDensities(actIdx);
    fprintf('Simulating actuation density %d of %d\n', actIdx, numActDens);
    
    for i=1:numGridSizes
        fprintf('Simulating network size %d of %d\n', i, numGridSizes);
        for j=1:numSimsPerPt
            fprintf('\tSim %d of %d\n', j, numSimsPerPt)
            locSizes{actIdx}(i,j) = get_ideal_locality(systems{actIdx,i,j}, params);
        end
    end
end

save('data/scan_d_vs_network_size.mat');

%% Plots
% Remember that our "locality" is actually different from paper
% paper: d=0 means only self communication
% us   : d=1 means only self communication

% TODO: plot average and standard deviations

load('data/scan_d_vs_network_size.mat');
figure(); hold on;
for actIdx=1:numActDens
    actDens = actDensities(actIdx);
    plot(gridSizes.^2, mean(locSizes{actIdx},2) - 1);
end

legend('100% Actuated', '80% Actuated', '60% Actuated'); % Change this if needed
xlabel('# Subsystems');
ylabel('Minimum locality size');
ylim([0 6]);