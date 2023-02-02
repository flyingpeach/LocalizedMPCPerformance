clear all; clc;
warning off;

%% User-specified parameters
numSimsPerPt  = 5;
T             = 15;gridSize      = 5;
specRads      = [0.5, 1.0, 1.5, 2.0, 2.5];
actDensities   = [1.0, 0.8, 0.6];
connectThresh  = 0.65;
Ts             = 0.2;

% Need to pick seeds such that grid is fully connected
seeds = [700, 703, 704, 705, 706]; % for gridSize = 5
plotTopology = false;

%% Generate and visualize plants
numNodes      = gridSize*gridSize;
numActDens   = length(actDensities);
systems      = cell(numActDens, numSimsPerPt);

for j=1:numSimsPerPt
    seed = seeds(j);
    rng(seed);

    [adjMtx, nodeCoords, susceptMtx, inertiasInv, dampings] = generate_grid_topology(gridSize, connectThresh, seed);
    
    if plotTopology
        figure(j);
        plot_graph(adjMtx, nodeCoords, 'k');
    end

    for actIdx = 1:numActDens
        actDens = actDensities(actIdx);
        numActs = round(actDens*numNodes);
            
        actuatedNodes = randsample(numNodes, numActs);
        systems{actIdx,j}    = generate_grid_plant(actuatedNodes, adjMtx, susceptMtx, inertiasInv, dampings, Ts);
    
        % Use custom communication structure for grid
        systems{actIdx,j}.AComm = adjust_grid_sys_locality(systems{actIdx,j}.A);

    end
end

%% Simulations
params       = MPCParams();
params.tFIR_ = T+1; % Code and paper use different conventions
numSpecRads  = length(specRads);
locSizes     = cell(numActDens, 1);

for actIdx=1:numActDens
    locSizes{actIdx} = zeros(numSpecRads, numSimsPerPt);
    actDens = actDensities(actIdx);
    fprintf('Simulating actuation density %d of %d\n', actIdx, numActDens);
    
    for i=1:numSpecRads
        fprintf('Simulating spectral radius size %d of %d\n', i, numSpecRads);
        for j=1:numSimsPerPt
            fprintf('\tSim %d of %d\n', j, numSimsPerPt)

            sys             = systems{actIdx,j}.copy(); % deepcopy
            specRadOriginal = max(abs(eig(sys.A)));
            sys.A           = sys.A / specRadOriginal * specRads(i);
            locSizes{actIdx}(i,j) = get_optimal_locality(sys, params);
        end
    end
end

save('data/d_vs_spec_radius.mat');

%% Plots
% Remember that our "locality" is actually different from paper
% paper: d=0 means only self communication
% us   : d=1 means only self communication

load('data/d_vs_spec_radius.mat');

figure(); hold on;
for actIdx=1:numActDens
    actDens = actDensities(actIdx);
    plot(specRads, mean(locSizes{actIdx},2) - 1);
end

legend('100% Actuated', '80% Actuated', '60% Actuated'); % Change this if needed
xlabel('Spectral radius');
ylabel('Minimum locality size');
ylim([0 5])