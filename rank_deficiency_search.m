clear all; clc;
warning off;

%% Search for rank deficiency
numSims = 2000;

gridSizeMin = 4;
gridSizeMax = 11;

actDensMin = 0.2;
actDensMax = 1.0;

specRadMin = 0.5;
specRadMax = 2.5;

horizonMin = 1;  % Corresponds to T in paper
horizonMax = 20;

Ts = 0.2;
rng(2023);

%% Generate plants
seeds     = randi([1, numSims*10], 1, numSims); % For plant generation
gridSizes = randi([gridSizeMin, gridSizeMax], 1, numSims);
actDens   = rand(1, numSims) * (actDensMax-actDensMin) + actDensMin;
systems   = cell(numSims);

skipSims = []; % Skip unconnected systems
for j=1:numSims
    seed     = seeds(j);
    gridSize = gridSizes(j);
    numNodes = gridSize*gridSize;
    numActs  = round(actDens(j)*numNodes);
    
    connectThresh = 0.65;
    if gridSize > 6 % Larger grids more likely to be disconnected; avoid
        connectThresh = 0.6;
    end

    [adjMtx, nodeCoords, susceptMtx, inertiasInv, dampings] = generate_grid_topology(gridSize, connectThresh, seed);
        
    actuatedNodes = randsample(numNodes, numActs);
    systems{j}    = generate_grid_plant(actuatedNodes, adjMtx, susceptMtx, inertiasInv, dampings, Ts);

    % Use custom communication structure for grid
    systems{j}.AComm = adjust_grid_sys_locality(systems{j}.A);
    
    connected = check_connected(systems{j});
    if ~connected
        %figure(j) % Uncomment to double-check
        %plot_graph(adjMtx, nodeCoords, 'k');
        skipSims = [skipSims j];
    end
end

trueSims    = setdiff(1:numSims, skipSims);
trueNumSims = length(trueSims);

%% Simulations
specRads  = rand(1, numSims) * (specRadMax-specRadMin) + specRadMin;
horizons  = randi([horizonMin, horizonMax], 1, numSims);
locSizes  = zeros(1, numSims);
rankDefs  = false(1, numSims);

for i=1:trueNumSims
    j = trueSims(i);
    fprintf('Simulation %d of %d\n', i, trueNumSims);

    params = MPCParams();
    params.tFIR_ = horizons(j) + 1; % Code and paper use different conventions
    
    sys             = systems{j}.copy(); % deepcopy
    specRadOriginal = max(abs(eig(sys.A)));
    sys.A           = sys.A / specRadOriginal * specRads(j);
    
    [locSizes(j), ~, ~, rankDefs(j)] = get_optimal_locality(sys, params);
    
    if rankDefs(j)
        fprintf('==Rank deficiency found at sim %d\n==', j);
    end
end

save('data/rank_deficiency_search.mat');

%% See if any rank deficiency was reported
find(rankDefs)
