clear all; clc;
warning off;

%% Search for rank deficiency
numSims = 1000;

gridSizeMin = 4;
gridSizeMax = 11;

actDensMin = 0.2;
actDensMax = 1.0;

specRadMin = 0.5;
specRadMax = 2.5;

horizonMin = 3;  % Corresponds to T in paper
horizonMax = 20;

Ts = 0.2;
rng(2022);

%% Generate plants
seeds     = randi([1, numSims*10], 1, numSims); % For plant generation
gridSizes = randi([gridSizeMin, gridSizeMax], 1, numSims);
actDens   = rand(1, numSims) * (actDensMax-actDensMin) + actDensMin;
systems   = cell(numSims);

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
end

%% Simulations
specRads  = rand(1, numSims) * (specRadMax-specRadMin) + specRadMin;
horizons  = randi([horizonMin, horizonMax], 1, numSims);
locSizes  = zeros(1, numSims);
rankDefs  = false(1, numSims);

for j=1:numSims
    fprintf('Simulation %d of %d\n', j, numSims);

    params = MPCParams();
    params.tFIR_ = horizons(j) + 1; % Code and paper use different conventions
    
    sys             = systems{j}.copy(); % deepcopy
    specRadOriginal = max(abs(eig(sys.A)));
    sys.A           = sys.A / specRadOriginal * specRads(j);
    
    % Scale tolerance with system/horizon size (effective size)
    % Will range from 1e-8 to 1e-5
    effSize = sys.Nx*(horizons(j)) + sys.Nu*(horizons(j)-1);        
    eps = 10^(2*log10(effSize) - 12);
    
    [locSizes(j), ~, ~, rankDefs(j)] = get_optimal_locality(sys, params, eps);
    
    if rankDefs(j)
        fprintf('==Rank deficiency found at sim %d\n==', j);
    end
end

save('data/rank_deficiency_search.mat');

%% See if any rank deficiency was reported
find(rankDefs)
