clear all; clc;

%% User-specified parameters
numSims = 100;

actDensMin = 0.7;
actDensMax = 1.0;

gridSize   = 3;
numNonzero = 2*gridSize*gridSize;
locality   = 2; % 1 is self-only
tFIR       = 3;

connectThresh = 0.6;
eps = 1e-8;

%% Simulations
numNodes = gridSize * gridSize; Ts = 0.2;

params = MPCParams();
params.locality_ = locality;
params.tFIR_     = tFIR;

rankRatios = zeros(numSims, 1);

for i=1:numSims
    if ~mod(i, 10)
        fprintf('Simulating seed %d of %d\n', i, numSims);
    end
    rng(i); % Seeded for reproducibility
    
    actDens = rand() * (actDensMax - actDensMin) + actDensMin;
    numActs        = round(actDens*numNodes);
    actuatedNodes  = randsample(numNodes, numActs);

    [adjMtx, nodeCoords, susceptMtx, inertiasInv, dampings] = generate_grid_topology(gridSize, connectThresh, i);
    sys = generate_grid_plant(actuatedNodes, adjMtx, susceptMtx, inertiasInv, dampings, Ts);
    
    x0 = zeros(sys.Nx, 1);
    x0(randsample(sys.Nx, numNonzero)) = rand(numNonzero, 1);

    [C1, C2, C3] = get_locality_subspace(sys, x0, params);
    rankRatios(i) = rank(C1*C3, eps) / rank(C1, eps);
end

mean(rankRatios)