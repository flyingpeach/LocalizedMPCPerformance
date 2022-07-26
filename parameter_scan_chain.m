clear all; clc;

%% User-specified parameters
numSims = 200;

rhoMax     = 2.0;
actDensMin = 0.5;
actDensMax = 1.0;

Nx         = 10;
numNonzero = 10;
locality   = 2; % 1 is self-only
tFIR       = 3;

eps = 1e-8;

%% Simulations
params = MPCParams();
params.locality_ = locality;
params.tFIR_     = tFIR;

rankRatios = zeros(numSims, 1);

for i=1:numSims
    if ~mod(i, 10)
        fprintf('Simulating seed %d of %d\n', i, numSims);
    end
    rng(i); % Seeded for reproducibility
    
    rho     = rand() * rhoMax;
    actDens = rand() * (actDensMax - actDensMin) + actDensMin;

    sys = LTISystem(); sys.Nx = Nx;
    generate_rand_chain(sys, rho, actDens);
    
    x0 = zeros(sys.Nx, 1);
    x0(randsample(sys.Nx, numNonzero)) = rand(numNonzero, 1);
    
    [C1, C2, C3] = get_locality_subspace(sys, x0, params);
    rankRatios(i) = rank(C1*C3, eps) / rank(C1, eps);
end

mean(rankRatios)