clear all; clc;

%% User-specified parameters
rho        = 2.0;
actDens    = 1.0;
Nx         = 5;
numNonzero = Nx;
locality   = 2; % 1 is self-only
tFIR       = 2;
eps = 1e-8;

params = MPCParams();
params.locality_ = locality;
params.tFIR_     = tFIR;
    
sys = LTISystem(); sys.Nx = Nx; 
generate_rand_chain(sys, rho, actDens);
sys.B2 = eye(sys.Nx);

x0 = zeros(sys.Nx, 1);
x0(randsample(sys.Nx, numNonzero)) = rand(numNonzero, 1);
    
[C1, C2, C3] = get_locality_subspace(sys, x0, params);
% C1 = Z2*X, C3 = (I - Fp*F)

rankRatios = rank(C1*C3, eps) / rank(C1, eps)