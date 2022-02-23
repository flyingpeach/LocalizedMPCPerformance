clear all; close all; clc;

%% User-tuned params
rng(420);

sys    = LTISystem();
sys.Nx = 8; sys.B1 = eye(sys.Nx);
alpha = 0.8; rho = 1.5; actDens = 0.7; 
generate_dbl_stoch_chain(sys, rho, actDens, alpha);

numNonzero = 1;
locality   = 4; % 1 is self-only 
tFIR       = 4;
tHorizon   = 20;

eps = 1e-8;

%% Setup
params = MPCParams();
params.locality_ = locality;
params.tFIR_     = tFIR;

idxs = randperm(sys.Nx);
x0   = zeros(sys.Nx, 1);
x0(idxs(1:numNonzero)) = 1;
w = zeros(sys.Nx, tHorizon); % noiseless

QSqrt = eye(sys.Nx);
RSqrt = eye(sys.Nu);

%% Space size calculations
[C1, C2, C3] = get_locality_subspace(sys, x0, params);

% Rank ratio of 1 indicates Yd = Y
rankRatio = rank(C1*C3, eps) / rank(C1, eps)

%% MPC Comparison
params.mode_  = MPCMode.Centralized;
params.tFIR_  = tFIR;
params.QSqrt_ = QSqrt;
params.RSqrt_ = RSqrt;

% Note: There's an implicit LB enforced by Gurobi which can cause
%       infeasibility for very unstable systems

[xLoc, uLoc, ~]   = sls_mpc(sys, x0, w, params, tHorizon);

params.locality_  = sys.Nx; % effectively no locality
[xGlob, uGlob, ~] = sls_mpc(sys, x0, w, params, tHorizon);

objGlob = get_cost_fn(params, xGlob, uGlob);
objLoc  = get_cost_fn(params, xLoc, uLoc);

costRatio = objLoc/objGlob


