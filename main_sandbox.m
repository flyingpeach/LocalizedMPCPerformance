clear all; clc;
warning off

%% Grid example
seed          = 420;
gridSize      = 5;
tFIR          = 10;
connectThresh = 0.65;
actDens       = 0.8;
Ts            = 0.2;

params       = MPCParams();
params.tFIR_ = tFIR;

numNodes      = gridSize * gridSize; 
numActs       = round(actDens*numNodes);
actuatedNodes = randsample(numNodes, numActs);
[adjMtx, nodeCoords, susceptMtx, inertiasInv, dampings] = generate_grid_topology(gridSize, connectThresh, seed);
sys = generate_grid_plant(actuatedNodes, adjMtx, susceptMtx, inertiasInv, dampings, Ts);

adjustLocality = true;
%plot_graph(adjMtx, nodeCoords, 'k')

%% Get locality
locality = get_ideal_locality(sys, params, adjustLocality)

%% Setup for control problem
x0     = rand(sys.Nx, 1);

params.locality_ = locality;
params.QSqrt_    = diag(rand(sys.Nx, 1));
params.RSqrt_    = diag(rand(sys.Nu, 1));

Nx = sys.Nx; Nu = sys.Nu; T = params.tFIR_;
nPhi  = Nx*T + Nu*(T-1);
IO    = [eye(Nx); zeros(Nx*(T-1), Nx)];
QSqrt = params.QSqrt_;
RSqrt = params.RSqrt_;
ZAB   = get_constraint_zab(sys, T);

%% Calculate global trajectory
fprintf('Calculating global trajectory\n');

cvx_begin quiet
variable Psi1(nPhi, Nx)

% Set up objective function
obj1 = 0;
for k = 1:T
    kx         = get_range(k, Nx); % rows of Psi representing Psix
    vect       = QSqrt*Psi1(kx, 1:Nx)*x0;
    obj1 = obj1 + vect'*vect;
end
for k = 1:T-1
    ku         = Nx*T + get_range(k, Nu); % rows of Psi representing Psiu
    vect       = RSqrt*Psi1(ku, 1:Nx)*x0;
    obj1 = obj1 + vect'*vect;
end

% Constraints
ZAB*Psi1 == IO; % dynamics constraints
minimize(obj1)
cvx_end

%% Calculate localized trajectory [Adapted from mpc_centralized]
fprintf('Calculating localized trajectory\n');

PsiSupp = get_sparsity_psi(sys, params, adjustLocality); 
PsiSupp = PsiSupp(:, 1:Nx);
suppSizePsi = sum(sum(PsiSupp));

cvx_begin quiet
expression Psi2(nPhi, Nx)
variable PsiSuppVals(suppSizePsi, 1)
Psi2(PsiSupp) = PsiSuppVals; 

% Set up objective function
obj2 = 0;
for k = 1:T
    kx         = get_range(k, Nx); % rows of Psi representing Psix
    vect       = QSqrt*Psi2(kx, 1:Nx)*x0;
    obj2 = obj2 + vect'*vect;
end
for k = 1:T-1
    ku         = Nx*T + get_range(k, Nu); % rows of Psi representing Psiu
    vect       = RSqrt*Psi2(ku, 1:Nx)*x0;
    obj2 = obj2 + vect'*vect;
end

% Constraints
ZAB*Psi2 == IO; % dynamics constraints
minimize(obj2)
cvx_end

objDiff  = norm(obj2 - obj1) / obj2
trajDiff = norm(Psi2*x0 - Psi1*x0) / norm(Psi1*x0)