% Sanity check for nullspace formulation of vectorized Phi
% We solve a per-iteration MPC problem in two ways: standard formulation,
% and nullspace formulation. The results should be identical.
clear; clc;

%% User-tuned params
rng(123);

sys    = LTISystem();
sys.Nx = 10; sys.B1 = eye(sys.Nx);
alpha = 0.8; rho = 1.5; actDens = 0.7; 
generate_dbl_stoch_chain(sys, rho, actDens, alpha);

numNonzero = 1;
locality   = 3; % 1 is self-only 
tFIR       = 4;
tHorizon   = 20;

%% Setup
params = MPCParams();
params.locality_ = locality;
params.tFIR_     = tFIR;

idxs = randperm(sys.Nx);
x0   = zeros(sys.Nx, 1);
x0(idxs(1:numNonzero)) = 10;

params.QSqrt_ = eye(sys.Nx);
params.RSqrt_ = eye(sys.Nu);

%% Standard formulation using matrix, x0 [Adapted from mpc_centralized]
Nx = sys.Nx; Nu = sys.Nu; T = params.tFIR_;

IO = [eye(Nx); zeros(Nx*(T-1), Nx)];

QSqrt = params.QSqrt_;
RSqrt = params.RSqrt_;

ZAB     = get_constraint_zab(sys, T);
PsiSupp = get_sparsity_psi(sys, params); 
PsiSupp = PsiSupp(:, 1:Nx);
suppSizePsi = sum(sum(PsiSupp));

nPhi = Nx*T + Nu*(T-1);

cvx_begin quiet
cvx_precision best

expression Psi(nPhi, Nx)
variable PsiSuppVals(suppSizePsi, 1)
Psi(PsiSupp) = PsiSuppVals; 

% Set up objective function
obj1 = 0;
for k = 1:T
    kx         = get_range(k, Nx); % rows of Psi representing Psix
    vect       = QSqrt*Psi(kx, 1:Nx)*x0;
    obj1 = obj1 + vect'*vect;
    end
for k = 1:T-1
    ku         = Nx*T + get_range(k, Nu); % rows of Psi representing Psiu
    vect       = RSqrt*Psi(ku, 1:Nx)*x0;
    obj1 = obj1 + vect'*vect;
end

% Constraints
ZAB*Psi == IO; % dynamics constraints
minimize(obj1)
cvx_end

traj1 = Psi*x0;

%% Nullspace formulation
[C1, C2, C3] = get_locality_subspace(sys, x0, params);

cvx_begin quiet
expression traj2(nPhi)
variable w(nPhi*Nx)

yp    = pinv(ZAB)*IO*x0;
yn    = [zeros(Nx,1); C1*C2 + C1*C3*w];
traj2 = yp + yn;

obj2 = 0;
for k = 1:T
    kx         = get_range(k, Nx); % rows of Psi representing Psix
    vect       = QSqrt*traj2(kx);
    obj2 = obj2 + vect'*vect;
    end
for k = 1:T-1
    ku         = Nx*T + get_range(k, Nu); % rows of Psi representing Psiu
    vect       = RSqrt*traj2(ku);
    obj2 = obj2 + vect'*vect;
end

minimize(obj2)
cvx_end

%% These should be very small
trajDiff = norm(traj1 - traj2) / norm(traj1)
objDiff  = norm(obj1 - obj2) / norm(obj1)

