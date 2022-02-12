clear all; close all; clc;

%% Setup plant + parameters
sys    = LTISystem;
sys.Nx = 10; sys.B1 = eye(sys.Nx);
alpha = 0.8; rho = 2; actDens = 0.7; 
generate_dbl_stoch_chain(sys, rho, actDens, alpha);

params = MPCParams();
params.tFIR_     = 2;
params.locality_ = 2;

x0 = zeros(sys.Nx, 1);
x0(randi(sys.Nx)) = 1;


%% Algorithm
Nx   = sys.Nx; Nu = sys.Nu; T = params.tFIR_;
nPhi = Nx*T + Nu*(T-1);

ZAB  = get_constraint_zab(sys, T);
IO   = [eye(Nx); zeros(Nx*(T-1), Nx)];
ZABp = pinv(ZAB);

Z1 = ZABp*IO;              Z1 = Z1(Nx+1:end, :);
Z2 = eye(nPhi) - ZABp*ZAB; Z2 = Z2(Nx+1:end, :);

PsiSupp = get_sparsity_psi(sys, params);
PsiSupp = PsiSupp(Nx+1:end, 1:Nx);

Z = []; X = [];
for i=1:Nx
    Z = blkdiag(Z, Z2);
    X = [X x0(i)*eye(nPhi)];
end

Ed            = eye((nPhi-Nx)*Nx);
nonZero       = find(PsiSupp);
Ed(nonZero,:) = [];

G = -Ed*vec(Z1); 
F = Ed*Z; Fp = pinv(F); IFF = eye(nPhi*Nx)-Fp*F;

fullRank  = rank(Z2*X)
localRank = rank(Z2*X*IFF)



