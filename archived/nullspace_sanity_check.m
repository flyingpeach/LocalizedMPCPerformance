% Sanity check for nullspace formulation of vectorized Phi
% We solve a per-iteration MPC problem in two ways: standard formulation,
% and nullspace formulation. The results should be identical.
clear; clc;

%% User-tuned params
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

%% Nullspace formulation (first version)
Zp = ZAB\IO;
Zp = Zp(Nx+1:end, :);
Zh = eye(nPhi) - ZAB\ZAB;
Zh = Zh(Nx+1:end, :);

Zblk = []; X = [];
for i=1:Nx
    Zblk = blkdiag(Zblk, Zh);
    X    = [X x0(i)*eye(nPhi)];
end

zeroIdx = find(~PsiSupp(Nx+1:end,:));
F       = Zblk(zeroIdx, :);
g       = -vec(Zp);
g       = g(zeroIdx);
yp      = ZAB\IO*x0;
IFF     = eye(size(F,2)) - F\F;

cvx_begin quiet
expression traj2(nPhi)
variable w(nPhi*Nx)

yn    = [zeros(Nx,1); Zh*X*(F\g) + Zh*X*IFF*w];
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

%% Nullspace formulation (second version)
suppIdx = find(PsiSupp);
numSupp = length(suppIdx);

Zblk = []; X = [];
for i=1:Nx
    Zblk = blkdiag(Zblk, ZAB);
    X    = [X x0(i)*eye(nPhi)];
end

H   = Zblk(:, suppIdx);
h   = vec(IO);
Xm  = X(:, suppIdx);
IHH = eye(numSupp) - H\H;
yp  = Xm*(H\h);

cvx_begin quiet
expression traj3(nPhi)
variable v(numSupp)

yn    = Xm*IHH*v;
traj3 = yp + yn;

obj3 = 0;
for k = 1:T
    kx         = get_range(k, Nx); % rows of Psi representing Psix
    vect       = QSqrt*traj3(kx);
    obj3 = obj3 + vect'*vect;
    end
for k = 1:T-1
    ku         = Nx*T + get_range(k, Nu); % rows of Psi representing Psiu
    vect       = RSqrt*traj3(ku);
    obj3 = obj3 + vect'*vect;
end

minimize(obj3)
cvx_end

%% These should be very small
trajDiff1 = norm(traj1 - traj2) / norm(traj1)
objDiff1  = norm(obj1 - obj2) / norm(obj1)

trajDiff2 = norm(traj1 - traj3) / norm(traj1)
objDiff2  = norm(obj1 - obj3) / norm(obj1)
