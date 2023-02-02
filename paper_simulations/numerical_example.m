clc; close all; clear;

A = [1 2 0;
     3 4 5;
     0 6 7];
B = [1 0;
     0 0;
     0 1];
Nx = 3;
Nu = 2;
T  = 1;

sys = LTISystem();
sys.A = A; sys.B2 = B;
sys.Nx = Nx; sys.Nu = Nu;

params       = MPCParams();
params.tFIR_ = T+1; % Code and paper use different conventions

locality = get_ideal_locality(sys, params)

%% Looking at individual matrices
eps = 1e-8; % zero-trimming

ZAB  = get_constraint_zab(sys, T);
IO   = [eye(Nx); zeros(Nx*(T-1), Nx)];
ZABp = pinv(ZAB);

a = ZABp*IO;
a(abs(a) < 1e-8) = 0;

b = eye(size(ZAB,2)) - ZABp*ZAB;
b(abs(b) < 1e-8) = 0;

PsiSupp  = get_sparsity_psi(sys, params);
PsiSupp  = PsiSupp(:, 1:Nx); % The rest of the matrix is for robust only

Zh = b(Nx+1:end,:);
Zhblk = blkdiag(Zh, Zh, Zh);
F = Zhblk([3, 5, 11, 14],:);

X = [eye(8), eye(8), eye(8)];

Zh*X*(eye(size(F,2)) - pinv(F)*F);

x0 = ones(Nx,1);

% Similar nonzeros whether we use pinv(H)*H or H\H
% Remember to edit get_local_subspace.m to see full untrimmed matrix
params.locality_ = 2;
G = get_local_subspace(sys, x0, params, eps);
full(G)

