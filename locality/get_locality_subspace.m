function [C1, C2, C3] = get_locality_subspace(sys, x0, params)
% We use tFIR and locality size from params
% Space of nonlocal trajectories: Image(C1)
% Space of local trajectories: Images(C1*C3)
% Space of reachable Y = C1*C2 + C1*C3*[free variable]

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
F = Ed*Z; Fp = pinv(F); 

C1 = Z2*X; C2 = Fp*G; C3 = eye(nPhi*Nx)-Fp*F;