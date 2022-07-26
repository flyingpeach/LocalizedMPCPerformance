function mtx = get_local_subspace(sys, x0, params)
% We use tFIR and locality size from params
% Space of local trajectories is proportional to the size of Image(mtx)
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

Ed            = speye((nPhi-Nx)*Nx);
nonZero       = find(PsiSupp);
Ed(nonZero,:) = [];
F             = sparse(Ed*Z);

mtx1 = F\F; % this will give rank deficiency warnings; it's ok
mtx  = Z2*X*(eye(nPhi*Nx)-mtx1);