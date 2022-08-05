function mtx = get_local_subspace(sys, x0, params, adjustLocality)
% We use tFIR and locality size from params
% Space of local trajectories is proportional to the size of Image(mtx)
% adjustLocality: adjust definition of locality to work with grid example

Nx   = sys.Nx; Nu = sys.Nu; T = params.tFIR_;
nPhi = Nx*T + Nu*(T-1);

ZAB  = get_constraint_zab(sys, T);

% Unneeded calculations for rank comparison; needed for subspace vector
% IO = [eye(Nx); zeros(Nx*(T-1), Nx)];
% Zp = ZAB\IO;              
% Zp = Zp(Nx+1:end, :);

Zh = speye(nPhi) - ZAB\ZAB; 
Zh = sparse(Zh(Nx+1:end, :));

if adjustLocality
    PsiSupp = get_sparsity_psi(sys, params, adjustLocality);
else
    PsiSupp = get_sparsity_psi(sys, params);
end
PsiSupp = PsiSupp(Nx+1:end, 1:Nx);

nz   = nPhi - Nx;
Zblk = sparse(nz*Nx, nPhi*Nx);
for i=1:Nx
    Zblk((i-1)*nz+1:i*nz,(i-1)*nPhi+1:i*nPhi) = Zh;
end

% This takes too much memory
% Zblk = kron(speye(Nx), Zh); % Block diagonal matrix of Nx blocks 

ZblkDensity = nnz(Zblk)/size(Zblk,1)/size(Zblk,2)

X = sparse(nPhi, Nx*nPhi);
for i=1:Nx
    idxEnd   = nPhi*i;
    idxStart = idxEnd - nPhi + 1;    
    X(:, idxStart:idxEnd) =  x0(i)*speye(nPhi);
end

zeroIdx = find(~PsiSupp);
F       = Zblk(zeroIdx,:);
mtx1    = F\F; % this will give rank deficiency warnings; it's ok
mtx     = Zh*X*(speye(nPhi*Nx)-mtx1);