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

EPS = 1e-8;
Zh(abs(Zh) < EPS) = 0; 

if adjustLocality
    PsiSupp = get_sparsity_psi(sys, params, adjustLocality);
else
    PsiSupp = get_sparsity_psi(sys, params);
end
PsiSupp = PsiSupp(Nx+1:end, 1:Nx);
zeroIdx = find(~PsiSupp);
nz      = nPhi - Nx;

FFs = cell(Nx, 1); % Carries blocks of F\F
for i=1:Nx
    startIdx = (i-1)*nz+1;
    endIdx   = i*nz;    
    zeroHere = zeroIdx(zeroIdx >= startIdx & zeroIdx <= endIdx) - (i-1)*nz;
    FFs{i}   = Zh(zeroHere,:)\Zh(zeroHere,:);
end

XFF = sparse(nPhi, Nx*nPhi); % X*(I-F\F)
for i=1:Nx
    idxEnd   = nPhi*i;
    idxStart = idxEnd - nPhi + 1;
    XFF(:, idxStart:idxEnd) = x0(i)*(speye(nPhi) - FFs{i});
end

mtx     = Zh*XFF;