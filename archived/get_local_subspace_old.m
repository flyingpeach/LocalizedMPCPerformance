function mtx = get_local_subspace_old(sys, x0, params, adjustLocality)
% We use tFIR and locality size from params
% Space of local trajectories is proportional to the size of Image(mtx)
% adjustLocality: adjust definition of locality to work with grid example

% NOTE: This function doesn't check if the locality constraints can
%       be satisfied; it only checks for the theoretical size of the
%       nullspace, if it exists.

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

IFFs = cell(Nx, 1); % Carries blocks of I - F\F
fprintf('Calculating I - Fp*F\n');
for i=1:Nx
    startIdx = (i-1)*nz+1;
    endIdx   = i*nz;    
    zeroHere = zeroIdx(zeroIdx >= startIdx & zeroIdx <= endIdx) - (i-1)*nz;
    IFFs{i}   = speye(nPhi) - Zh(zeroHere,:)\Zh(zeroHere,:);
    IFFs{i}(abs(IFFs{i}) < EPS) = 0;
end

fprintf('Calculating X*(I-Fp*F)\n');
XFF = sparse(nPhi, Nx*nPhi); % X*(I-F\F)
for i=1:Nx
    idxEnd   = nPhi*i;
    idxStart = idxEnd - nPhi + 1;
    XFF(:, idxStart:idxEnd) = x0(i)*IFFs{i};
end

fprintf('Calculating entire matrix\n');
mtx     = Zh*XFF;

% Implementing sparse zero-column finder because matlab doesn't do it
zeroCols = false(size(mtx,2), 1);
for i=1:size(mtx,2)
    if isempty(find(mtx(:,i), 1))
        zeroCols(i) = true;
    end
end
mtx(:,zeroCols) = [];


