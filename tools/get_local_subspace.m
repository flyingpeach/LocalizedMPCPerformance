function mtx = get_local_subspace(sys, x0, params, adjustLocality, eps)
% We use tFIR and locality size from params
% Space of local trajectories is proportional to the size of Image(mtx)
% adjustLocality: adjust definition of locality to work with grid example
% eps           : how to determine whether solution exists

Nx   = sys.Nx; Nu = sys.Nu; T = params.tFIR_;
nPhi = Nx*T + Nu*(T-1);
ZAB  = sparse(get_constraint_zab(sys, T));
IO   = [eye(Nx); zeros(Nx*(T-1), Nx)];
h    = vec(IO);

if adjustLocality
    PsiSupp = get_sparsity_psi(sys, params, adjustLocality);
else
    PsiSupp = get_sparsity_psi(sys, params);
end
PsiSupp  = PsiSupp(:, 1:Nx);
suppIdx  = find(PsiSupp);
suppSize = length(suppIdx);

mtx = [];
for i=1:Nx
    idxStart = (i-1)*nPhi+1;
    idxEnd   = i*nPhi;
    myIdx    = suppIdx(suppIdx >= idxStart & suppIdx <= idxEnd) - (i-1)*nPhi;
    nCols    = length(myIdx);
    
    % Check if there is a solution
    Hi       = ZAB(:,myIdx);
    hi       = h((i-1)*Nx*T+1:i*Nx*T);
    testSol  = Hi\hi;
    
    % Note: using max instead of 2-norm to make this check less
    %       dimension-dependent
    if max(abs(Hi*testSol - hi)) > eps % No solution exists
        mtx = 0;
        return;
    end
    
    IHHi     = speye(nCols) - Hi\Hi;
    zeroCols = false(nCols,1); % Get rid of zero columns
    for j=1:nCols
        if isempty(find(IHHi(:,j),1))
            zeroCols(j) = true;
        end
    end
    IHHi(:,zeroCols) = [];
    
    xi  = x0(i)*speye(nPhi); 
    mtx = [mtx xi(Nx+1:end,myIdx)*IHHi];
end