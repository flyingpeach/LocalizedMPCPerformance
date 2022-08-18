function mtx = get_local_subspace(sys, x0, params, adjustLocality)
% We use tFIR and locality size from params
% Space of local trajectories is proportional to the size of Image(mtx)
% adjustLocality: adjust definition of locality to work with grid example

Nx   = sys.Nx; Nu = sys.Nu; T = params.tFIR_;
nPhi = Nx*T + Nu*(T-1);
ZAB  = sparse(get_constraint_zab(sys, T));

if adjustLocality
    PsiSupp = get_sparsity_psi(sys, params, adjustLocality);
else
    PsiSupp = get_sparsity_psi(sys, params);
end
PsiSupp = PsiSupp(:, 1:Nx);
suppIdx = find(PsiSupp);

% TODO: there is a faster way to construct these
H = []; 
for i=1:Nx
    H = blkdiag(H, ZAB);
end

X = [];
for i=1:Nx
    X = [X x0(i)*speye(nPhi)];
end
X = X(Nx+1:end, suppIdx);

H = H(:, suppIdx);


IO = [eye(Nx); zeros(Nx*(T-1), Nx)];
h = vec(IO);
% Check if there is even a solution!
EPS = 1e-8;
test = H\h;
if norm(H*test - h, 'fro') > EPS % no solution exists
    mtx = 0;
    return; 
end


mtx1 = speye(length(suppIdx)) - H\H;
mtx  = X * mtx1;

% Zh = speye(nPhi) - ZAB\ZAB; 
% Zh = sparse(Zh(Nx+1:end, :));
% 
% EPS = 1e-8;
% Zh(abs(Zh) < EPS) = 0; 
% 
% nz      = nPhi - Nx;
% 
% IFFs = cell(Nx, 1); % Carries blocks of I - F\F
% fprintf('Calculating I - Fp*F\n');
% for i=1:Nx
%     startIdx = (i-1)*nz+1;
%     endIdx   = i*nz;    
%     zeroHere = zeroIdx(zeroIdx >= startIdx & zeroIdx <= endIdx) - (i-1)*nz;
%     IFFs{i}   = speye(nPhi) - Zh(zeroHere,:)\Zh(zeroHere,:);
%     IFFs{i}(abs(IFFs{i}) < EPS) = 0;
% end
% 
% fprintf('Calculating X*(I-Fp*F)\n');
% XFF = sparse(nPhi, Nx*nPhi); % X*(I-F\F)
% for i=1:Nx
%     idxEnd   = nPhi*i;
%     idxStart = idxEnd - nPhi + 1;
%     XFF(:, idxStart:idxEnd) = x0(i)*IFFs{i};
% end
% 
% fprintf('Calculating entire matrix\n');
% mtx     = Zh*XFF;
% 
% % Implementing sparse zero-column finder because matlab doesn't do it
% zeroCols = false(size(mtx,2), 1);
% for i=1:size(mtx,2)
%     if isempty(find(mtx(:,i), 1))
%         zeroCols(i) = true;
%     end
% end
% mtx(:,zeroCols) = [];
% 
% 
