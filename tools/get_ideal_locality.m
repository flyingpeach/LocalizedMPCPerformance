function locality = get_ideal_locality(sys, params, adjustLocality, varargin)
% locality: size of local communication region (d-hop neighbors) for which
%           the size of the trajectory space is unchanged compared to 
%           global communication, assuming x0 is dense
%           convention: d=1 means self-communication only
% params : MPCParams(); the locality_ field will be populated with the 
%          ideal locality
% adjustLocality: adjust definition of locality to work with grid example
% eps           : how to determine whether solution exists fpr loc subspace


if length(varargin) > 0
    eps = varargin{1}; 
else
    eps = 1e-8; % Default value
end

% Values of x0 need to be nonzero; specific value doesn't matter
x0 = ones(sys.Nx, 1);

maxLoc = sys.Nx;
if adjustLocality
    maxLoc = sys.Nx/2; % Grid; two states per node
end

for locality=2:maxLoc
    fprintf('Checking locality size %d\n', locality);
    params.locality_ = locality;
    rankRatio        = get_rank_ratio(sys, x0, params, adjustLocality, eps);
    if rankRatio == 1
        break;
    end
end