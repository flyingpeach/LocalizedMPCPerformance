function locality = get_ideal_locality(sys, params)
% locality: size of local communication region (d-hop neighbors) for which
%           the size of the trajectory space is unchanged compared to 
%           global communication, assuming x0 is dense
%           convention: d=1 means self-communication only
% params : MPCParams(); the locality_ field will be populated with the 
%          ideal locality

% Values of x0 need to be nonzero; specific value doesn't matter
x0 = ones(sys.Nx, 1);

for locality=2:sys.Nx
    fprintf('Checking locality size %d\n', locality);
    params.locality_ = locality;
    rankRatio        = get_rank_ratio(sys, x0, params);
    if rankRatio == 1
        break;
    end
end