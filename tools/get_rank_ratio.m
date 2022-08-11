function rankRatio = get_rank_ratio(sys, x0, params, adjustLocality)
% adjustLocality: adjust definition of locality to work with grid example

EPS = 1e-8; % hardcoded
mtx = get_local_subspace(sys, x0, params, adjustLocality);

fprintf('Calculating rank\n');
n     = (sys.Nu*(params.tFIR_-1)); % Max possible rank
svals = svds(mtx, n); % Largest singular values

zeroSv = find(svals < EPS, 1); % = rank+1
if isempty(zeroSv)
    rankRatio = 1; % largest n singular values are all nonzero
else
    rankRatio = (zeroSv-1)/n;
end

% For comparison
% rankRatioFull = rank(full(mtx), EPS) / (sys.Nu*(params.tFIR_-1))