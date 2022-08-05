function rankRatio = get_rank_ratio(sys, x0, params, adjustLocality)
% adjustLocality: adjust definition of locality to work with grid example

EPS = 1e-8; % hardcoded

mtx       = get_local_subspace(sys, x0, params, adjustLocality);
rankRatio = rank(full(mtx), EPS) / (sys.Nu*(params.tFIR_-1));