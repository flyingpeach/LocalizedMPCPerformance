function rankRatio = get_rank_ratio(sys, x0, params)
EPS = 1e-8; % hardcoded

mtx       = get_local_subspace(sys, x0, params);
rankRatio = rank(mtx, EPS) / (sys.Nu*(params.tFIR_-1));