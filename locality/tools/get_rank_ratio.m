function rankRatio = get_rank_ratio(sys, x0, params)
EPS = 1e-8; % hardcoded

[C1, C2, C3] = get_local_subspace(sys, x0, params);
rankRatio    = rank(C1*C3, EPS) / (sys.Nu*(params.tFIR_-1));