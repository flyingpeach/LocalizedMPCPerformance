function rankRatio = get_rank_ratio(sys, x0, params, adjustLocality)
% adjustLocality: adjust definition of locality to work with grid example

mtx = get_local_subspace(sys, x0, params, adjustLocality);

fprintf('Calculating rank\n');
rankRatio = rank(full(mtx)) / (sys.Nu*(params.tFIR_-1))