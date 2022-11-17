function [rankRatio, parTime, rankTime] = get_rank_ratio(sys, x0, params, eps)
% adjustLocality: adjust definition of locality to work with grid example

[mtx, parTime] = get_local_subspace(sys, x0, params, eps);

tic;
rankRatio = rank(full(mtx)) / (sys.Nu*(params.tFIR_-1));
rankTime  = toc;

fprintf('Rank ratio: %.2f\n', rankRatio);

