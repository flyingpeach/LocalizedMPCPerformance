%% Setup
% Remember that our "locality" is actually different from paper
% paper: d=0 means only self communication
% code : d=1 means only self communication
clear; clc;

color{1} = [0 204 0]/255; % Green
color{2} = [255 51 204]/255; % Pink

fontsizeLabel  = 14;
fontsizeLegend = 12;
fontsizeAxis   = 12;
fontsizeTitle  = 14;

LineWidth  = 4;
MarkerSize = LineWidth*8;

%% Figure: runtime

load('data/scan_time_vs_network_size.mat');
parMeans  = mean(parTimes, 2)';
parStds   = std(parTimes, 0, 2)';
rankMeans = mean(rankTimes, 2)';
rankStds  = std(rankTimes, 0, 2)';

figure(1); box on;
% HACK to overcome MATLAB bug
semilogy(linspace(0,100,5), linspace(1,1.1,5), 'Color', [1 1 1], 'HandleVisibility', 'off');

% Mean runtimes
hold on;
semilogy(gridSizes.^2, rankMeans, '.-', 'LineWidth', LineWidth, ... 
         'MarkerSize', MarkerSize, 'Color', color{2}, ...
         'MarkerEdge', color{2}, 'MarkerFaceColor', color{2});
semilogy(gridSizes.^2, parMeans, '.-', 'LineWidth', LineWidth, ... 
         'MarkerSize', MarkerSize, 'Color', color{1}, ...
         'MarkerEdge', color{1}, 'MarkerFaceColor', color{1});

% Standard deviations of runtimes
parFill = fill([gridSizes.^2,fliplr(gridSizes.^2)], ...
     [parMeans - parStds, fliplr(parMeans + parStds)], color{1});
alpha(.4);
set(parFill, 'EdgeAlpha', 0);
rankFill = fill([gridSizes.^2,fliplr(gridSizes.^2)], ...
     [rankMeans - rankStds, fliplr(rankMeans + rankStds)], color{2});
alpha(.4);
set(rankFill, 'EdgeAlpha', 0);

% Legend
leg = legend('Rank determination (non-parallelized)', 'Matrix construction (parallelized)');
set(leg, 'Interpreter', 'latex'); 
set(leg, 'Fontsize', fontsizeLegend);

% Figure settings
set(gca, 'FontSize',fontsizeAxis)
title({'\textbf{Algorithm runtime}'}, 'Interpreter', 'latex', 'Fontsize', fontsizeTitle);
xlabel({'\# of subsystems in the network'}, 'Interpreter', 'latex', 'Fontsize', fontsizeLabel);
ylabel({'Avg. total runtime (s)'}, 'Interpreter', 'latex', 'Fontsize', fontsizeLabel);
xlim([5 125]);
ylim([3e-4 1e1]);