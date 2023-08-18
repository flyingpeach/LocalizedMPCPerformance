%% Setup
% Remember that our "locality" is actually different from paper
% paper: d=0 means only self communication
% code : d=1 means only self communication
clear; clc;

colors{1} = [0 204 0]/255; % Green
colors{2} = [0 153 255]/255; % Blue
colors{3} = [255 51 204]/255; % Pink

fontsizeLabel  = 14;
fontsizeLegend = 12;
fontsizeAxis   = 12;
fontsizeTitle  = 14;

LineWidth  = 4;
MarkerSize1 = LineWidth*2;
MarkerSize2 = LineWidth*3.5;

Marker1 = 'o-';
Marker2 = 'x-';

figure(1);

%% Runtime vs. network size
load('data/runtime_vs_network_size.mat');
parMeans  = mean(parTimes, 2)';
parStds   = std(parTimes, 0, 2)';
rankMeans = mean(rankTimes, 2)';
rankStds  = std(rankTimes, 0, 2)';

subplot(1,2,1); box on;
% HACK to overcome MATLAB bug
semilogy(linspace(0,100,5), linspace(1,1.1,5), 'Color', [1 1 1], 'HandleVisibility', 'off');

% Mean runtimes
hold on; % Have to do this after first semilogy otherwise linear plot
semilogy(gridSizes.^2, rankMeans, Marker1, 'LineWidth', LineWidth, ... 
         'MarkerSize', MarkerSize1, 'Color', colors{3}, ...
         'MarkerEdge', colors{3}, 'MarkerFaceColor', colors{3});
semilogy(gridSizes.^2, parMeans, Marker2, 'LineWidth', LineWidth, ... 
         'MarkerSize', MarkerSize2, 'Color', colors{1}, ...
         'MarkerEdge', colors{1}, 'MarkerFaceColor', colors{1});

% Standard deviations of runtimes
parFill = fill([gridSizes.^2,fliplr(gridSizes.^2)], ...
     [parMeans - parStds, fliplr(parMeans + parStds)], colors{1});
alpha(.4);
set(parFill, 'EdgeAlpha', 0);
rankFill = fill([gridSizes.^2,fliplr(gridSizes.^2)], ...
     [rankMeans - rankStds, fliplr(rankMeans + rankStds)], colors{3});
alpha(.4);
set(rankFill, 'EdgeAlpha', 0);

% Legend
leg = legend('Rank determination', 'Matrix construction');
set(leg, 'Interpreter', 'latex'); 
set(leg, 'Fontsize', fontsizeLegend);

% Figure settings
set(gca, 'FontSize',fontsizeAxis)
xlabel({'\# of subsystems in the network'}, 'Interpreter', 'latex', 'Fontsize', fontsizeLabel);
ylabel({'Avg. total runtime (s)'}, 'Interpreter', 'latex', 'Fontsize', fontsizeLabel);
xlim([14 123]);
ylim([3e-4 15]);

%% Runtime vs. horizon size
load('data/runtime_vs_horizon_size.mat');
parMeans  = mean(parTimes, 2)';
parStds   = std(parTimes, 0, 2)';
rankMeans = mean(rankTimes, 2)';
rankStds  = std(rankTimes, 0, 2)';

subplot(1,2,2); box on;
% HACK to overcome MATLAB bug
semilogy(linspace(0,100,5), linspace(1,1.1,5), 'Color', [1 1 1], 'HandleVisibility', 'off');

% Mean runtimes
hold on; % Have to do this after first semilogy otherwise linear plot
semilogy(Ts, rankMeans, Marker1, 'LineWidth', LineWidth, ... 
         'MarkerSize', MarkerSize1, 'Color', colors{3}, ...
         'MarkerEdge', colors{3}, 'MarkerFaceColor', colors{3});
semilogy(Ts, parMeans, Marker2, 'LineWidth', LineWidth, ... 
         'MarkerSize', MarkerSize2, 'Color', colors{1}, ...
         'MarkerEdge', colors{1}, 'MarkerFaceColor', colors{1});

% Standard deviations of runtimes
parFill = fill([Ts,fliplr(Ts)], ...
     [parMeans - parStds, fliplr(parMeans + parStds)], colors{1});
alpha(.4);
set(parFill, 'EdgeAlpha', 0);
rankFill = fill([Ts,fliplr(Ts)], ...
     [rankMeans - rankStds, fliplr(rankMeans + rankStds)], colors{3});
alpha(.4);
set(rankFill, 'EdgeAlpha', 0);

% Figure settings
set(gca, 'FontSize',fontsizeAxis)
xlabel({'Horizon length'}, 'Interpreter', 'latex', 'Fontsize', fontsizeLabel);
xlim([4 31]);
ylim([3e-4 15]);