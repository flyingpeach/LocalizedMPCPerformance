%% Setup
% Remember that our "locality" is actually different from paper
% paper: d=0 means only self communication
% code : d=1 means only self communication
clear; clc;

colors{1} = [0 204 0]/255; % Green
colors{2} = [0 153 255]/255; % Blue
colors{3} = [255 51 204]/255; % Pink

fontsizeLabel  = 14;
fontsizeLegend = 10;
fontsizeAxis   = 12;
fontsizeTitle  = 14;

LineWidth  = 4;
MarkerSize = LineWidth*8;

figure(2);

%% Actuation density
load('data/scan_d_vs_act_dens.mat');
actMeans = mean(locSizes,2)' - 1;
actStds  = std(locSizes, 0, 2)';
color    = [0 0 0];
subplot(1,3,1); box on; hold on;
plot(actDensities, actMeans, '.-', 'LineWidth', LineWidth, ... 
    'MarkerSize', MarkerSize, 'Color', color, ...
    'MarkerEdge', color, 'MarkerFaceColor', color);

color    = [0 0 0];
actFill = fill([actDensities, fliplr(actDensities)], ...
       [actMeans - actStds, fliplr(actMeans + actStds)], color);
alpha(.4);
set(actFill, 'EdgeAlpha', 0);

set(gca, 'FontSize',fontsizeAxis)
xlabel({'Actuation density'}, 'Interpreter', 'latex', 'Fontsize', fontsizeLabel);
ylabel({'Minimum locality size'}, 'Interpreter', 'latex', 'Fontsize', fontsizeLabel);
xlim([0.1 1.1]);
ylim([0.8 6.2]);


%% Network sizes
load('data/scan_d_vs_network_size.mat');
netMeans = {};
netStds  = {};

for actIdx=1:numActDens
    netMeans{actIdx} = mean(locSizes{actIdx},2)' - 1;
    netStds{actIdx}  = std(locSizes{actIdx}, 0, 2)';
end

subplot(1,3,2); box on; hold on;

for actIdx=numActDens:-1:1
    color = colors{actIdx};
    plot(gridSizes.^2, netMeans{actIdx}, '.-', 'LineWidth', LineWidth, ... 
        'MarkerSize', MarkerSize, 'Color', color, ...
        'MarkerEdge', color, 'MarkerFaceColor', color);
end

for actIdx=numActDens:-1:1
    color = colors{actIdx};
    horFill{actIdx} = fill([gridSizes.^2, fliplr(gridSizes.^2)], ...
           [netMeans{actIdx} - netStds{actIdx}, fliplr(netMeans{actIdx} + netStds{actIdx})], color);
    alpha(.4);
    set(horFill{actIdx}, 'EdgeAlpha', 0);
end

set(gca,'FontSize',fontsizeAxis)
xlabel({'\# of subsystems', 'in the network'}, 'Interpreter', 'latex', 'Fontsize', fontsizeLabel);
xlim([10 126]);
ylim([0.8 6.2]);

%% Horizon size
load('data/scan_d_vs_horizon_size.mat');
horMeans = {};
horStds  = {};

for actIdx=1:numActDens
    horMeans{actIdx} = mean(locSizes{actIdx},2)' - 1;
    horStds{actIdx}  = std(locSizes{actIdx}, 0, 2)';
end

subplot(1,3,3); box on; hold on;

for actIdx=numActDens:-1:1
    color = colors{actIdx};
    plot(Ts, horMeans{actIdx}, '.-', 'LineWidth', LineWidth, ... 
        'MarkerSize', MarkerSize, 'Color', color, ...
        'MarkerEdge', color, 'MarkerFaceColor', color);
end

for actIdx=numActDens:-1:1
    color = colors{actIdx};
    horFill{actIdx} = fill([Ts, fliplr(Ts)], ...
           [horMeans{actIdx} - horStds{actIdx}, fliplr(horMeans{actIdx} + horStds{actIdx})], color);
    alpha(.4);
    set(horFill{actIdx}, 'EdgeAlpha', 0, 'HandleVisibility', 'off');
end

set(gca,'FontSize',fontsizeAxis)
xlabel({'Horizon size'}, 'Interpreter', 'latex', 'Fontsize', fontsizeLabel);
xlim([3 32]);
ylim([0.8 6.2]);

leg = legend('60\% Actuated', '80\% Actuated', '100\% Actuated');
set(leg, 'Fontsize', fontsizeLegend);
set(leg, 'Interpreter', 'latex'); 