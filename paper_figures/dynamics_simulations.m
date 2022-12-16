clear all; clc;
warning off;

%% User-specified parameters
numSims        = 20;
tFIR           = 15;
actDensity     = 1.0;
gridSize       = 11;
connectThresh  = 0.6;
Ts             = 0.2;
tHorizon       = 20;

% Need to pick seeds such that grid is fully connected
% These work for gridSize = 11, connectThresh=0.6
seeds = [700, 702, 707, 710, 712, ...
         715, 717, 718, 719, 720, ...
	     727, 730, 731, 733, 734, ...
	     735, 738, 740, 741, 742];

plotTopology = false;

% Generate and visualize plants
numNodes    = gridSize*gridSize;
systems     = cell(numSims, 1);

for i=1:numSims
    seed = seeds(i);
    rng(seed);

    [adjMtx, nodeCoords, susceptMtx, inertiasInv, dampings] = generate_grid_topology(gridSize, connectThresh, seed);
    
    if plotTopology
        figure(i);
        plot_graph(adjMtx, nodeCoords, 'k');
    end

    numActs       = round(actDensity*numNodes);
    actuatedNodes = randsample(numNodes, numActs);    
    systems{i}    = generate_grid_plant(actuatedNodes, adjMtx, susceptMtx, inertiasInv, dampings, Ts);
        
    % Use custom communication structure for grid
    systems{i}.AComm = adjust_grid_sys_locality(systems{i}.A);
end

%% Check locality sizes
% params       = MPCParams();
% params.tFIR_ = tFIR;
% locSizes     = zeros(1, numSims);
% 
% for i=1:numSims
%     fprintf('Checking system %d of %d\n', i, numSims);
%     locSizes(i) = get_ideal_locality(systems{i}, params);
% end
% 
% % Remember that our "locality" is actually different from paper
% % paper: d=0 means only self communication
% % us   : d=1 means only self communication
% 
% locSizes - 1

%% HACK: remove this, uncomment previous section
% we have run this already with locSizes
locSizes = 2*ones(1, numSims);

%% Dynamics
pList  = cell(numSims,1);
xsGlob = cell(numSims,1);
usGlob = cell(numSims,1);
xsLoc  = cell(numSims,1);
usLoc  = cell(numSims,1);

objGlobs  = nan(numSims,1);
objLocs   = nan(numSims,1);
costDiffs = nan(numSims,1);

for i=1:numSims
    fprintf('Running sim %d of %d \n', i, numSims);
    sys = systems{i};
    params = MPCParams();
    params.tFIR_  = tFIR;    
    
    params.QSqrt_ = diag(rand(sys.Nx, 1)); % [0,1]
    params.RSqrt_ = diag(rand(sys.Nu, 1)); % [0,1]
    
    params.stateConsMtx_ = eye(sys.Nx);
    params.stateUB_      = 20 * ones(sys.Nx, 1); % Freq constraints, loose
    params.stateUB_(1:2:sys.Nx) = 4; % Phase constraints
    params.stateLB_ = -params.stateUB_;

    params.locality_ = locSizes(i);    
    pList{i}         = params; % Save this info
    
    x0 = 4*(rand(sys.Nx,1) - 0.5); % [-2, 2]
    
    % Global MPC
    fprintf('Calculating global MPC\n')
    [xsGlob{i}, usGlob{i}] = global_mpc(sys, x0, params, tHorizon);

    % Local MPC
    fprintf('Calculating local MPC\n')
    [xsLoc{i}, usLoc{i}] = local_mpc(sys, x0, params, tHorizon);
    
    % Can do this in postprocessing, calculating here because
    % we want to see live progress/debug if needed
    objGlobs(i)  = get_cost_fn(params, xsGlob{i}, usGlob{i});
    objLocs(i)   = get_cost_fn(params, xsLoc{i}, usLoc{i});
    costDiffs(i) = (objLocs(i) - objGlobs(i))/objGlobs(i);
    fprintf('Cost diff: %.4e\n\n', costDiffs(i)); % Should be small    
end

save('data/dynamics_simulations.mat');

%% Plots

costDiffs
max(costDiffs)

% Compare differences of states and inputs
% Choose one to plot if needed?

