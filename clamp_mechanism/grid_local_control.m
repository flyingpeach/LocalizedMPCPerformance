clear; clc;

% User-defined params
specRad      = 1.5;
seed         = 421;
statePenSqrt = 5;
inputPenSqrt = 1;
Ts = 40;

rng(seed);

gridSize = 4;
Nx       = gridSize * gridSize;
nodeCoords = zeros(Nx, 2);
for i=1:gridSize
    theseRows = (i-1)*gridSize+1:i*gridSize;
    nodeCoords(theseRows, 1) = i*ones(gridSize, 1);
    nodeCoords(theseRows, 2) = 1:gridSize;
end 

adjMtx = eye(Nx);

% Hand-generated for test purposes
links = [1 2; 2 3; 3 4; 6 7; 7 8; 9 10; 10 11; 11 12; 13 14; 14 15;
         1 5; 2 6; 3 7; 4 8; 6 10; 9 13; 10 14; 11 15; 12 16];
for i=1:length(links)
    adjMtx(links(i,1), links(i,2)) = 1;
    adjMtx(links(i,2), links(i,1)) = 1;
end

% Visualize graph to double-check
% plot_graph(adjMtx, nodeCoords, 'k'); axis equal;

suppASize = sum(sum(adjMtx));
suppA     = rand(suppASize, 1);
A         = zeros(Nx, Nx);
A(adjMtx ~=0) = suppA;

% Scale 
A = A ./ max(abs(eig(A))) * specRad;

%% Design naive controller
numPatches = 2;
patches    = cell(numPatches, 1); % exclude clamped nodes
patches{1} = [1 2 3 4 5 7 8];
patches{2} = [9 10 11 12 13 14 15 16];

clamps = [6];  

ctrllers    = cell(numPatches, 1);
ctrllers{1} = [2 8];
ctrllers{2} = [12 14];

numClamps = length(clamps);
numCtrls  = 0;
inputs    = clamps;
for i=1:numPatches
    numCtrls = numCtrls + length(ctrllers{i});
    inputs = [inputs ctrllers{i}];
end

Nu = numClamps + numCtrls;
B  = zeros(Nx, Nu);
for i=1:Nu
    B(inputs(i),i) = 1;
end

Kc        = zeros(Nu, Nx);
for i=1:length(clamps)
    Kc(i,:) = -A(clamps(i),:);
end

% LQR patches
for i=1:numPatches
    patchIdx  = patches{i};
    patchSize = length(patches{i});
    ALoc      = A(patchIdx, patchIdx);

    numLocCtrls = length(ctrllers{i});
    BLoc      = zeros(patchSize, numLocCtrls);
    for j=1:numLocCtrls
        actIdxLoc = find(patches{i} == ctrllers{i}(j));
        BLoc(actIdxLoc, j) = 1;
    end
    
    QLoc = statePenSqrt^2*eye(patchSize);
    RLoc = inputPenSqrt^2*eye(numLocCtrls);
    S    = idare(ALoc, BLoc, QLoc, RLoc);
    Kc(numClamps+(2*i-1):numClamps+2*i, patchIdx) = -(BLoc'*S*BLoc + RLoc)\(BLoc'*S*ALoc);
end

%% Simulate naive controller
[Rc, Mc] = get_cl_maps_from_k(A, B, Kc);
Tc       = length(Rc);
  
% Cost
cost_c = 0;
for k=1:Tc
    cost_c = cost_c + norm([statePenSqrt*Rc{k} inputPenSqrt*Mc{k}'], 'fro').^2;
end

%% Design LQR controller
Q = statePenSqrt^2*eye(Nx);
R = inputPenSqrt^2*eye(Nu);

S = idare(A, B, Q, R);
cost_l = 0;
for i = 1:Nx
    ei     = zeros(Nx,1);
    ei(i)  = 1;
    cost_l = cost_l + ei'*S*ei;
end

Kl = -(B'*S*B + R)\(B'*S*A);

%% Use SLS to find closed-loop map (no comm constraints)
% Set up closed-loop constraints (for Phix)
PhixSupp  = false(Nx, Nx);
PhiuSupp  = false(Nu, Nx);

PhixSupp(clamps, :) = true; % Clamps can be affected by both patches
PhixSupp(patches{1}, [patches{1} 6]) = true;
PhixSupp(patches{2}, [patches{2} 6]) = true;

% Order: 6, 2, 8, 12, 14
PhiuSupp(1,:) = true;
PhiuSupp(2:3, [patches{1} 6]) = true;
PhiuSupp(4:5, [patches{2} 6]) = true;

% Set up parameters for SLS functions
sys = LTISystem();
sys.A = A; sys.B2  = B;
sys.Nu = Nu;

sys.Nx = Nx;
sys.Nw  = sys.Nx;
sys.Nz  = sys.Nu + sys.Nx;
sys.B1  = eye(sys.Nx); 
sys.C1  = [statePenSqrt*speye(sys.Nx); sparse(sys.Nu, sys.Nx)];
sys.D11 = sparse(sys.Nz, sys.Nw);
sys.D12 = [inputPenSqrt*sparse(sys.Nx, sys.Nu); speye(sys.Nu)];
sys.sanity_check();

slsParams       = SLSParams();
slsParams.T_    = Ts;
slsParams.add_objective(SLSObjective.H2, 1);

clMaps = state_fdbk_sls_custom_sparsity(sys, slsParams, PhixSupp, PhiuSupp);

Rs = clMaps.R_; Ms = clMaps.M_;

% Cost
cost_s = 0;
for k=1:Ts
    cost_s = cost_s + norm([statePenSqrt*Rs{k} inputPenSqrt*Ms{k}'], 'fro').^2;
end

%% Output costs
cost_clamp_normalized = cost_c/cost_l
cost_sls_normalized   = cost_s/cost_l

%% Plot disturbance response
distNode = 11;

x0 = zeros(Nx, 1);  x0(distNode) = 1;
T  = min(Ts, Tc);
xs_c = zeros(Nx, T); us_c = zeros(Nx, T);
xs_l = zeros(Nx, T); us_l = zeros(Nx, T);
xs_s = zeros(Nx, T); us_s = zeros(Nx, T);

for k=1:T
     xs_c(:,k) = Rc{k}(:,distNode); us_c(:,k) = B*Mc{k}(:,distNode);
     xs_l(:,k) = (A+B*Kl)^(k-1)*x0; us_l(:,k) = B*Kl*xs_l(:,k);  
     xs_s(:,k) = Rs{k}(:,distNode); us_s(:,k) = B*Ms{k}(:,distNode);
end

plot_heat_map(xs_c, us_c, 'Clamped');
plot_heat_map(xs_l, us_l, 'Optimal');
plot_heat_map(xs_s, us_s, 'SLS');
