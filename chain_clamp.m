clear; clc;

%% Plant setup
Nx = 9;
alpha = 0.2;
rho   = 1.2; % stability
junk  = 1;

sys = LTISystem; sys.Nx = Nx;
generate_dbl_stoch_chain(sys, rho, junk, alpha); 
A = full(sys.A);

clamps    = [2, 8];
numClamps = length(clamps); 
ctrls     = [5];
numCtrls  = length(ctrls);
ctrlRad   = 2; 

inputs = [clamps ctrls];
Nu     = length(inputs);

B = zeros(Nx, Nu);
for i=1:Nu
    B(inputs(i),i) = 1;
end

%% Clamp + controller design
Kc = zeros(Nu, Nx);

for i=1:numClamps
    Kc(i,:) = -A(clamps(i), :);
end

% LQR penalties
statePen = 1;
inputPen = 1;

for i=1:numCtrls
    patchIdx = ctrls(i)-ctrlRad:ctrls(i)+ctrlRad;
    ALoc     = A(patchIdx, patchIdx);
    BLoc     = zeros(2*ctrlRad + 1, 1);
    BLoc(ctrlRad+1) = 1;
    
    QLoc = statePen*eye(size(ALoc));
    RLoc = inputPen;
    S    = idare(ALoc, BLoc, QLoc, RLoc);
    Kc(numClamps+i, patchIdx) = -(BLoc'*S*BLoc + RLoc)\(BLoc'*S*ALoc);
end

%% Simulation
[Rc, Mc] = get_cl_maps_from_k(A, B, Kc);
Tc       = length(Rc);
 
% Cost
cost_c = 0;
for k=1:Tc
    cost_c = cost_c + norm([sqrt(statePen)*Rc{k} sqrt(inputPen)*Mc{k}'], 'fro').^2;
end

% Compare to true LQR cost
S = idare(A, B, statePen*eye(Nx), inputPen*eye(Nu));
cost_k = 0;
for i = 1:Nx
    ei     = zeros(Nx,1);
    ei(i)  = 1;
    cost_k = cost_k + ei'*S*ei;
end

cost_k
cost_c_norm  = cost_c / cost_k

%% Plot disturbance response
distNode = 3;

xs_c = zeros(Nx, Tc);
us_c = zeros(Nx, Tc);
for k=1:Tc
     xs_c(:,k) = Rc{k}(:,distNode);
     us_c(:,k) = B*Mc{k}(:,distNode);
end

plot_heat_map(xs_c, us_c, '');