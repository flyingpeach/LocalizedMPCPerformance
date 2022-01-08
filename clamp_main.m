clear; clc; close;

%% User-determined parameters
patchSize  = 3; % Odd number
numPatches = 3;
alpha      = 0.2; % how much neighbors affect each other
rho        = 1.2; % stability

% LQR weights
statePen   = 1;
inputPen   = 1;

%% Plant setup for chain
numClamps = numPatches-1; % edges of chain don't need clamps
numCtrls  = numPatches;
Nx        = numPatches*patchSize + numClamps;
ctrlRad   = floor(patchSize/2);

clamps = patchSize+1:patchSize+1:Nx;
ctrls  = ceil(patchSize/2):patchSize+1:Nx;

junk  = 1;
sys = LTISystem(); sys.Nx = Nx;
generate_dbl_stoch_chain(sys, rho, junk, alpha); 
A = full(sys.A);

inputs = [clamps ctrls];
Nu     = length(inputs);
B = zeros(Nx, Nu);
for i=1:Nu
    B(inputs(i),i) = 1;
end

%% Synthesize different controllers to compare
ctrller_naive_clamp; % get Kc, cost_c, Rc, Mc
ctrller_lqr; % get Kl, cost_l

cost_clamp_normalized = cost_c/cost_l

%% Plot disturbance response
distNode = 5;

x0   = zeros(Nx, 1);  x0(distNode) = 1;
xs_c = zeros(Nx, Tc); us_c = zeros(Nx, Tc);
xs_l = zeros(Nx, Tc); us_l = zeros(Nx, Tc);

for k=1:Tc
     xs_c(:,k) = Rc{k}(:,distNode); us_c(:,k) = B*Mc{k}(:,distNode);
     xs_l(:,k) = (A+B*Kl)^(k-1)*x0;     us_l(:,k) = B*Kl*xs_l(:,k);     
end

simcost_l = norm([statePen*xs_l inputPen*us_l], 'fro').^2;
simcost_c = norm([statePen*xs_c inputPen*us_c], 'fro').^2;
simcost_c_norm = simcost_c / simcost_l

plot_heat_map(xs_c, us_c, 'Clamped');
plot_heat_map(xs_l, us_l, 'Optimal');

