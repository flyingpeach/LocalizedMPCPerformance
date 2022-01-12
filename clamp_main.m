clear; clc; close;

%% User-determined parameters
patchSize  = 3; % Odd number
numPatches = 3;
alpha      = 0.2; % how much neighbors affect each other
rho        = 1.2; % stability

% SLS-only params
T = 80;

% LQR weights
statePenSqrt = 1;
inputPenSqrt = 1;

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
ctrller_lqr; % get Kl, cost_l
ctrller_naive_clamp; % get Kc, cost_c, Rc, Mc
ctrller_sls; % get cost_s, Rs, Ms

cost_clamp_normalized = cost_c/cost_l
cost_sls_normalized   = cost_s/cost_l

%% Plot disturbance response
distNode = 3;

x0 = zeros(Nx, 1);  x0(distNode) = 1;
T  = min(Ts, Tc);
xs_c = zeros(Nx, T); us_c = zeros(Nx, T);
xs_l = zeros(Nx, T); us_l = zeros(Nx, T);
xs_s = zeros(Nx, T); us_s = zeros(Nx, T);

for k=1:T
     xs_l(:,k) = (A+B*Kl)^(k-1)*x0;     us_l(:,k) = B*Kl*xs_l(:,k);     

     xs_c(:,k) = Rc{k}(:,distNode); us_c(:,k) = B*Mc{k}(:,distNode);
     xs_s(:,k) = Rs{k}(:,distNode); us_s(:,k) = B*Ms{k}(:,distNode);

end

simcost_l = norm([statePenSqrt*xs_l inputPenSqrt*us_l], 'fro').^2;
simcost_c = norm([statePenSqrt*xs_c inputPenSqrt*us_c], 'fro').^2;
simcost_s = norm([statePenSqrt*xs_s inputPenSqrt*us_s], 'fro').^2;
simcost_c_norm = simcost_c / simcost_l
simcost_s_norm = simcost_s / simcost_l

plot_heat_map(xs_l, us_l, 'Optimal');
plot_heat_map(xs_c, us_c, 'Clamped');
plot_heat_map(xs_s, us_s, 'SLS');

