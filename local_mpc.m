function [xs, us] = local_mpc(sys, x0, params, tHorizon)
% params   : MPCParams() object containing horizon, objective, constraints
%            note: terminal cost/constr not accounted for in this case
% tHorizon : how long to run MPC for
% Adapted from mpc_centralized; specialized for nominal only

params.sanity_check_cent();

% For convenience
Nx = sys.Nx; Nu = sys.Nu; T = params.tFIR_;
nPhi  = Nx*T + Nu*(T-1);
QSqrt = params.QSqrt_;
RSqrt = params.RSqrt_;

xs = zeros(Nx, tHorizon);
us = zeros(Nu, tHorizon-1);

xs(:,1) = x0;

PsiSupp = get_sparsity_psi(sys, params); 
PsiSupp = PsiSupp(:, 1:Nx);
suppSizePsi = sum(sum(PsiSupp));
ZAB         = get_constraint_zab(sys, T);
IO          = [eye(Nx); zeros(Nx*(T-1), Nx)];
[H, h]      = get_constraint_h(sys, params);

for t=1:tHorizon-1
    fprintf('Calculating time %d of %d\n', t+1, tHorizon);
    obj = 0; % objective function
        
    cvx_begin quiet
    expression Psi(nPhi, Nx)
    variable PsiSuppVals(suppSizePsi, 1)
    Psi(PsiSupp) = PsiSuppVals;     

    for k=1:T-1
        ku   = Nx*T + get_range(k, Nu); % rows of Psi representing Psiu
        vect = RSqrt*Psi(ku, 1:Nx)*x0;
        obj  = obj + vect'*vect;
    end
    for k=1:T
        kx   = get_range(k, Nx); % rows of Psi representing Psix
        vect = QSqrt*Psi(kx, 1:Nx)*x0;
        obj  = obj + vect'*vect;
    end
    
    % Dynamics constraints
    % Note: have to use this formulation b/c otherwise  cvx reports infeasible
    %       for ZAB*Psi == IO, even if we set low precision
    EPS = 1e-8;
    norm(ZAB*Psi - IO, 'fro') <= EPS;
    
    % State/input constraints
    H*Psi(:, 1:Nx)*x0 <= h;
    
    minimize(obj)
    cvx_end
    
    uNext = Psi(Nx*T+1:Nx*T+Nu,:)*xs(:,t);
    
    if isnan(uNext)
        error('MPC infeasible');
    end
    
    us(:,t)   = uNext(:,1); % apply control action
    xs(:,t+1) = sys.A*xs(:,t) + sys.B2*us(:,t); % propagate dynamics
end