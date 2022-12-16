function [xs, us] = global_mpc(sys, x0, params, tHorizon)
% params   : MPCParams() object containing horizon, objective, constraints
%            note: terminal cost/constr not accounted for in this case
% tHorizon : how long to run MPC for

params.sanity_check_cent();

% For convenience
Nx = sys.Nx; Nu = sys.Nu; T = params.tFIR_;
QSqrt = params.QSqrt_;
RSqrt = params.RSqrt_;

xs = zeros(Nx, tHorizon);
us = zeros(Nu, tHorizon-1);

xs(:,1) = x0;

for t=1:tHorizon-1
    obj = 0; % objective function
    
    cvx_begin quiet
    variable xs_t(Nx, T) % forecasted trajectory
    variable us_t(Nu, T-1)
    
    xs_t(:,1) == xs(:,t); % current state

    for k=1:T-1
        obj = obj + us_t(:,k)'*RSqrt*RSqrt*us_t(:,k);
        
        % Dynamics constraints
        xs_t(:,k+1) == sys.A*xs_t(:,k) + sys.B2*us_t(:,k);
        
        % Input constraints
        if params.has_input_cons()
            params.inputConsMtx_ * us_t(:,k) <= params.inputUB_;
            params.inputConsMtx_ * us_t(:,k) >= params.inputLB_;
        end
    end
    
    for k=1:T
        obj = obj + xs_t(:,k)'*QSqrt*QSqrt*xs_t(:,k);
        
        % State constraints
        if params.has_state_cons()
            params.stateConsMtx_ * xs_t(:,k) <= params.stateUB_;
            params.stateConsMtx_ * xs_t(:,k) >= params.stateLB_;           
        end
    end        

    minimize(obj)
    cvx_end
    if isnan(us_t(:,1))
        error('MPC infeasible');
    end
    
    us(:,t)   = us_t(:,1); % apply control action
    xs(:,t+1) = sys.A*xs(:,t) + sys.B2*us(:,t); % propagate dynamics
end