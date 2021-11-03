clear; close all; clc; 

% user-specified
d = 3;

% specify system matrices
sys    = LTISystem;
sys.Nx = 10; sys.Nw = sys.Nx; 

% generate sys.A, sys.B2
alpha = 0.2; rho = 1; actDens = 1;
generate_dbl_stoch_chain(sys, rho, actDens, alpha); 

% objectives
QSqrt = eye(sys.Nx);
RSqrt = eye(sys.Nu);

sys.Nz  = sys.Nu + sys.Nx;
sys.B1  = eye(sys.Nx); 
sys.C1  = [sparse(QSqrt); sparse(sys.Nu, sys.Nx)];
sys.D11 = sparse(sys.Nz, sys.Nw);
sys.D12 = [sparse(sys.Nx, sys.Nu); sparse(RSqrt)];
sys.sanity_check();

% params for standard sls
slsParams    = SLSParams();
slsParams.T_ = 5;
slsParams.add_objective(SLSObjective.H2, 1); % H2 objective, coefficient = 1
slsParams.add_constraint(SLSConstraint.Locality, d); % locality = 3-hops

slsOut = state_fdbk_sls(sys, slsParams);
slsR = slsOut.R_;
slsM = slsOut.M_;

% params for h2 (uses MPC parameters)
h2params.locality_ = d;
h2params.QSqrt_ = QSqrt;
h2params.RSqrt_ = RSqrt;

[h2R, h2M] = h2_sls(sys, h2params);