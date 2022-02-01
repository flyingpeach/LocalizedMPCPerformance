function clMaps = state_fdbk_sls_custom_sparsity(sys, params, RSupp, MSupp, varargin)
% System level synthesis with state feedback
% Returns 
%    clMaps   : SLSOutputs containing system responses and other info
% Inputs
%    sys      : LTISystem containing system matrices
%    params   : SLSParams containing parameters
%    varargin expects clMapsIn (for two-step SLS only)
clMapsIn = [];
try 
    clMapsIn = varargin{1};
end
if isempty(clMapsIn) % Otherwise we are not solving for closed-loop maps
    fprintf('Finding closed loop maps\n\n');
end
params.print()
params.sanity_check()

T = params.T_;

cvx_precision low
cvx_begin

expression Rs(sys.Nx, sys.Nx, T)
expression Ms(sys.Nu, sys.Nx, T)

if params.approx_
    expression Deltas(sys.Nx, sys.Nx, T)
end

% populate decision variables for ease-of-use
R = cell(T, 1); 
M = cell(T, 1);
for t = 1:T
    R{t} = Rs(:,:,t); M{t} = Ms(:,:,t); 
end
if params.approx_
   Delta = cell(T, 1);
   for t = 1:T
       Delta{t} = Deltas(:,:,t);
   end
end

% enforce sparsity constraints
suppSizeR = sum(sum(RSupp));
suppSizeM = sum(sum(MSupp));
variable RMSupps(suppSizeR*T + suppSizeM*T)
spot = 0;
for t = 1:T
    R{t}(RSupp) = RMSupps(spot+1:spot+suppSizeR);
    spot = spot + suppSizeR;

    M{t}(MSupp) = RMSupps(spot+1:spot+suppSizeM);
    spot = spot + suppSizeM;
end

objective = get_total_objective(sys, params, R, M, clMapsIn);

% achievability  / approx achievability constraints
R{1} == eye(sys.Nx);

if params.approx_
    for t=1:T-1
        Delta{t} = R{t+1} - sys.A*R{t} - sys.B2*M{t};
    end
    Delta{T} = - sys.A*R{T} - sys.B2*M{T};
    % regularization for stability
    objective = objective + params.stabCoeff_ * get_stab_norm(Delta);
else
    R{T} == zeros(sys.Nx, sys.Nx);
    for t=1:T-1
        R{t+1} == sys.A*R{t} + sys.B2*M{t};
    end
end

minimize(objective);
cvx_end

% outputs
clMaps              = CLMaps();
clMaps.acts_        = get_acts_rfd(sys, M); % rfd actuator selection
clMaps.R_           = R;
clMaps.M_           = M;
clMaps.solveStatus_ = cvx_status;

if strcmp(cvx_status, 'Solved')
    fprintf(['Solved!', '\n\n']);
else
    sls_warning(['Solver exited with status ', cvx_status]);
end
end


% local functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function acts = get_acts_rfd(sys, M)
tol = 1e-4;

acts = [];
    for i=1:sys.Nu
        if norm(vec(M{1}(i,:)),2) > tol
            acts = [acts; i];
        end
    end
end