%% Use SLS to find closed-loop map (no comm constraints)
% Set up closed-loop constraints (for Phix)
PhixSupp  = false(Nx, Nx);
PhiuSupp  = true(Nu, Nx);
clampsExt = [0 clamps Nx+1]; % "clamps" at end of chain
lastClamp = clampsExt(1);
nextClamp = clampsExt(2);

for i=1:Nx
    if ismember(i, clamps) 
        lastClampIdx = find(clampsExt==i);
        nextClamp    = clampsExt(lastClampIdx+1);        
        
        % A disturbance at the clamp can spread toward the next 2 clamps
        PhixSupp(lastClamp+1:nextClamp-1, i) = true;

        lastClamp    = i;
    else
        PhixSupp(lastClamp+1:nextClamp-1, i) = true;
    end
end

% Set up parameters for SLS functions
sys.B2  = B;
sys.Nu = Nu;

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

%% Simulation 
% Cost
cost_s = 0;
for k=1:Ts
    cost_s = cost_s + norm([statePenSqrt*Rs{k} inputPenSqrt*Ms{k}'], 'fro').^2;
end
