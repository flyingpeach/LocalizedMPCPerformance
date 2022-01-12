%% Custom SLS design
% Same idea as clamp: do not let disturbances through

% First, determine locality structure
RSupp = false(Nx, Nx);

clampsExt    = [0 clamps Nx+1]; % "clamps" at end of chain
lastClamp    = clampsExt(1);
nextClamp    = clampsExt(2);

for i=1:Nx
    if ismember(i, clamps) 
        lastClampIdx = find(clampsExt==i);
        nextClamp    = clampsExt(lastClampIdx+1);        
        
        % A disturbance at the clamp can spread toward the next 2 clamps
        RSupp(lastClamp+1:nextClamp-1, i) = true;

        lastClamp    = i;
    else
        RSupp(lastClamp+1:nextClamp-1, i) = true;
    end
end

MSupp     = true(Nu, Nx);
suppSizeR = sum(sum(RSupp));
suppSizeM = sum(sum(MSupp));

% The following is adapted from state_fdbk_sls.m from the toolbox
cvx_begin quiet

expression R(Nx, Nx, T) 
expression M(Nu, Nx, T)

% populate decision variables for ease-of-use
Rs = cell(T, 1); 
Ms = cell(T, 1);
for k = 1:T
    Rs{k} = R(:,:,k); Ms{k} = M(:,:,k); 
end

variable RMSupps(suppSizeR*T + suppSizeM*T)
spot = 0;
for t = 1:T
    Rs{t}(RSupp) = RMSupps(spot+1:spot+suppSizeR);
    spot = spot + suppSizeR;

    Ms{t}(MSupp) = RMSupps(spot+1:spot+suppSizeM);
    spot = spot + suppSizeM;
end
    
objective = 0;
for k=1:T
    vect = vec([statePenSqrt*Rs{k}; inputPenSqrt*Ms{k}]);
    objective = objective + vect'*vect;
end

% achievability constraints
Rs{1} == eye(Nx);
Rs{T} == zeros(Nx, Nx);
for t=1:T-1
    Rs{t+1} == A*Rs{t} + B*Ms{t};
end

minimize(objective);
cvx_end

%% Simulation
Ts = length(Rs);
 
% Cost
cost_s = 0;
for k=1:Ts
    cost_s = cost_s + norm([statePenSqrt*Rs{k} inputPenSqrt*Ms{k}'], 'fro').^2;
end