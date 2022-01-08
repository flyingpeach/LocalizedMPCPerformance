%% Clamp + controller design
Kc = zeros(Nu, Nx); % K for clamped design
for i=1:numClamps
    Kc(i,:) = -A(clamps(i), :);
end

% LQR penalties
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