%% Design
S = idare(A, B, statePen*eye(Nx), inputPen*eye(Nu));
cost_l = 0;
for i = 1:Nx
    ei     = zeros(Nx,1);
    ei(i)  = 1;
    cost_l = cost_l + ei'*S*ei;
end

Kl = -(B'*S*B + inputPen*eye(Nu))\(B'*S*A);
