%% LQR design
Q = statePenSqrt^2*eye(Nx);
R = inputPenSqrt^2*eye(Nu);

S = idare(A, B, Q, R);
cost_l = 0;
for i = 1:Nx
    ei     = zeros(Nx,1);
    ei(i)  = 1;
    cost_l = cost_l + ei'*S*ei;
end

Kl = -(B'*S*B + R)\(B'*S*A);
