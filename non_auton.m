function [Tx] = non_auton(t, x)
global N dim lam T;

T(3*dim*N+1:4*dim*N, 2*dim*N+1:3*dim*N) = lam * exp(t);

Tx = T * x;
end
