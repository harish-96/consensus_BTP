function [Tx] = non_auton_one_var(t, x)
global N dim freqs A;
t
A_local = A * 0.1*sin(freqs*t);
D = diag(sum(A_local));
L = D - A_local;
lam = kron(L, eye(dim));
T = kron([0 1 0;
          -1 -1 -1;
          0 1 0], eye(dim*N));

T(dim*N+1:2*dim*N, 1:dim*N) = -lam;
T(dim*N+1:2*dim*N, 2*dim*N+1:3*dim*N) = -kron(D, eye(dim));
T(2*dim*N+1:3*dim*N, dim*N+1:2*dim*N) = kron(D, eye(dim));

Tx = T*x;
end
