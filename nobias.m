%Parameters

global N dim sigma gamma lam T;

N = 10;
dim = 3;
sigma = 10;
gamma = 1;

%Create an undirected, acyclic and connected graph
A = zeros(N);
s1 = 2:N; s2 = [1];
for i=2:N
    id = randi([1 length(s2)]);
    s2 = [s2, i];
    A(i, s2(id)) = 1;
    A(s2(id), i) = 1;
end
% A = [0 1 0 0 0; 1 0 1 0 0; 0 1 0 1 1; 0 0 1 0 0; 0 0 1 0 0];
G = graph(A);
plot(G);graph(A)

% Compute the matrices associated with the graph
D = diag(sum(A));
L = D - A;
lam = kron(L, eye(dim));
% Initialize random measurement bias
beta = rand(dim*N, 1);
k = kron(D, eye(dim))*beta;

% Initialize the positions, velocities and control inputs randomly
pos0 = 10 * rand(dim*N, 1);
vel0 = 1*rand(dim*N, 1);

% Setting up the matrix differential equation. dX/dt = T*X

T = kron([0 1; 0 -1], eye(dim*N));

T(dim*N+1:2*dim*N, 1:dim*N) = -lam;

% Solving the equations using ode45
auton = @(t, x) T*x;
init_cond = [pos0;vel0];

tmax = 10; nTime = 100; dt = tmax/nTime; 
t = linspace(0, tmax, nTime);
[t, sol] = ode45(@(t, x)auton(t,x), t, init_cond);

% Computing the positions from the velocity history
pos = zeros(nTime, dim*N);
pos(1, :) = pos0;
for t_step=2:nTime
    pos(t_step, :) = pos(t_step-1, :) + dt * sol(t_step-1, dim*N+1:end);
end

% Plotting the results
% figure;
% plot(t, sol(:, 1:dim:dim*N) - sol(:, dim*N+1:dim:2*dim*N));
% figure;
for i=1:dim
    figure;
    plot(t, pos(:, i:dim:end))
    str = sprintf('x%d position',i);
    title(str)
end
