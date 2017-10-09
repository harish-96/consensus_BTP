%Parameters

global N dim sigma gamma lam T;

N = 6;
dim = 3;
sigma = 1;
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
beta = rand(dim*N, 1) - 0.5;
k = rand(dim*N, 1) - 0.5;

% Initialize the positions, velocities and control inputs randomly
pos0 = 10 * rand(dim*N, 1) - 5;
% pos0(end-dim+1:end) = 1;
vel0 = 1 * rand(dim*N, 1) - 0.5;
beta_ad0 = rand(dim*N, 1) - 0.5;
k_ad0 = rand(dim*N, 1) - 0.5;
y0 = lam*pos0 + kron(D, eye(dim))*beta_ad0;

% Error in initial parameter estimates: Tilde quantities
k_t0 = k - k_ad0;
beta_t0 = beta - beta_ad0;

% Setting up the matrix differential equation. dX/dt = T*X
%(beta, k, velocity, position)
T = kron([0 0 0 sigma;
          sigma*gamma -sigma*gamma -gamma sigma;
          -1 1 -1 -1;
          0 0 0 -sigma], eye(dim*N));

T(3*dim*N+1:4*dim*N, 2*dim*N+1:3*dim*N) = lam;
% T(end-dim+1:end, end-4*dim+1:end) = kron([0 0 0 1], eye(dim));
% T(end-N*dim-dim+1:end-N*dim, :) = 0;

% Solving the equations using ode45
auton = @(t, x) T*x;
init_cond = [beta_t0;k_t0;vel0;y0];

tmax = 5; nTime = 100; dt = tmax/nTime; 
t = linspace(0, tmax, nTime);
[t, sol] = ode45(@(t, x)auton(t,x), t, init_cond);

% Computing the positions from the velocity history
pos = zeros(nTime, dim*N);
pos(1, :) = pos0;
for t_step=2:nTime
    pos(t_step, :) = pos(t_step-1, :) + dt * sol(t_step-1, 2*dim*N+1:3*dim*N);
end

% Plotting the results
figure;
plot(t, sol(:, 1:dim:dim*N) - sol(:, dim*N+1:dim:2*dim*N));
figure;
for i=1:dim
%     figure;
    plot(t, pos(:, i:dim:end))
    str = sprintf('x%d position',i);
    title(str)
end
