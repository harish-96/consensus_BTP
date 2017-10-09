%Parameters

global N dim freqs A;

N = 5;
dim = 3;

%Create an undirected, acyclic and connected graph
A = zeros(N);
s1 = 2:N; s2 = [1];
for i=2:N
    id = randi([1 length(s2)]);
    s2 = [s2, i];
    A(i, s2(id)) = 1;
    A(s2(id), i) = 1;
end

G = graph(A);
plot(G);

% Compute the matrices associated with the graph
freqs = 5*rand(N,N)*2*pi;
D = diag(sum(A));
L = D - A;
lam = kron(L, eye(dim));

% Initialize random measurement bias
beta = rand(dim*N, 1) - 0.5;

% Initialize the positions, velocities and control inputs randomly
pos0 = 10 * rand(dim*N, 1) - 5;
vel0 = 1 * rand(dim*N, 1) - 0.5;
beta_ad0 = rand(dim*N, 1) - 0.5;
% y0 = lam*pos0;

% Error in initial parameter estimates: Tilde quantities
beta_t0 = beta - beta_ad0;

% Solving the equations using ode45
% auton = @(t, x) T*x;
init_cond = [pos0', vel0', beta_t0']';

tmax = 1000; nTime = 10000; dt = tmax/nTime; 
tvals = linspace(0, tmax, nTime);
[tvals, sol] = ode45(@(t1, x)non_auton_one_var(t1,x), tvals, init_cond);

% Computing the positions from the velocity history
pos = zeros(nTime, dim*N);
pos(1, :) = pos0;
for t_step=2:nTime
    pos(t_step, :) = pos(t_step-1, :) + dt * sol(t_step-1, dim*N+1:2*dim*N);
end

% Plotting the results
figure;
plot(tvals, sol(:, 1:dim:dim*N)+sol(:, 2*dim*N+1:dim:3*dim*N));
ylabel('y+$\widetilde\beta$', 'Interpreter', 'latex')
figure;
plot(tvals, sol(:, dim*N+1:dim:2*dim*N))
ylabel('velocity')
for i=1:dim
    figure;
    plot(tvals, pos(:, i:dim:end))
    str = sprintf('x%d position',i);
    title(str)
end

% Computing the bias
% y(inf) = -beta_t(inf). So, beta = beta_ad(inf) - y(inf)

y_inf = sol(end, 1:dim*N)';
beta_ad_inf=(trapz(-sol(:, dim*N+1:2*dim*N), 1)*dt)' + beta_ad0;
beta_estimate = beta_ad_inf - y_inf;
display(['The error in the bias estimate is ', num2str(norm(beta_estimate - beta))])