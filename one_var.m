%Parameters
global N dim lam T D;
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
D = diag(sum(A));
L = D - A;
lam = kron(L, eye(dim));

% Initialize random measurement bias
beta = rand(dim*N, 1) - 0.5;

% Initialize the positions, velocities and control inputs randomly
pos0 = 10 * rand(dim*N, 1) - 5;
vel0 = 1 * rand(dim*N, 1) - 0.5;
beta_ad0 = rand(dim*N, 1) - 0.5;
y0 = lam*pos0;

% Error in initial parameter estimates: Tilde quantities
beta_t0 = beta - beta_ad0;

% Setting up the matrix differential equation. dX/dt = T*X
%(xi, vel, beta_t)
T = kron([0 1 0;
          -1 -1 -1;
          0 1 0], eye(dim*N));
T(dim*N+1:2*dim*N, 1:dim*N) = -lam;
T(dim*N+1:2*dim*N, 2*dim*N+1:3*dim*N) = -kron(D, eye(dim));
T(2*dim*N+1:3*dim*N, dim*N+1:2*dim*N) = kron(D, eye(dim));

% Solving the equations using ode45
auton = @(t, x) T*x;
init_cond = [y0', vel0', beta_t0']';

tmax = 50; nTime = 5000; dt = tmax/nTime; 
t = linspace(0, tmax, nTime);
[t, sol] = ode45(@(t, x)auton(t,x), t, init_cond);
% [t, sol] = ode45(@(t, x)non_auton_one_var(t,x), t, init_cond);

% Computing the positions from the velocity history
pos = zeros(nTime, dim*N);
pos(1, :) = pos0;
for t_step=2:nTime
    pos(t_step, :) = pos(t_step-1, :) + dt * sol(t_step-1, dim*N+1:2*dim*N);
end

% Plotting the results
figure;
plot(t, L*sol(:, 1:dim:dim*N)'+D*sol(:, 2*dim*N+1:dim:3*dim*N)');
ylabel('$\Lambda\xi+\Delta\widetilde\beta$', 'Interpreter', 'latex')
figure;
plot(t, sol(:, dim*N+1:dim:2*dim*N))
ylabel('velocity')
for i=1:dim
    figure;
    plot(t, pos(:, i:dim:end))
    str = sprintf('x%d position',i);
    title(str)
end

% Computing the bias
% y(inf) = -beta_t(inf). So, \Delta beta = \Delta beta_ad(inf) - y(inf)

y_inf = lam * sol(end, 1:dim*N)';
beta_ad_inf = trapz(kron(D, eye(dim))*-sol(:, dim*N+1:2*dim*N)', 2)*dt + beta_ad0;
beta_estimate = beta_ad_inf - kron(diag(1./diag(D)), eye(dim))*y_inf;

display(['The error in the bias estimate is ', num2str(norm(beta_estimate - beta))])