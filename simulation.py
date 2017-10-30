import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import networkx as nx

dim = 3
N = 5
alpha = 10 * np.eye(N)
beta_order = 10
pos_order = 1
vel_order = 1
exp_order = 1

def generateGraph(N):
    A = np.zeros((N, N));
    range(1, N); s2 = [0];
    for i in range(1, N): 
        idx = np.random.randint(0, len(s2))
        s2.append(i)
        A[i, s2[idx]] = 1 
    A[s2[idx], i] = 1
    G2 = nx.from_numpy_matrix(A)
    nx.draw_networkx(G2)

    D = np.diag(np.sum(A, axis=0))
    L = D - A
    return A, L, D

A, L, D = generateGraph(N)

T = np.zeros((3*N, 3*N))
T[:N, N:2*N] = np.eye(N)
T[N:2*N, N:2*N] = -(L + alpha)#-np.eye(N)
T[N:2*N, :N] = -(L + alpha)
T[N:2*N, 2*N:] = -(D + alpha)
T[2*N:, N:2*N] = (D + alpha)

def diffeq(x, t):
    T1 = T.copy()
    T1[2*N:, N:2*N] *= np.sin(t)
    T[N:2*N, 2*N:] *= np.sin(t)
    T1 = np.kron(T1, np.eye(dim))
    return np.dot(T1, x)

pos0 = np.random.uniform(-0.5, 0.5, dim*N) * pos_order
vel0 = np.random.uniform(-0.5, 0.5, dim*N) * vel_order
beta = np.random.uniform(-0.5, 0.5, dim*N) * beta_order
init_cond = np.concatenate([pos0, vel0, beta])

time = np.linspace(0, 100, 1000)
sol = odeint(diffeq, init_cond, time)

x = sol[:, :N*dim:dim]
y = sol[:, 1:N*dim:dim]
z = sol[:, 2:N*dim:dim]
vx = sol[:, N*dim:2*N*dim:dim]

plt.figure()
plt.plot(time, x, linewidth=0.7)
plt.ylabel('x')
plt.figure()
plt.plot(time, y, linewidth=0.7)
plt.ylabel('y')
plt.figure()
plt.plot(time, z, linewidth=0.7)
plt.ylabel('z')
plt.figure()
plt.plot(time, vx, linewidth=0.7)
plt.ylabel('x-velocity')
plt.show()
