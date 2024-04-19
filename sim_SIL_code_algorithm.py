import numpy as np
import array
x_m_prev = np.array([[0],[0]])
A_m_prev = np.array([[0, 0],[0, 0]])
B_m_prev = np.array([[0],[0]])

def adaptive_estimation(current, speed, u):
    global x_m_prev, A_m_prev, B_m_prev
    dt = 0.01
    P = np.array([[1, 1],[1, 3]])
    Gamma_a = np.array([[100, 0],[0, 100]])
    Gamma_b = np.array([[100, 0],[0, 100]])
    x = np.array([[current],[speed]])
    x_m = x_m_prev + (np.dot(A_m_prev, x_m_prev) + np.dot(B_m_prev, u)) * dt
    e = x - x_m
    A_m = A_m_prev + np.dot(np.dot(np.dot(Gamma_a, P), e), x_m.T)*dt
    B_m = B_m_prev + np.dot(np.dot(Gamma_b, P), e)*u*dt
    
    x_m_prev = x_m
    A_m_prev = A_m
    B_m_prev = B_m
    out = np.concatenate((A_m, B_m), axis = 1)
    return out