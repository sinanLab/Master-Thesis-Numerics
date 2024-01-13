# Define function to perform optimization
import numpy as np
from scipy.optimize import minimize
from scipy.integrate import odeint
def optimize_model(real_data, t_fit, n_iter):

    # Define model function
    def model(y, t, kappa_h, kappa_r, rho_1, rho_2, rho_3, chi_h, chi_r, tau, omega_1, omega_2, omega_3, xi_h, gamma_1, gamma_2, gamma_3, psi_h, psi_r):
        S_h, E_h, A_h, I_h, Q_h, R_h, S_r, E_r, I_r = y
        N_h = 1425624477 # 2:04/ 5-7-2023 / https://www.worldometers.info/world-population/China-population/
        
        # Define initial conditions
        E_h0 = 5    # initial number of infected
        A_h0 = 5
        I_h0 = 1
        Q_h0 = 1    # initial number of recovered
        R_h0 = 0
        S_h0 = N_h - (E_h0+A_h0+I_h0+Q_h0+R_h0)

        E_r0 = 100
        I_r0 = 50
        S_r0 = S_h0*10e-2   # this condition is given from this paper: https://doi.org/10.3390/math11051121
        N_r = S_r0+E_r0+I_r0
        
        dS_hdt = kappa_h*N_h -rho_1*I_r*S_h/N_h - rho_2*I_h*S_h/N_h - chi_h*S_h + tau*Q_h
        dE_hdt = rho_1*I_r*S_h/N_h + rho_2*I_h*S_h/N_h - (omega_1 + omega_2 + omega_3 + chi_h)*E_h
        dA_hdt = omega_1*E_h - (xi_h+chi_h+gamma_1)*A_h
        dI_hdt = omega_2*E_h + xi_h*A_h - (chi_h + psi_h + gamma_2)*I_h
        dQ_hdt = omega_3*E_h - (tau+gamma_3+psi_h+chi_h)*Q_h
        dR_hdt = gamma_1*A_h +gamma_2*I_h+gamma_3*Q_h-chi_h*R_h
        dS_rdt = kappa_r*N_r - rho_3*S_r*I_r/N_r - chi_r*S_r
        dE_rdt = rho_3*S_r*I_r/N_r - (chi_r + omega_3)*E_r
        dI_rdt = omega_3*E_r - (chi_r + psi_r)*I_r
        return [dS_hdt, dE_hdt, dA_hdt, dI_hdt, dQ_hdt, dR_hdt, dS_rdt, dE_rdt, dI_rdt]
        

    N_h = 1425624477 # 2:04/ 5-7-2023 / https://www.worldometers.info/world-population/China-population/
            
    # Define initial conditions
    E_h0 = 5    # initial number of infected
    A_h0 = 5
    I_h0 = 1
    Q_h0 = 1    # initial number of recovered
    R_h0 = 0
    S_h0 = N_h - (E_h0+A_h0+I_h0+Q_h0+R_h0)

    E_r0 = 100
    I_r0 = 50
    S_r0 = S_h0*10e-2   # this condition is given from this paper: https://doi.org/10.3390/math11051121
    N_r = S_r0+E_r0+I_r0
    y0 = S_h0, E_h0, A_h0, I_h0, Q_h0, R_h0, S_r0, E_r0, I_r0

    kappa_h_guesses = np.random.rand(n_iter)
    kappa_r_guesses = np.random.rand(n_iter)
    rho_1_guesses = np.random.rand(n_iter)
    rho_2_guesses = np.random.rand(n_iter)
    rho_3_guesses = np.random.rand(n_iter)
    chi_h_guesses = np.random.rand(n_iter)
    chi_r_guesses = np.random.rand(n_iter)
    tau_guesses = np.random.rand(n_iter)
    omega_1_guesses = np.random.rand(n_iter)
    omega_2_guesses = np.random.rand(n_iter)
    omega_3_guesses = np.random.rand(n_iter)
    xi_h_guesses = np.random.rand(n_iter)
    gamma_1_guesses = np.random.rand(n_iter)
    gamma_2_guesses = np.random.rand(n_iter)
    gamma_3_guesses = np.random.rand(n_iter)
    psi_h_guesses = np.random.rand(n_iter)
    psi_r_guesses = np.random.rand(n_iter)

    # Define function to fit data
    def fit_model(kappa_h, kappa_r, rho_1, rho_2, rho_3, chi_h, chi_r, tau, omega_1, omega_2, omega_3, xi_h, gamma_1, gamma_2, gamma_3, psi_h, psi_r):
        return odeint(model, y0, t_sim, args=(kappa_h, kappa_r, rho_1, rho_2, rho_3, chi_h, chi_r, tau, omega_1, omega_2, omega_3, xi_h, gamma_1, gamma_2, gamma_3, psi_h, psi_r))[:,3][t_fit]
    
    # Define function to compute least squares error
    def least_squares_error(params):
        kappa_h, kappa_r, rho_1, rho_2, rho_3, chi_h, chi_r, tau, omega_1, omega_2, omega_3, xi_h, gamma_1, gamma_2, gamma_3, psi_h, psi_r = params
        prediction = fit_model(kappa_h, kappa_r, rho_1, rho_2, rho_3, chi_h, chi_r, tau, omega_1, omega_2, omega_3, xi_h, gamma_1, gamma_2, gamma_3, psi_h, psi_r)
        er = np.sum((real_data - prediction)**2)
        return er
    
    # Use optimization algorithm to find best-fit parameters from all initial guesses
    best_error = float('inf')
    rightbound = 4
    bounds = [(0.001, rightbound),(0.001, rightbound),(0.001, rightbound),(0.001, rightbound),(0.001, rightbound),(0.001, rightbound),(0.001, rightbound),(0.001, rightbound),(0.001, rightbound),(0.001, rightbound),(0.001, rightbound),(0.001, rightbound),(0.001, rightbound),(0.001, rightbound),(0.001, rightbound),(0.001, rightbound),(0.001, rightbound)]
    for i in range(n_iter):
        init_guess = [kappa_h_guesses[i], kappa_r_guesses[i], rho_1_guesses[i], rho_2_guesses[i], rho_3_guesses[i], chi_h_guesses[i], chi_r_guesses[i], tau_guesses[i], omega_1_guesses[i], omega_2_guesses[i], omega_3_guesses[i], xi_h_guesses[i], gamma_1_guesses[i], gamma_2_guesses[i], gamma_3_guesses[i], psi_h_guesses[i], psi_r_guesses[i]]
        res = minimize(least_squares_error, init_guess, bounds=bounds)
        if res.fun < best_error:
            best_error = res.fun
            best_params = res.x
    
    # Simulate SIR model with best-fit parameters
    solution = odeint(model, y0, t_sim, args=(best_params[0], best_params[1], best_params[2], best_params[3], best_params[4], best_params[5], best_params[6], best_params[7], best_params[8], best_params[9], best_params[10], best_params[11], best_params[12], best_params[13], best_params[14], best_params[15], best_params[16]))
    S_h, E_h, A_h, I_h, Q_h, R_h, S_r, E_r, I_r = solution.T
    
    return best_params, I_h