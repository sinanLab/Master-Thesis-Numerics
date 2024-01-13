clear all; clc;
t0 = 0; tf = 10; h = 0.4; N = ceil(tf/h); t = t0:h:tf;
S_h = zeros(length(t), 1); E_h = zeros(length(t), 1); I_h = zeros(length(t), 1); H = zeros(length(t), 1);
R_h = zeros(length(t), 1); S_v = zeros(length(t), 1); E_v = zeros(length(t), 1); I_v = zeros(length(t), 1);

S_h(1) = 150; E_h(1) = 40; I_h(1) = 15; H(1) = 10; R_h(1) = 9; 
S_v(1) = 90; E_v(1) = 70; I_v(1) = 50;
%% parameters
% Pre-used
alpha_1 = 10.1; nu = 0.01; theta_1 = 0.009; theta_2 = 0.001; delta_1 = 0.004; delta_2 = 0.01;  eta_1 = 0.7; eta_2 = 0.03;
gamma_1 = 0.025;  alpha_2 = 8.23; 
% Assumed
gamma_2 = 0.63; rho = 0.10; Omega = 0.28; hbar = 0.55; 

%% Step size
wh = (1-exp(-h));

%% NSFD Scheme (Explicite Scheme)
for i = 1:N
    t(i+1) = t(i)+h;
    S_h(i+1) = (S_h(i)+wh*alpha_1)./(1+nu*theta_1*wh*I_v(i)+delta_1);
    E_h(i+1) = (E_h(i)+nu*theta_1*wh*S_h(i+1)*I_v(i))./(1+wh*(eta_1+rho+delta_1));
    I_h(i+1) = (I_h(i)+eta_1*wh*E_h(i+1))./(1+wh*(gamma_1+delta_1+Omega+hbar));
    H(i+1) = (H(i)+hbar*wh*I_h(i+1))./(1+wh*(delta_1+gamma_2));
    R_h(i+1) = (R_h(i)+wh*(gamma_1*I_h(i)+gamma_2*H(i)))./(1+delta_1*wh);
    S_v(i+1) = (S_v(i)+alpha_2*wh)./(1+wh*nu*theta_2*I_h(i)+delta_2);
    E_v(i+1) = (E_v(i)+nu*theta_2*wh*I_h(i)*S_v(i+1))./(1+wh*(delta_2+eta_2));
    I_v(i+1) = (I_v(i)+eta_2*wh*E_v(i+1))./(1+delta_2*wh);
end
%% Visualization
close all; clf;
figure(1)
plot(t, S_h, 'r-', 'LineWidth',2)
xlabel('Time [t]', 'Interpreter', 'latex', 'FontSize', 13)
ylabel('Susceptible Human Population', 'Interpreter', 'latex', 'FontSize', 13)

figure(2)
plot(t, E_h, 'r-', 'LineWidth',0.5)
xlabel('Time [t]', 'Interpreter', 'latex', 'FontSize', 13)
ylabel('Exposed Human Population', 'Interpreter', 'latex', 'FontSize', 13)
xlim([0, 20])

figure(3)
plot(t, I_h, 'r-', 'LineWidth',2)
xlabel('Time [t]', 'Interpreter', 'latex', 'FontSize', 13)
ylabel('Infected Human Population', 'Interpreter', 'latex', 'FontSize', 13)
xlim([0, 1000])

figure(4)
plot(t, H, 'r-', 'LineWidth',2)
xlabel('Time [t]', 'Interpreter', 'latex', 'FontSize', 13)
ylabel('Hospitilized Human Population', 'Interpreter', 'latex', 'FontSize', 13)
xlim([0, 90])

figure(5)
plot(t, R_h, 'r-', 'LineWidth',2)
xlabel('Time [t]', 'Interpreter', 'latex', 'FontSize', 13)
ylabel('Recovered Human Population', 'Interpreter', 'latex', 'FontSize', 13)
xlim([0, 1700])

figure(6)
plot(t, S_v, 'r-', 'LineWidth',2)
xlabel('Time [t]', 'Interpreter', 'latex', 'FontSize', 13)
ylabel('Susceptible Vector Population', 'Interpreter', 'latex', 'FontSize', 13)
xlim([0, 80])

figure(7)
plot(t, E_v, 'r-', 'LineWidth',2)
xlabel('Time [t]', 'Interpreter', 'latex', 'FontSize', 13)
ylabel('Exposed Vector Population', 'Interpreter', 'latex', 'FontSize', 13)
xlim([0, 170])

figure(8)
plot(t, I_v, 'r-', 'LineWidth',2)
xlabel('Time [t]', 'Interpreter', 'latex', 'FontSize', 13)
ylabel('Infected Vector Population', 'Interpreter', 'latex', 'FontSize', 13)
xlim([0, 700])