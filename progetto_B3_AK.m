close all;
clear all;
clc;

%% 1.0, initialisation of symbols and necessary values
x           = [sym('x_1') ; sym('x_2')];
u           = sym('u');
beta        = sym('beta');
k           = sym('k');

% ch
l           = 4;
J_0         = sym('J_0');
J           = sym('J_', [1, l]);
phi         = sym('phi_', [1, l]);

J_x         = J_0;
for i = 1:l
    J_x     = J_x + J(i)*cos(i*x(1)+phi(i));
end

%% 1.1, functions dot_x = f(x, u) ; y = h(x, u)
f_xu        = [x(2); 1/J_x * (u - beta*x(2) - k*x(1))];
h_xu        = x(1);

disp("f(x, u) = ");
disp(f_xu);
disp("h(x, u) = ");
disp(h_xu);

% values at equilibrium
x_eq        = [17/36*pi; 0];
u_eq        = 85/9*pi;

disp("x_eq = ");
disp(x_eq);
disp("u_eq = ");
disp(u_eq); 

% system parametres
J_0_param   = 2;
J_param     = [0.5 0.01 0.01 0.08];
phi_param   = [-1 2 1.5 -2];
beta_param  = 0.2;
k_param     = 20;

symbols     = [reshape(x.', 1, []) u J_0 J phi beta k];
values      = [reshape(x_eq.', 1, []) u_eq J_0_param J_param phi_param beta_param k_param];

% values check at equilibrium
h_xu_eq     = subs(h_xu, symbols, values);
f_xu_eq     = subs(f_xu, symbols, values);

disp("h(x_eq, u_eq) = ");
disp(h_xu_eq);

disp("f(x_eq, u_eq) = ");
disp(f_xu_eq);

%% 1.2, linearised matrixes
A           = [diff(f_xu(1), x(1)) diff(f_xu(1), x(2)) ; diff(f_xu(2), x(1)) diff(f_xu(2), x(2))];
A_eq        = double(subs(A, symbols, values));
disp("A_eq = ");
disp(round(A_eq, 4));

B           = [diff(f_xu(1), u) ; diff(f_xu(2), u)];
B_eq        = double(subs(B, symbols, values));
disp("B_eq = ");
disp(round(B_eq, 4));

C           = [diff(h_xu, x(1)) diff(h_xu, x(2))];
C_eq        = double(subs(C, symbols, values));
disp("C_eq = ");
disp(round(C_eq, 4));

D           = [diff(h_xu, u)];
D_eq        = double(subs(D, symbols, values));
disp("D_eq = ");
disp(round(D_eq, 4));

%%  2, transfer function
eigs        = double(eig(A_eq));
disp("eigenvalues = ");
disp(round(eigs, 4));

s           = tf('s');
G           = C_eq*inv(s*eye(2) - A_eq)*B_eq + D_eq
system      = ss(A_eq, B_eq, C_eq, D_eq);
pp          = pole(G);
zz          = zero(G);
pzmap(G);
grid on;

disp("poles = ");
disp(pp);
disp("zeroes = ");
disp(zz);

%% 3.0, project specs
% 1)
w           = 3;
d           = 3;
e_star      = 0.05;

% 2)
pm_target   = 45;

% 3)
S_star     = 8;

% 4)
epsilon     = 1;
T_star     = 0.1;

% 5)
A_d         = 50;
omega_d_min = 1e-4;
omega_d_MAX = 0.1;

% 6)
A_n         = 60;
omega_n_min = 1.5 * 1e4;
omega_n_MAX = 1e6;

%% 3.1, static regulator 
mu_s_error  = (w + d) / e_star / abs(evalfr(G, j*0));

mu_s_d      = 10 ^ (A_d / 20) / abs(evalfr(G, j*omega_d_MAX));
R_s         = max(mu_s_error,mu_s_d);

% extended system
G_e         = R_s * G;

% Bode(G_e) with specs
figure(2);
hold on;

% overshoot_star => phase margin
xi_star     = abs(log(S_star/100)) / sqrt(pi^2 + log(S_star/100)^2);
pm          = max(xi_star * 100, pm_target);

% specs on d
Bnd_d_x     = [omega_d_min; omega_d_MAX; omega_d_MAX; omega_d_min];
Bnd_d_y     = [A_d; A_d; -150; -150];
patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);

% specs on n
Bnd_n_x     = [omega_n_min; omega_n_MAX; omega_n_MAX; omega_n_min];
Bnd_n_y     = [-A_n; -A_n; 100; 100];
patch(Bnd_n_x, Bnd_n_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);

% specs on settling time (minimum critical pulse)
omega_Ta_min    = 1e-4; % lower bound for plotting
omega_Ta_max    = 300/(pm * T_star); % omega_c >= 100/(pm*T^*)
Bnd_Ta_x        = [omega_Ta_min; omega_Ta_max; omega_Ta_max; omega_Ta_min];
Bnd_Ta_y        = [0; 0; -150; -150];
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);

% legend colours
Legend_mag      = ["A_d"; "A_n"; "\omega_{c,min}"; "G(j\omega)"];
legend(Legend_mag);

% Bode plot with stability margins
omega_plot_min  = 1e-4;
omega_plot_max  = 1e5;
margin(G_e,{omega_plot_min,omega_plot_max});
grid on; zoom on;

% overshoot specs (phase margin)
omega_c_min     = omega_Ta_max;
omega_c_max     = omega_n_min;

phi_up          = pm - 180;
phi_low         = -270; % lower bound for plot

Bnd_Mf_x        = [omega_c_min; omega_c_max; omega_c_max; omega_c_min];
Bnd_Mf_y        = [phi_up; phi_up; phi_low; phi_low];
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);

% legend colours
Legend_arg      = ["L(j\omega)"; "M_f"];
legend(Legend_arg);

%% 3.2, dynamic regulator
pm_star         = pm + 8;
omega_c_star    = 2500 / (T_star * pm);

mag_omega_c_star_dB = abs(evalfr(G_e, j*omega_c_star));
arg_omega_c_star    = rad2deg(angle(evalfr(G_e, j*omega_c_star)));

M_star          = 1/mag_omega_c_star_dB;
phi_star        = pm_star - 180 - arg_omega_c_star;

tau             = (M_star - cos(phi_star*pi/180))/(omega_c_star*sin(phi_star*pi/180));
alpha_tau       = (cos(phi_star*pi/180) - 1/M_star)/(omega_c_star*sin(phi_star*pi/180));
alpha           = alpha_tau / tau;

if min(tau,alpha) < 0
    fprintf('Errore: parametri rete anticipatrice negativi');
    return;
end

R_hf             = 1/(1 + s/3e3);
R_d              = (1 + tau * s)/(1 + alpha * tau * s)*R_hf;
R                = R_s * R_d;
L                = R * G; % funzione di anello

figure(3);
hold on;

% specs on amplitude
patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_n_x, Bnd_n_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);
legend(Legend_mag);

% Bode plot with stability margins
margin(L,{omega_plot_min,omega_plot_max});
grid on; zoom on;

% specs on phase
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);
legend(Legend_arg);

%% 4.1, performance check on closed loop

% complementary sensitivity function
F               = L/(1+L);

% step response
figure(4);

T_simulation    = 2*T_star;
[y_step,t_step] = step(w*F, T_simulation);
plot(t_step,y_step,'b');
grid on, zoom on, hold on;

LV = evalfr(w*F,0);

% overshoot constraint
patch([0,T_simulation,T_simulation,0],[LV*(1+S_star/100),LV*(1+S_star/100),LV*2,LV*2],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);

% 1% settling time constraint
patch([T_star,T_simulation,T_simulation,T_star],[LV*(1-0.01),LV*(1-0.01),0,0],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
patch([T_star,T_simulation,T_simulation,T_star],[LV*(1+0.01),LV*(1+0.01),LV*2,LV*2],'g','FaceAlpha',0.1,'EdgeAlpha',0.1);

ylim([0,LV*2]);

Legend_step = ["Risposta al gradino"; "Vincolo sovraelongazione"; "Vincolo tempo di assestamento"];
legend(Legend_step);

%% 4.2, output noise check

% sensitivity function
S               = 1/(1+L);
figure(5);

% output noise simulation
tt              = 0:1e-2:2e2;
dd              = 0;
for k = 1:3
    dd          = dd + d*sin(0.01*k*tt);
end
y_d             = lsim(S,dd,tt);
hold on, grid on, zoom on
plot(tt,dd,'m')
plot(tt,y_d,'b') 
grid on
legend('d(t)','y_d(t)')

%% 4.3, input noise check

% complementary sensitivity function 
figure(6);

% input noise simulation
tt              = 0:1e-5:2*1e-3;
nn              = 0;
for k = 1:3
    nn          = nn + 3*sin(1e5*k*tt);
end
y_n             = lsim(F,nn,tt);
hold on, grid on, zoom on
plot(tt,nn,'m')
plot(tt,y_n,'b')
grid on
legend('n(t)','y_n(t')
