%% SCRIPT : autopilote PID avec action intégrale (boucle for t = t0:dt:tf)
clear all; close all; clc;

%% constantes
Re = 6371210;
g0 = 9.80665;

%% temps
dt = 0.01; t0 = 0; tf = 30;

%% historiques
u_hist = []; v_hist = []; w_hist = [];
phi_hist = []; theta_hist = []; psi_hist = [];
p_hist = []; q_hist = []; r_hist = [];
t_hist = [];
x_hist = []; y_hist = []; z_hist = []; % historiques position

%% paramètres système
m0 = 100; 
m_fuel= 0.25*m0;
mc = 0.5;
tc = m_fuel/mc;
A = 50; B = 300; C = 300;   % Ixx, Iyy, Izz
S = 2; dm = 0.15;

%% coefficients aérodynamiques (ton modèle)
Cl_deltaA = 0.5; Cl_p = 0.5;
Cm_alpha = 0.5; Cm_deltaE = 0.5; Cm_q = 0.5;
Cn_beta = 0.5; Cn_deltaR = 0.5; Cn_r = 0.5;
Cx0 = 0.5;
Cy_beta = 0.5; Cy_deltaR = 0.5; Cy_r = 0.5;
Cz_alpha = 0.5; Cz_deltaE = 0.5; Cz_q = 0.5;
Vm = 200;

eps_val = 1e-8;      % protection numérique

%% conditions initiales
u = 0; v = 0; w = 0;
p = 0; q = 0; r = 0;
phi = 0; theta = 0; psi = 0;

%% position initiale
x = 0; y = 0; z = 0;

%% autopilote PID+I : gains initiaux (point de départ)
% Roulis (phi)
kp_phi = 12;
ki_phi = 0;    % intégrale roulis (à réduire si instabilité)
kd_phi = 0.5;

% Tangage (theta)
kp_theta = 12;
ki_theta = 0;
kd_theta = 9;

% Lacet (psi)
kp_psi = 12;
ki_psi = 0;
kd_psi = 6;

% saturations
delta_max = deg2rad(30);


%% consignes (modifiable)
phi_r = 0;      % consigne roulis
theta_r =- pi/6;           % consigne tangage
psi_r = pi/12;            % consigne lacet

%% boucle temporelle (sans index i)
for t = t0 : dt : tf

 
    %% gravité + densité
    g = g0 * (Re / (Re + h))^2;
    rho = 1.222 * (1 - 0.0000225 * h)^(4.256);
    Q = 0.5 * rho * Vm^2;

  %% alpha / beta
    alpha = w/( u + eps_val);
    beta  = v/(u + eps_val);

    %% ---------- Contrôleurs PID+I (loi) ----------
    % erreurs
   

    % commande non-saturée (PID)

    %% ---------- Aérodynamique & propulsion ----------
    Cx = Cx0;
    Cy = Cy_beta * beta + Cy_deltaR * deltaR + Cy_r * r * dm / (2*Vm);



    CL = Cl_deltaA * deltaA + Cl_p * p * dm / (2*Vm);
    CM = Cm_alpha * alpha + Cm_deltaE * deltaE + Cm_q * q * dm / (2*Vm);



    % poussée & masse
    if (t < tc)
        Px = 100;
        m = m0 - mc * t;
    else
        Px = 0;
        m = m0 - mc * tc;
    end
% force et moments
    Fx = Q * S * Cx + Px;
  

    L  = Q * S * dm * CL;
  

    %% ---------- Dynamique de translation ----------
      %% ---------- Dynamique position ----------
    % matrice rotation corps -> inertiel (ZYX)
    R = [cos(theta)*cos(psi), sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi), cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi);
         cos(theta)*sin(psi), sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi), cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi);
         -sin(theta),         sin(phi)*cos(theta),                              cos(phi)*cos(theta)];

    V_inert = R * [u; v; w];

    x = x + V_inert(1) * dt;
    y = y + V_inert(2) * dt;
    z = z + V_inert(3) * dt;
    
    u_dot = g * sin(theta) + Fx/m - q*w + v*r;
    v_dot = -g * sin(phi) * cos(theta) + Fy/m - r*u + p*w;
    w_dot = -g * cos(theta) * cos(phi) + Fz/m - p*v + q*u;

    u = u + u_dot * dt;
  

    %% ---------- Cinématique Euler ----------
    cphi = cos(phi); sphi = sin(phi);
    ctheta = cos(theta);
    if abs(ctheta) < eps_val
        ctheta = eps_val * sign(ctheta + eps_val);
    end

    phi_dot   = p + sphi * tan(theta) * q + cphi * tan(theta) * r;
    theta_dot = cphi * q - sphi * r;
    psi_dot   = (r * cphi + q * sphi) / ctheta;

    %% ---------- Dynamique rotation ----------
    p_dot = ((B - C) * q * r + L) / A;


    p = p + p_dot * dt;
 

    phi   = phi + phi_dot * dt;
    
  

    %% stockage historiques
    u_hist = [u_hist u];
 

    phi_hist = [phi_hist phi];
 

    p_hist = [p_hist p];

   
    x_hist = [x_hist x];

    t_hist = [t_hist t];

end

%% Affichage
figure;
plot(t_hist,u_hist,'r',t_hist,v_hist,'g',t_hist,w_hist,'b','LineWidth',1.5);
xlabel('Temps (s)'); ylabel('m/s'); grid on; legend('u','v','w'); title('Vitesses corps');

figure;
plot(t_hist,p_hist,'r',t_hist,q_hist,'g',t_hist,r_hist,'b','LineWidth',1.5);
xlabel('Temps (s)'); ylabel('rad/s'); grid on; legend('p','q','r'); title('Vitesses angulaires');

figure;
plot(t_hist,phi_hist,'r',t_hist,theta_hist,'g',t_hist,psi_hist,'b','LineWidth',1.5);
xlabel('Temps (s)'); ylabel('rad'); grid on; legend('\phi','\theta','\psi'); title('Angles Euler');

figure;
plot3(x_hist, y_hist, z_hist,'b','LineWidth',1.5); grid on; axis equal;
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title('Trajectoire du missile');
