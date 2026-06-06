% =========================================================================
% PEC_MATLAB_B - TRANSIENT SIMULATION STUDY & VORTEX SHEDDING PHENOMENA
% =========================================================================

rho_ar = 1.225;      
visc_ar = 1.5e-5;    
St_base = 0.2;       

D_topo = d_ext_mm / 1000; 
f1 = frequencias(1);      

V_crit = (f1 * D_topo) / St_base;
Re = (V_crit * D_topo) / visc_ar;

if Re < 3e5 
    St = 0.20; Cl_base = 0.6;
elseif Re >= 3e5 && Re < 3e6 
    St = 0.25; Cl_base = 0.2;
else 
    St = 0.28; Cl_base = 0.1;
end

if torre == 1 
    fator_rugosidade = 1.3; 
elseif torre == 2 
    fator_rugosidade = 1.15; 
else
    fator_rugosidade = 1.0;  
end

Cl = Cl_base * fator_rugosidade;

V_crit = (f1 * D_topo) / St;
omega_s = 2 * pi * f1; 

v1 = V(:, 1);
phi1 = v1 / v1(end-1); 

M_eq = phi1' * M_red * phi1; 
K_eq = (2 * pi * f1)^2 * M_eq; 

switch material
    case 'Aço'
        zeta = 0.02;      
        E = 210e9;       
        dens = 7850;     
    case 'Alumínio'
        zeta = 0.015;    
    case 'Betão'
        zeta = 0.05;     
    otherwise
        zeta = 0.02;     
end

C_eq = 2 * zeta * sqrt(M_eq * K_eq); 

F_L = 0.5 * rho_ar * V_crit^2 * D_topo * Cl;
F_modal = F_L; 

dt = 1 / (20 * f1); 
t_final = 30;       
t = 0:dt:t_final;
n_steps = length(t);

beta = 0.25;
gamma = 0.5;

u = zeros(1, n_steps); 
v = zeros(1, n_steps); 
a = zeros(1, n_steps); 

a(1) = (F_modal * sin(0) - C_eq*v(1) - K_eq*u(1)) / M_eq;

for i = 1:n_steps-1
    P_next = F_modal * sin(omega_s * t(i+1));
    
    c1 = 1 / (beta * dt^2);
    c2 = 1 / (beta * dt);
    c3 = (1 / (2 * beta)) - 1;
    
    K_hat = K_eq + c1 * M_eq + (gamma / (beta * dt)) * C_eq;
    
    term_M = M_eq * (c1 * u(i) + c2 * v(i) + c3 * a(i));
    term_C = C_eq * ((gamma / (beta * dt)) * u(i) + (gamma / beta - 1) * v(i) + dt * (gamma / (2 * beta) - 1) * a(i));
    P_hat = P_next + term_M + term_C;
    
    u(i+1) = P_hat / K_hat;
    
    a_next = c1 * (u(i+1) - u(i)) - c2 * v(i) - c3 * a(i);
    v_next = v(i) + (1 - gamma) * dt * a(i) + gamma * dt * a_next;
    
    a(i+1) = a_next;
    v(i+1) = v_next;
end

figure('Color', 'w', 'Name', 'Análise de Vórtices - Newmark');

subplot(2,1,1);
plot(t, u * 1000, 'b', 'LineWidth', 1.2); 
title(['Resposta Dinâmica no Topo (V_{crit} = ', num2str(V_crit, '%.2f'), ' m/s)']);
ylabel('Deslocamento (mm)');
grid on;

subplot(2,1,2);
plot(t, v, 'r');
title('Velocidade no Topo');
xlabel('Tempo (s)');
ylabel('Velocidade (m/s)');
grid on;

fprintf('\n=================================================');
fprintf('\n--- RESULTADOS DA ANÁLISE DINÂMICA (SDOF) ---');
fprintf('\nVelocidade Crítica: %.2f m/s', V_crit);
fprintf('\nReynolds: %.2e', Re);
fprintf('\nDeslocamento Máximo (Steady-State): %.2f mm\n', max(abs(u(round(end/2):end))) * 1000);

U = 15;             
rho_ar = 1.225;     

if torre == 1 && var_diam == 1
    D_referencia = diams_seccoes(end) / 1000; 
elseif torre == 2 && tipo_geo == 1
    D_referencia = d_topo_mm / 1000;          
else
    D_referencia = d_ext_mm / 1000;           
end

f1 = frequencias(1); 
V_crit = (f1 * D_referencia) / St_base;
Re = (V_crit * D_referencia) / visc_ar;

if Re < 3e5, St = 0.20; Cl = 0.6;
elseif Re >= 3e5 && Re < 3e6, St = 0.25; Cl = 0.2;
else, St = 0.28; Cl = 0.1; 
end

V_crit = (f1 * D_referencia) / St; 
fs = (St * U) / D_referencia;      

fprintf('\n--- RESULTADOS DA ANÁLISE TRANSIENTE (MDOF) ---\n');
fprintf('Frequência de Vortex Shedding (fs): %.3f Hz\n', fs);
fprintf('Relação fs / f1: %.2f', fs / res_f(end));
fprintf('\n=================================================\n');

dt = 0.005;                         
t_final = 30;                        
t = 0:dt:t_final;                   
n_steps = length(t);

omega1 = 2 * pi * frequencias(1);
omega2 = 2 * pi * frequencias(2);
alpha = zeta * (2 * omega1 * omega2) / (omega1 + omega2);
beta_rayleigh = zeta * 2 / (omega1 + omega2); 

C_red = alpha * M_red + beta_rayleigh * K_red;

F_amp = 0.5 * rho_ar * U^2 * (d_ext_mm/1000) * Cl * L_total; 
F_t = F_amp * sin(2 * pi * fs * t);

gamma = 0.5;
beta_n = 0.25; 
a0 = 1/(beta_n*dt^2); a1 = gamma/(beta_n*dt); a2 = 1/(beta_n*dt);
a3 = 1/(2*beta_n) - 1; a4 = gamma/beta_n - 1; a5 = dt/2*(gamma/beta_n - 2);

u = zeros(size(K_red, 1), n_steps);
v = zeros(size(K_red, 1), n_steps);
acc = zeros(size(K_red, 1), n_steps);

K_eff = K_red + a0*M_red + a1*C_red;

for i = 1:n_steps-1
    Force = zeros(size(K_red, 1), 1); 
    
    for n = 1:NUMEL
        idx_node = 2*n - 1; 
        z_perc = n / NUMEL;
        
        if torre == 1 
            idx_s = sum(z_perc > posicoes_flange) + 1;
            d_local = diams_seccoes(idx_s) / 1000;
        elseif torre == 2 && tipo_geo == 1 
            d_local = (d_ext_mm - (d_ext_mm - d_topo_mm) * z_perc) / 1000;
        else
            d_local = d_ext_mm / 1000;
        end
        
        F_local_amp = 0.5 * rho_ar * U^2 * d_local * Cl * (L_total/NUMEL);
        
        Force(idx_node) = F_local_amp * sin(2 * pi * fs * t(i+1));
    end 
    
    R_eff = Force + M_red*(a0*u(:,i) + a2*v(:,i) + a3*acc(:,i)) + ...
            C_red*(a1*u(:,i) + a4*v(:,i) + a5*acc(:,i));
    
    u(:,i+1) = K_eff \ R_eff;
    
    acc(:,i+1) = a0*(u(:,i+1) - u(:,i)) - a2*v(:,i) - a3*acc(:,i);
    v(:,i+1) = v(:,i) + (1-gamma)*dt*acc(:,i) + gamma*dt*acc(:,i+1);
end

figure('Color', 'w', 'Name', 'Força e Deslocamento no topo');
subplot(2,1,1);
plot(t, F_t, 'r'); title('Força de Sustentação no Topo (Vortex Shedding)');
ylabel('Força (N)'); grid on;

subplot(2,1,2);
plot(t, u(end-1, :)*1000); title('Deslocamento Transiente do Topo');
xlabel('Tempo (s)'); ylabel('Deslocamento (mm)'); grid on;

if ~exist('diams_seccoes', 'var')
    diams_seccoes = ones(1, NUMEL) * d_ext_mm;
end

save('dados_simulacao.mat', 'u', 'fs','zeta', 't_final','acc', 'v', 't', 'dt', 'L_total', 'NUMEL', 'E', 'thick_mm', 'diams_seccoes', 'K_red', 'M_red', 'C_red', 'zeta');
fprintf('-> Dados guardados em dados_simulacao.mat\n');