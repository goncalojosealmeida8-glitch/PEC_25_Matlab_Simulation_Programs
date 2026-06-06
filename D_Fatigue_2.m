% ==================================================================================
% PEC_MATLAB_D - CONTINUOUS STUDIES ON FATIGUE ASSESSMENT (OPERATIONAL SERVICE LIFE)
% ==================================================================================

clear; clc;

if exist('propriedades_torre.mat', 'file')
    load('propriedades_torre.mat');
else
    error('Corre primeiro o Script 1 para gerar as matrizes da torre!');
end

V_vetor = 1:1:25; 
danos_por_velocidade = zeros(length(V_vetor), 1);
sigma_max_por_velocidade = zeros(length(V_vetor), 1);

dt = 0.01; t_final = 30; t = 0:dt:t_final;
time_steps = length(t);
gamma = 0.5; beta = 0.25;
a0 = 1/(beta*dt^2); a1 = gamma/(beta*dt); a2 = 1/(beta*dt);
a3 = 1/(2*beta) - 1; a4 = gamma/beta - 1;
a5 = dt/2*(gamma/beta - 2); a6 = dt*(1-gamma); a7 = gamma*dt;
alpha_ray = 0.05; beta_ray = 0.001; 
C_red = alpha_ray * M_red + beta_ray * K_red;

K_eff = K_red + a0*M_red + a1*C_red;
fprintf('-> A iniciar varrimento para %d velocidades...\n', length(V_vetor));

for j = 1:length(V_vetor)
    U = V_vetor(j);
    
    fs = (0.2 * U) / (diams_seccoes(end)/1000); 
    
    u_sim = zeros(size(K_red,1), time_steps);
    v_sim = zeros(size(K_red,1), time_steps);
    acc_sim = zeros(size(K_red,1), time_steps);
    
    for i = 1:time_steps-1
        Force = zeros(size(K_red,1), 1);
        for e = 1:NUMEL
             idx_node = 2*e - 1;
             d_local = diams_seccoes(e)/1000;
             F_amp = 0.5 * 1.225 * U^2 * d_local * 0.6 * (L_total/NUMEL);
             Force(idx_node) = F_amp * sin(2 * pi * fs * t(i+1));
        end
        
        R_eff = Force + M_red*(a0*u_sim(:,i) + a2*v_sim(:,i) + a3*acc_sim(:,i)) + ...
                C_red*(a1*u_sim(:,i) + a4*v_sim(:,i) + a5*acc_sim(:,i));
        
        u_sim(:,i+1) = K_eff \ R_eff;
        
        acc_sim(:,i+1) = a0*(u_sim(:,i+1) - u_sim(:,i)) - a2*v_sim(:,i) - a3*acc_sim(:,i);
        v_sim(:,i+1) = v_sim(:,i) + a6*acc_sim(:,i) + a7*acc_sim(:,i+1);
    end
    
    idx_base = 1:4;
    u_el = u_sim(idx_base, :);
    curv = (1/(L_total/NUMEL)^2) * ([-6 + 12*0.5, (L_total/NUMEL)*(-4+6*0.5), 6-12*0.5, (L_total/NUMEL)*(-2+6*0.5)] * u_el);
    sigma_max_t = E_pa * curv * (diams_seccoes(1)/2000); 
    
    delta_sigma = max(sigma_max_t) - min(sigma_max_t);
    sigma_max_por_velocidade(j) = delta_sigma / 1e6; 
    
    if sigma_max_por_velocidade(j) > 1e-3
        NR = 2e6 * (80 / sigma_max_por_velocidade(j))^3;
        n_ciclos = fs * t_final;
        danos_por_velocidade(j) = n_ciclos / NR;
    end
end

k = 2; C = 8; 
prob_vento = (k/C) * (V_vetor/C).^(k-1) .* exp(-(V_vetor/C).^k);

fator_anual = (365 * 24 * 3600) / t_final;
dano_anual = sum(danos_por_velocidade .* prob_vento' * fator_anual);

vida_util_anos = 1 / dano_anual;

figure('Color', 'w');

subplot(1,2,1);
plot(V_vetor, sigma_max_por_velocidade, 'r-o', 'LineWidth', 1.5);
grid on; title('Sensibilidade à Velocidade');
xlabel('Velocidade do Vento (m/s)'); ylabel('\sigma_{max} (MPa)');

subplot(1,2,2);
N_curva = logspace(4, 8, 100);
S_curva = 80 * (2e6 ./ N_curva).^(1/3); 
loglog(N_curva, S_curva, 'k--'); hold on;
scatter(2e6 * (80./sigma_max_por_velocidade).^3, sigma_max_por_velocidade, 'filled');
title('Verificação Curva S-N'); xlabel('Ciclos (N)'); ylabel('\Delta\sigma (MPa)');
grid on;

fprintf('\nESTIMATIVA DE VIDA ÚTIL: %.1f Anos\n', vida_util_anos);
