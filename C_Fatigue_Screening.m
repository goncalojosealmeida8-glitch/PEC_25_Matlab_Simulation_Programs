% =========================================================================
% PEC_MATLAB_C - PRELIMINARY FATIGUE DEGRADATION SCREENING
% =========================================================================

clear; clc;

if exist('dados_simulacao.mat', 'file')
    load('dados_simulacao.mat');
    fprintf('-> Dados da simulação carregados com sucesso.\n');
else
    error('ERRO: O ficheiro "dados_simulacao.mat" não existe. Adiciona o comando "save" ao transient_study.m e corre-o primeiro.');
end

indices_trans = 1:2:size(u, 1);
u_nodes = u(indices_trans, :);
a_nodes = acc(indices_trans, :);
z_coord = linspace(0, L_total, size(u_nodes, 1));

H_plot = z_coord; 

fprintf('-> A calcular histórico de tensões por elemento...\n');

n_gdls = size(u, 1); 
n_elems = NUMEL;
le = L_total / n_elems;
time_steps = size(u, 2);
sigma_hist = zeros(n_elems, time_steps);

if length(diams_seccoes) ~= n_elems
    old_diams = diams_seccoes;
    diams_seccoes_full = interp1(linspace(0, 1, length(old_diams)), old_diams, linspace(0, 1, n_elems), 'nearest');
else
    diams_seccoes_full = diams_seccoes;
end

for t_idx = 1:time_steps
    for e = 1:n_elems
        idx = (2*e-1):(2*e+2);
        if idx(end) > n_gdls, break; end
        
        u_el = u(idx, t_idx);
        
        d2u_dz2 = (1/le^2) * ([-6 + 12*0.5, le*(-4+6*0.5), 6-12*0.5, le*(-2+6*0.5)] * u_el);
        
        d_local = diams_seccoes_full(e) / 1000; 
        
        sigma_hist(e, t_idx) = E * d2u_dz2 * (d_local / 2);
    end
end

I_elems = zeros(n_elems, 1);
M_hist = zeros(n_elems, time_steps);

if ~exist('thick_mm', 'var'), thick_mm = 10; end 

for e = 1:n_elems
    d_ext_m = diams_seccoes(min(e, length(diams_seccoes))) / 1000;
    r_ext_m = d_ext_m / 2;
    r_int_m = r_ext_m - (thick_mm / 1000); 
    I_elems(e) = (pi/4) * (r_ext_m^4 - r_int_m^4);
    
    idx = (2*e-1):(2*e+2);
    if idx(end) <= n_gdls
        for t_idx = 1:time_steps
            u_el = u(idx, t_idx);
            d2u_dz2 = (1/le^2) * ([-6 + 12*0.5, le*(-4+6*0.5), 6-12*0.5, le*(-2+6*0.5)] * u_el);
            M_hist(e, t_idx) = E * I_elems(e) * d2u_dz2;
        end
    end
end

figure('Color', 'w', 'Name', 'Dashboard de Integridade Estrutural', 'Units', 'normalized', 'Position', [0.15 0.1 0.7 0.8]);

subplot('Position', [0.10, 0.55, 0.36, 0.38]);
max_a = max(abs(a_nodes), [], 2);
plot(z_coord, max_a, 'k', 'LineWidth', 2);
title('Aceleração Máxima'); xlabel('Altura (m)'); ylabel('m/s^2');
grid on;

subplot('Position', [0.54, 0.55, 0.36, 0.38]);
sigma_max_perfil = max(abs(sigma_hist), [], 2) / 1e6; 
z_elems = linspace(le/2, L_total - le/2, n_elems);
plot(z_elems, sigma_max_perfil, 'b', 'LineWidth', 2);
hold on;
[sig_max_global, idx_s] = max(sigma_max_perfil);
plot(z_elems(idx_s), sig_max_global, 'ro', 'MarkerFaceColor', 'r');
title('Tensão Máxima (MPa)'); xlabel('Altura (m)'); ylabel('MPa');
grid on;

subplot('Position', [0.32, 0.08, 0.36, 0.38]);
M_max_plot = max(abs(M_hist), [], 2) / 1000; 
plot(z_elems, M_max_plot, 'm', 'LineWidth', 2);
hold on;
[M_max_global, idx_m] = max(M_max_plot);
plot(z_elems(idx_m), M_max_global, 'ko', 'MarkerFaceColor', 'k');
title('Momento Fletor (kN.m)'); xlabel('Altura (m)'); ylabel('kN.m');
grid on;

fprintf('\n--- RESULTADOS FINAIS ---\n');
fprintf('Aceleração Max: %.2f m/s^2\n', max(max_a));
fprintf('Tensão Max Flexão: %.2f MPa\n', sig_max_global);
fprintf('Momento Fletor Max: %.2f kN.m\n', M_max_global);
fprintf('-------------------------\n');

```
