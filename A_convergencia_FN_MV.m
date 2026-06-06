%PEC_A - Análise Modal e Convergência de Elementos Finitos
clear; clc;

L_total   = input('Total Lenght (m): ');      
d_ext_mm  = input('External diameter (mm): ');
thick_mm  = input('Thickness (mm): ');
n_sim     = input('Número de simulações (n): ');
numel_base = input('Nº de Elementos inicial: ');
numel_inc  = input('Incremento de elementos: ');

disp(' ')
disp('Material:')
disp('-----------------------------------------')
disp('0 - Stainless Steel')
disp('1 - Aluminium')
disp('2 - PVC')
disp('3 for Other')
disp('-----------------------------------------')
material = input(' < ');
   
if isempty(material) || material < 0 || material > 3
    disp('--------------------------')
    disp('ERROR IN INPUT ')
    disp('Material has been set to 0 (Steel)')
    disp('--------------------------')
    material = 0;
end

if material == 0
    dens  = 7850;        
    Y_MOD = 210;         
    n2 = 'Aço';
elseif material == 1
    dens  = 2700;        
    Y_MOD = 70;          
    n2 = 'Alumínio';
elseif material == 2
    dens  = 1300;          
    Y_MOD = 2.41;          
    n2 = 'PVC';
else
    dens  = input('Material density (kg/m3): ');
    Y_MOD = input('Young Modulus (GPa): ');
    n2 = input('Nome do Material: ');
end
 
E = Y_MOD * 1e9; 

disp(' ')
disp('Tipo de Torre:')
disp('-----------------------------------------')
disp('0 - Lisa')
disp('1 - Flangeada')
disp('2 - Encaixe')
disp('-----------------------------------------')
torre = input(' < ');

posicoes_flange = []; 
massa_f = 0; 
J_flange = 0;
n_encaixes = 0;
n1 = ''; n3 = ''; n4 = '';

if isempty(torre) || torre < 0 || torre > 2
    disp('--------------------------')
    disp('ERROR IN INPUT ')
    disp('Torre has been set to 0 (Lisa)')
    disp('--------------------------')
    torre = 0;
end

if torre == 0
    posicoes_flange = [];
    n1 = 'Lisa';

elseif torre == 1
    n1 = 'Flangeada';
    disp(' ')
    n_flanges = input('\n Quantas flanges deseja adicionar? ');
    d_ext_f = input('Diâmetro Externo das Flanges (mm): ');
    massa_f = input('Massa de cada flange (Kg): ');

    posicoes_flange = zeros(1, n_flanges);
    for p = 1:n_flanges
        posicoes_flange(p) = input(['Posição da flange ' num2str(p) ' (0 a 1 da altura): ']);
    end
    
    disp(' ')
    disp('Haverá variação de diâmetro entre as secções?')
    disp('0 - Não (Diâmetro constante em toda a torre)')
    disp('1 - Sim (Definir diâmetro para cada secção)')
    var_diam = input(' < ');
    
    if var_diam == 1
        n3 = 'com variação de diâmetro entre secções';
        diams_seccoes = zeros(1, n_flanges + 1);
        
        diams_seccoes(1) = d_ext_mm; 
        fprintf('Diâmetro da secção 1 (Base): %.0f mm (definido no início)\n', d_ext_mm);
        
        for s = 2:n_flanges + 1
            diams_seccoes(s) = input(['Diâmetro da secção ' num2str(s) ' (mm): ']);
        end
    else
        n3 = 'sem variação de diâmetro entre secções';
        diams_seccoes = ones(1, n_flanges + 1) * d_ext_mm;
    end
    
    r_ext_f = (d_ext_f / 2) / 1000;
    r_int_f = ((diams_seccoes(1)/2) - thick_mm) / 1000; 
    J_flange = 0.5 * massa_f * (r_ext_f^2 + r_int_f^2);

    if isempty(posicoes_flange)
        posicoes_flange = []; 
        fprintf('Nenhuma flange adicionada.\n');
    end
    

elseif torre == 2
    n1 = 'de Encaixe';
    
    disp(' ')
    disp('A geometria base da torre é:')
    disp('0 - Telescópica (Diâmetro constante)')
    disp('1 - Cónica (Variação linear de diâmetro)')
    tipo_geo = input(' < ');
    
    if tipo_geo == 1
        n4 = 'Telescópica (Diâmetro constante)';
        d_topo_mm = input('Diâmetro no Topo (mm): ');
    else
        n4 = 'Cónica (Variação linear de diâmetro)';
    end

    n_encaixes = input('Quantos encaixes deseja adicionar? ');
    pos_encaixe_perc = zeros(1, n_encaixes);
    L_overlap = zeros(1, n_encaixes);
    thick_int_mm = zeros(1, n_encaixes);
    
    for j = 1:n_encaixes
        disp(' ')
        fprintf('\n--- Dados do Encaixe nº%d ---\n', j);
        pos_encaixe_perc(j) = input('Posição do início do encaixe (0 a 1 da altura): ');
        L_overlap(j) = input('Comprimento da sobreposição/overlap (m): ');
        thick_int_mm(j) = input('Espessura do tubo interno no encaixe (mm): ');
    end
end

r_ext = (d_ext_mm / 1000) / 2;
r_int = r_ext - (thick_mm / 1000);
Area = pi * (r_ext^2 - r_int^2);
I_moment = (pi/4) * (r_ext^4 - r_int^4); 
rho_a = dens * Area;

r_int_f = r_ext; 

res_f = zeros(n_sim, 1); 
res_n = zeros(n_sim, 1);


for k = 1:n_sim
    NUMEL = numel_base + (k-1)*numel_inc;
    res_n(k) = NUMEL;
    le = L_total / NUMEL;
    
    NDOF = 2 * (NUMEL + 1);
    K_global = zeros(NDOF, NDOF);
    M_global = zeros(NDOF, NDOF);
    
    for i = 1:NUMEL
        z_mid = (i - 0.5) * le;
        z_perc = z_mid / L_total; 
        
        if torre == 1
            idx_s = sum(z_perc > posicoes_flange) + 1;
            d_atual_mm = diams_seccoes(idx_s);
            
        elseif torre == 2 && tipo_geo == 1
            d_atual_mm = d_ext_mm - (d_ext_mm - d_topo_mm) * z_perc;
            
        else
            d_atual_mm = d_ext_mm; 
        end
        
        r_ext_at = (d_atual_mm / 1000) / 2;
        r_int_at = r_ext_at - (thick_mm / 1000);
        I_atual = (pi/4) * (r_ext_at^4 - r_int_at^4);
        Area_atual = pi * (r_ext_at^2 - r_int_at^2);
        
        if torre == 2
            for j = 1:n_encaixes
                z_enc_ini = pos_encaixe_perc(j) * L_total;
                z_enc_fim = z_enc_ini + L_overlap(j);
                
                if z_mid >= z_enc_ini && z_mid <= z_enc_fim
                    r_ext_in = r_int_at; 
                    r_int_in = r_ext_in - (thick_int_mm(j) / 1000);
                    
                    I_atual = I_atual + (pi/4) * (r_ext_in^4 - r_int_in^4);
                    Area_atual = Area_atual + pi * (r_ext_in^2 - r_int_in^2);
                    break; 
                end
            end
        end
        
        rho_a_atual = dens * Area_atual;
        ke = (E * I_atual / le^3) * [12, 6*le, -12, 6*le; 6*le, 4*le^2, -6*le, 2*le^2; -12, -6*le, 12, -6*le; 6*le, 2*le^2, -6*le, 4*le^2];
        me = (rho_a_atual * le / 420) * [156, 22*le, 54, -13*le; 22*le, 4*le^2, 13*le, -3*le^2; 54, 13*le, 156, -22*le; -13*le, -3*le^2, -22*le, 4*le^2];
        
        idx = (2*i-1):(2*i+2);
        K_global(idx, idx) = K_global(idx, idx) + ke;
        M_global(idx, idx) = M_global(idx, idx) + me;
    end

    for p = 1:length(posicoes_flange)
        no_flange = round(posicoes_flange(p) * NUMEL) + 1;
        idx_trans = 2 * no_flange - 1;
        idx_rot = 2 * no_flange;
        M_global(idx_trans, idx_trans) = M_global(idx_trans, idx_trans) + massa_f;
        M_global(idx_rot, idx_rot) = M_global(idx_rot, idx_rot) + J_flange;
    end

    K_red = K_global(3:end, 3:end);
    M_red = M_global(3:end, 3:end);
    [V, D] = eig(K_red, M_red);
    frequencias = sort(sqrt(diag(D))) / (2*pi);
    
    res_f(k) = frequencias(1); 
end

if torre == 0
    tol = 1e-7; 
    min_sim_analise = 5;
    
elseif torre == 1
    if var_diam == 1
        tol = 1e-2;
        min_sim_analise = 15; 
    else
        tol = 1e-4;
        min_sim_analise = 10;
    end
    
elseif torre == 2
    if tipo_geo == 1
        tol = 1e-5;
        min_sim_analise = 20; 
    else
        tol = 1e-3;
        min_sim_analise = 20;
    end
end

idx_conv = n_sim; 
for k = 2:n_sim
    variacao = abs(res_f(k) - res_f(k-1)) / res_f(k-1);
    if variacao < tol && k >= min_sim_analise
        idx_conv = k;
        break;
    end
end

for k = 2:n_sim
    variacao = abs(res_f(k) - res_f(k-1)) / res_f(k-1);
    if variacao < tol && k > min_sim_analise
        idx_conv = k;
        break;
    end
end
numel_otimo = res_n(idx_conv);

fprintf('\n=================================================');
fprintf('\n   RESULTADOS DA ANÁLISE MODAL (Ponto Ótimo)');
fprintf('\n=================================================');
fprintf('\n Torre Selecionada: %s %s %s', n1, n3, n4);
fprintf('\n Material Selecionado: %s', n2);
fprintf('\n-------------------------------------------------');
for m = 1:3
    fprintf('\n Frequência Natural Modo %d: %.4f Hz', m, frequencias(m));
end
fprintf('\n=================================================\n');

n_nos = NUMEL + 1;
NDOF_total = 2 * n_nos;
V_completo = zeros(NDOF_total, size(V, 2));
V_completo(3:end, :) = V; 

indices_trans = 1:2:NDOF_total; 
H_plot = linspace(0, L_total, n_nos);

figure('Color', 'w', 'Name', 'Convergência e Modos de Vibração');

subplot(2, 3, 1:3); 
plot(res_n, res_f, '-b.', 'LineWidth', 1); hold on;
plot(numel_otimo, res_f(idx_conv), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
grid on;
title(['Estudo de Convergência (Ponto Ótimo: ' num2str(numel_otimo) ' elementos)']);
xlabel('Número de Elementos'); 
ylabel('1ª Freq. Natural (Hz)');
legend('Simulações', 'Ponto de Convergência');

for m = 1:3
    subplot(2, 3, m + 3);
    modo_lateral = V_completo(indices_trans, m);
    modo_norm = modo_lateral / max(abs(modo_lateral));
    
    plot(modo_norm, H_plot, 'b', 'LineWidth', 2); hold on;
    plot(-modo_norm, H_plot, 'r--', 'LineWidth', 1);
    plot([0 0], [0 L_total], 'k:'); 
    
    title(['Modo ' num2str(m) ' (' num2str(frequencias(m), 4) ' Hz)']);
    xlabel('Amplitude Rel.');
    if m == 1, ylabel('Altura (m)'); end
    grid on;
    axis([-1.2 1.2 0 L_total]);
end

fprintf('\n--- A iniciar exportação para Excel ---\n');

nomeArquivoExcel = 'Resultados_Frequencias.xlsx';
nomeFolha = 'Comparação de Torres';

if ~exist('frequencias', 'var') || isempty(frequencias)
    error('ERRO: A variável "frequencias" não existe ou está vazia!');
end

try
    nome_exibicao = strtrim(sprintf('%s %s %s', n1, n3, n4));
    
    r_ext_m = (d_ext_mm / 2) / 1000;
    r_int_m = r_ext_m - (thick_mm / 1000);
    area_m2 = pi * (r_ext_m^2 - r_int_m^2);
    massa_tubo = area_m2 * L_total * dens;
    massa_final = massa_tubo + (length(posicoes_flange) * massa_f);

    novaLinha = {nome_exibicao, n2, frequencias(1), frequencias(2), frequencias(3), massa_final};
    cabecalhos = {'Tipo_de_Torre', 'Material', 'Modo_1_Hz', 'Modo_2_Hz', 'Modo_3_Hz', 'Massa_Total_kg'};
    
    disp('-> Dados preparados com sucesso.');
catch ME
    fprintf('ERRO ao preparar dados: %s\n', ME.message);
end

try
    if exist(nomeArquivoExcel, 'file')
        disp('-> O ficheiro já existe. A tentar anexar...');
        opts = detectImportOptions(nomeArquivoExcel, 'Sheet', nomeFolha);
        tabelaAtual = readtable(nomeArquivoExcel, opts, 'Sheet', nomeFolha);
        
        novaTabela = cell2table(novaLinha, 'VariableNames', tabelaAtual.Properties.VariableNames);
        tabelaFinal = [tabelaAtual; novaTabela];
        
        writetable(tabelaFinal, nomeArquivoExcel, 'Sheet', nomeFolha);
    else
        disp('-> O ficheiro não existe. A criar novo...');
        T = cell2table(novaLinha, 'VariableNames', cabecalhos);
        writetable(T, nomeArquivoExcel, 'Sheet', nomeFolha);
    end
    fprintf('-> SUCESSO: Dados gravados em %s\n', fullfile(pwd, nomeArquivoExcel));
catch ME
    fprintf('ERRO CRÍTICO na gravação: %s\n', ME.message);
    if contains(ME.message, 'permission', 'IgnoreCase', true)
        disp('Sugerido: FECHA O EXCEL! O ficheiro está bloqueado.');
    end
end

if ~exist('diams_seccoes', 'var')
    diams_seccoes = repmat(d_ext_mm, 1, NUMEL); 
end

E_pa = Y_MOD * 1e9; 

save('propriedades_torre.mat', 'M_red', 'K_red', 'frequencias', 'L_total', 'thick_mm', 'diams_seccoes', 'E_pa', 'dens', 'NUMEL');