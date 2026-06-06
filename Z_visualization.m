% =========================================================================
% PEC_MATLAB_E - 3D TOWER REPRESENTATION (MONOPOLE / SLEEP-JOINT / FLANGED)
% =========================================================================

figure('Name', 'Torre 3D Universal', 'Color', 'w');
hold on; view(3); camlight; lighting gouraud;

z_nodes = linspace(0, L_total, NUMEL + 1);
raios_m = diams_seccoes / 2000; 

if length(diams_seccoes) ~= NUMEL
    old_diams = diams_seccoes;
    diams_seccoes = interp1(linspace(0,1,length(old_diams)), old_diams, linspace(0,1,NUMEL), 'nearest');
end

raios_nos = zeros(1, NUMEL + 1);

raios_nos(1:NUMEL) = diams_seccoes(:) / 2000;

raios_nos(NUMEL + 1) = raios_nos(NUMEL);

for e = 1:NUMEL
    r_b = raios_nos(e);     
    r_t = raios_nos(e+1);   
    
    [X, Y, Z] = cylinder([r_b, r_t], 30);
    
    z_base = z_nodes(e);
    z_topo = z_nodes(e+1);
    Z = Z * (z_topo - z_base) + z_base;
    
    surf(X, Y, Z, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    
    if e > 1 && (diams_seccoes(e) > diams_seccoes(e-1) * 1.02)
        [Xf, Yf, Zf] = cylinder(r_b * 1.5, 30);
        Zf = Zf * 0.05 + z_base; 
        surf(Xf, Yf, Zf, 'FaceColor', [0.2 0.2 0.2], 'EdgeColor', 'none');
    end
end

xlabel('X (m)'); ylabel('Y (m)'); zlabel('Altura (m)');
title('Representação 3D Final');

pbaspect([1 1 3]); 
xlim([-2, 2]);     
ylim([-2, 2]);
zlim([0, L_total]);
axis vis3d;        
rotate3d on;
grid on;
h_light = camlight('headlight');
set(gcf, 'WindowButtonMotionFcn', @(src, event) camlight(h_light, 'headlight'));