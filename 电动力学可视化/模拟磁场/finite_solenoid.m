%% 有限长螺线管磁场分布模拟
clc; clear; 

% 参数配置
L_sol = 1.0;          % 螺线管长度
R = 0.3;              % 螺线管半径
N_turns = 25;         % 总匝数
res = 0.02;           % 分辨率
N_show = 12;          % 可视化线圈数
gray_color = [0.4 0.4 0.4];

% 建立空间网格
[z_grid, r_grid] = meshgrid(-1.2:res:1.2, -0.8:res:0.8);
R_dist = abs(r_grid);

% 磁场分量叠加计算
Bz_tot = zeros(size(z_grid)); 
Br_tot = zeros(size(z_grid));
z_coils = linspace(-L_sol/2, L_sol/2, N_turns);

for i = 1:N_turns
    z_off = z_coils(i);
    alpha = sqrt(R^2 + R_dist.^2 + (z_grid-z_off).^2 - 2*R*R_dist + eps);
    beta  = sqrt(R^2 + R_dist.^2 + (z_grid-z_off).^2 + 2*R*R_dist + eps);
    k_sq = 4*R*R_dist ./ (beta.^2 + eps); 
    [K, E] = ellipke(k_sq);
    
    Bz_add = (1./(2*pi*beta)) .* (K + (R^2 - R_dist.^2 - (z_grid-z_off).^2)./(alpha.^2) .* E);
    Br_add = ((z_grid-z_off)./(2*pi*R_dist.*beta)) .* (-K + (R^2 + R_dist.^2 + (z_grid-z_off).^2)./(alpha.^2) .* E);
    
    Br_add(R_dist < 1e-6) = 0;
    Bz_tot = Bz_tot + Bz_add;
    Br_tot = Br_tot + Br_add;
end
Br_tot = sign(r_grid) .* Br_tot;

% 磁感应强度幅值处理
Bmag = sqrt(Bz_tot.^2 + Br_tot.^2); 
Bmag(Bmag > 15) = 15;

% 可视化绘图
figure('Color', 'w', 'Name', 'Solenoid Magnetic Field', 'Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.5]);
hold on;

% 绘制磁场强度填充背景
contourf(z_grid, r_grid, Bmag, 100, 'LineStyle', 'none'); 
colormap(jet);
hcb = colorbar;
ylabel(hcb, 'Magnetic Field Intensity |B|');
title('有限长螺线管磁场仿真');

% 绘制磁场流线
sy_seeds = [-1.1 * ones(1, 12), -1.1 * ones(1, 12)];
sz_seeds = [linspace(-0.7, 0.7, 12), linspace(-0.25, 0.25, 12)];
hL = streamline(z_grid, r_grid, Bz_tot, Br_tot, sy_seeds, sz_seeds);
set(hL, 'Color', [0.1 0.1 0.1], 'LineWidth', 0.8);

% 绘制归一化磁场矢量箭头
M = 5; 
u = Bz_tot(1:M:end, 1:M:end) ./ (Bmag(1:M:end, 1:M:end) + eps);
v = Br_tot(1:M:end, 1:M:end) ./ (Bmag(1:M:end, 1:M:end) + eps);
quiver(z_grid(1:M:end, 1:M:end), r_grid(1:M:end, 1:M:end), u, v, 0.4, 'r', 'LineWidth', 0.8);

% 绘制线圈实体标识
z_show = linspace(-L_sol/2, L_sol/2, N_show);
for i = 1:N_show
    zc = z_show(i);
    plot(zc, R, 'o', 'MarkerSize', 7, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', gray_color, 'LineWidth', 0.5);
    plot(zc, R, 'w.', 'MarkerSize', 2);
    plot(zc, -R, 'o', 'MarkerSize', 7, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', gray_color, 'LineWidth', 0.5);
    text(zc, -R, '×', 'Color', 'w', 'FontSize', 7, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end

% 绘制螺线管轮廓边框
rectangle('Position', [-L_sol/2, -R, L_sol, 2*R], 'EdgeColor', [0.2 0.2 0.2], 'LineStyle', '--', 'LineWidth', 0.8);

% 电流方向图例框
legend_x = 0.75; legend_y = 0.55; legend_w = 0.35; legend_h = 0.2;
rectangle('Position', [legend_x, legend_y, legend_w, legend_h], 'FaceColor', [1 1 1 0.8], 'EdgeColor', [0.5 0.5 0.5]);
plot(legend_x + 0.05, legend_y + 0.7*legend_h, 'o', 'MarkerSize', 8, 'MarkerFaceColor', gray_color);
plot(legend_x + 0.05, legend_y + 0.7*legend_h, 'w.', 'MarkerSize', 3);
text(legend_x + 0.1, legend_y + 0.7*legend_h, ' 电流流出 (Out)');
plot(legend_x + 0.05, legend_y + 0.3*legend_h, 'o', 'MarkerSize', 8, 'MarkerFaceColor', gray_color);
text(legend_x + 0.05, legend_y + 0.3*legend_h, '×', 'Color', 'w', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(legend_x + 0.1, legend_y + 0.3*legend_h, ' 电流流入 (In)');

% 坐标轴配置
axis equal;
axis([-1.2 1.2 -0.8 0.8]);
xlabel('z (轴向)'); ylabel('r (径向)');
set(gca, 'Layer', 'top', 'Box', 'on', 'TickDir', 'out');