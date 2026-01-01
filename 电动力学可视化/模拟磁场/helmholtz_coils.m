%% 亥姆霍兹线圈模拟 
clc; clear; 

% 参数配置
R = 0.5;            % 线圈半径
d = R;              % 线圈间距(亥姆霍兹条件 d=R)
res = 0.02;         % 网格分辨率
gray_color = [0.4 0.4 0.4];

% 建立空间网格(z轴向, r径向)
[z_grid, r_grid] = meshgrid(-1.0:res:1.0, -0.8:res:0.8);
R_dist = abs(r_grid); 

% 计算双线圈磁场分量叠加
[Bz1, Br1] = get_loop_field(z_grid, R_dist, R, -d/2);
[Bz2, Br2] = get_loop_field(z_grid, R_dist, R, d/2);

Bz_tot = Bz1 + Bz2;
Br_tot = sign(r_grid) .* (Br1 + Br2); 

% 磁感应强度幅值处理
Bmag = sqrt(Bz_tot.^2 + Br_tot.^2);
Bmag(Bmag > 10) = 10; 

% 可视化绘图
figure('Color', 'w', 'Name', 'Helmholtz Coils'); 
clf; hold on; axis equal;

% 绘制磁场强度填充背景
contourf(z_grid, r_grid, Bmag, 100, 'LineStyle', 'none'); 
colormap(jet);
hcb = colorbar;
ylabel(hcb, '磁感应强度 |B|');

% 绘制磁场流线
sy_seeds = -0.9 * ones(1, 16); 
sz_seeds = linspace(-0.5, 0.5, 16); 
hL = streamline(z_grid, r_grid, Bz_tot, Br_tot, sy_seeds, sz_seeds);
set(hL, 'Color', [0.1 0.1 0.1], 'LineWidth', 0.8);

% 绘制归一化红色磁场矢量箭头
M = 4; 
u = Bz_tot(1:M:end, 1:M:end) ./ (Bmag(1:M:end, 1:M:end) + eps);
v = Br_tot(1:M:end, 1:M:end) ./ (Bmag(1:M:end, 1:M:end) + eps);
quiver(z_grid(1:M:end, 1:M:end), r_grid(1:M:end, 1:M:end), u, v, 0.4, 'r', 'LineWidth', 0.8);

% 绘制线圈截面与电流方向标识
coil_z = [-d/2, d/2];
for cz = coil_z
    plot(cz, R, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', gray_color, 'LineWidth', 1);
    plot(cz, R, 'w.', 'MarkerSize', 4);
    plot(cz, -R, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', gray_color, 'LineWidth', 1);
    text(cz, -R, '×', 'Color', 'w', 'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end

% 电流方向图例框
legend_x = 0.65; legend_y = 0.65;
legend_w = 0.32; legend_h = 0.12;
rectangle('Position', [legend_x, legend_y, legend_w, legend_h], 'FaceColor', [1 1 1 0.9], 'EdgeColor', [0.6 0.6 0.6], 'LineWidth', 0.5);

plot(legend_x + 0.04, legend_y + 0.7*legend_h, 'o', 'MarkerSize', 7, 'MarkerFaceColor', gray_color, 'MarkerEdgeColor', 'k');
plot(legend_x + 0.04, legend_y + 0.7*legend_h, 'w.', 'MarkerSize', 3);
text(legend_x + 0.08, legend_y + 0.7*legend_h, '电流流出 (Out)', 'FontSize', 8.5);

plot(legend_x + 0.04, legend_y + 0.3*legend_h, 'o', 'MarkerSize', 7, 'MarkerFaceColor', gray_color, 'MarkerEdgeColor', 'k');
text(legend_x + 0.04, legend_y + 0.3*legend_h, '×', 'Color', 'w', 'FontSize', 7, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(legend_x + 0.08, legend_y + 0.3*legend_h, '电流流入 (In)', 'FontSize', 8.5);

% 坐标轴与标题配置
title('亥姆霍兹线圈磁场仿真 (轴对称截面)');
xlabel('z (轴向)'); ylabel('r (径向)');
axis([-1.0 1.0 -0.8 0.8]);
set(gca, 'Layer', 'top', 'Box', 'on', 'TickDir', 'out');

%% 单匝线圈磁场计算函数
function [Bz, Br] = get_loop_field(z_grid, R_dist, R, z0)
    alpha = sqrt(R^2 + R_dist.^2 + (z_grid-z0).^2 - 2*R*R_dist + eps);
    beta  = sqrt(R^2 + R_dist.^2 + (z_grid-z0).^2 + 2*R*R_dist + eps);
    k_sq = 4*R*R_dist ./ (beta.^2 + eps); 
    [K, E] = ellipke(k_sq);
    
    Bz = (1./(2*pi*beta)) .* (K + (R^2 - R_dist.^2 - (z_grid-z0).^2)./(alpha.^2) .* E);
    Br = ((z_grid-z0)./(2*pi*R_dist.*beta)) .* (-K + (R^2 + R_dist.^2 + (z_grid-z0).^2)./(alpha.^2) .* E);
    Br(R_dist < 1e-6) = 0;
end