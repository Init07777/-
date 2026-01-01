%% 圆电流环磁场模拟 
clc; clear; close all;

% 物理参数与全局配置
R = 0.6; I = 1; res = 0.02;
gray_color = [0.4 0.4 0.4];

% 建立网格并计算径向距离
[y_grid, z_grid] = meshgrid(-1.5:res:1.5, -1.2:res:1.2); 
r_dist = abs(z_grid); 

% 磁场核心计算（椭圆积分法）
alpha = sqrt(R^2 + r_dist.^2 + y_grid.^2 - 2*R*r_dist + eps);
beta  = sqrt(R^2 + r_dist.^2 + y_grid.^2 + 2*R*r_dist + eps);
k_sq = 4*R*r_dist ./ (beta.^2 + eps);
[K, E] = ellipke(k_sq);

By_plot = (I./(2*pi*beta)) .* (K + (R^2 - r_dist.^2 - y_grid.^2)./(alpha.^2 + eps) .* E);
Br_val = (I.*y_grid./(2*pi*r_dist.*beta + eps)) .* (-K + (R^2 + r_dist.^2 + y_grid.^2)./(alpha.^2 + eps) .* E);
Bz_plot = sign(z_grid) .* Br_val; 

% 磁感应强度幅值计算与截断
Bmag = sqrt(By_plot.^2 + Bz_plot.^2); 
Bmag(Bmag > 4) = 4; 

% 创建画布并初始化绘图
figure('Color', 'w', 'Units', 'normalized', 'Position', [0.1, 0.1, 0.5, 0.5]);
clf; hold on;

% 绘制磁感应强度填充背景与色标
contourf(y_grid, z_grid, Bmag, 100, 'LineStyle', 'none'); 
colormap(jet);
hcb = colorbar; 
ylabel(hcb, '磁感应强度 |B|');

% 绘制磁场流线
hL = streamline(y_grid, z_grid, By_plot, Bz_plot, -1.4*ones(1,15), linspace(-1,1,15));
set(hL, 'Color', [0.1 0.1 0.1], 'LineWidth', 0.8);

% 绘制归一化磁场矢量箭头
M = 6;
u = By_plot(1:M:end, 1:M:end) ./ (Bmag(1:M:end, 1:M:end) + eps);
v = Bz_plot(1:M:end, 1:M:end) ./ (Bmag(1:M:end, 1:M:end) + eps);
quiver(y_grid(1:M:end, 1:M:end), z_grid(1:M:end, 1:M:end), u, v, 0.4, 'r', 'LineWidth', 0.8);

% 绘制线圈截面标识
plot(0, R, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', gray_color, 'LineWidth', 1);
plot(0, R, 'w.', 'MarkerSize', 4);
plot(0, -R, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', gray_color, 'LineWidth', 1);
text(0, -R, '×', 'Color', 'w', 'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

% 绘制紧凑图例框
leg_w = 0.75; leg_h = 0.35; leg_x0 = 0.65; leg_y0 = 0.75;
patch('XData', [leg_x0, leg_x0+leg_w, leg_x0+leg_w, leg_x0], ...
      'YData', [leg_y0, leg_y0, leg_y0+leg_h, leg_y0+leg_h], ...
      'FaceColor', 'w', 'FaceAlpha', 0.85, 'EdgeColor', [0.5 0.5 0.5]);
plot(leg_x0+0.1, leg_y0+0.75*leg_h, 'o', 'MarkerSize', 7, 'MarkerFaceColor', gray_color, 'MarkerEdgeColor', 'k');
plot(leg_x0+0.1, leg_y0+0.75*leg_h, 'w.', 'MarkerSize', 3);
text(leg_x0+0.22, leg_y0+0.75*leg_h, '电流流出 (Out)', 'FontSize', 9, 'VerticalAlignment', 'middle');
plot(leg_x0+0.1, leg_y0+0.25*leg_h, 'o', 'MarkerSize', 7, 'MarkerFaceColor', gray_color, 'MarkerEdgeColor', 'k');
text(leg_x0+0.1, leg_y0+0.25*leg_h, '×', 'Color', 'w', 'FontSize', 7, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(leg_x0+0.22, leg_y0+0.25*leg_h, '电流流入 (In)', 'FontSize', 9, 'VerticalAlignment', 'middle');

% 坐标轴与标题设置
title('单圆电流环磁场仿真 (轴截面)'); 
xlabel('z (轴向)'); ylabel('r (径向)');
axis equal; axis([-1.5 1.5 -1.2 1.2]); 
box on;