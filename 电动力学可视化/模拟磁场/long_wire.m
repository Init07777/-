%% 无限长直导线磁场模拟 
clc; clear; 

% 参数配置
I = 5;            % 导线电流大小
mu0 = 4*pi*1e-7;  % 真空磁导率
k = 2;            % 可视化比例系数
gray_color = [0.4 0.4 0.4];

% 生成计算网格并转换极坐标
res = 0.01;
[X, Y] = meshgrid(-1:res:1, -1:res:1);
[theta, rho] = cart2pol(X, Y);

% 磁感应强度计算与幅值截断
Bmag = k ./ (rho + eps);
Bmag(Bmag > 8) = 8;

% 磁场矢量分量计算
Bx = -Bmag .* sin(theta);   
By = Bmag .* cos(theta);

% 可视化绘图
figure('Color', 'w', 'Name', 'Straight Wire Field', 'Units', 'normalized', 'Position', [0.1, 0.1, 0.5, 0.6]);
clf; hold on; box on;

% 绘制磁场强度填充背景
contourf(X, Y, Bmag, 100, 'LineStyle', 'none');
colormap(jet);
hcb = colorbar;
ylabel(hcb, '磁感应强度 |B|');

% 绘制环形磁感流线
seed_rho = [0.15, 0.3, 0.45, 0.6, 0.75, 0.9];
seed_x = seed_rho;
seed_y = zeros(size(seed_x));
hL = streamline(X, Y, Bx, By, seed_x, seed_y);
set(hL, 'Color', [0.1 0.1 0.1], 'LineWidth', 0.8);

% 绘制归一化红色磁场矢量箭头
M = 15;
u = Bx(1:M:end, 1:M:end) ./ (Bmag(1:M:end, 1:M:end) + eps);
v = By(1:M:end, 1:M:end) ./ (Bmag(1:M:end, 1:M:end) + eps);
quiver(X(1:M:end, 1:M:end), Y(1:M:end, 1:M:end), u, v, 0.4, 'r', 'LineWidth', 0.8);

% 绘制导线横截面标识
plot(0, 0, 'o', 'MarkerSize', 12, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', gray_color, 'LineWidth', 1.2);
plot(0, 0, 'w.', 'MarkerSize', 4);

% 绘制电流方向图例框
leg_w = 0.4; leg_h = 0.12; 
leg_x0 = 0.55; leg_y0 = 0.85;
patch('XData', [leg_x0, leg_x0+leg_w, leg_x0+leg_w, leg_x0], ...
      'YData', [leg_y0, leg_y0, leg_y0+leg_h, leg_y0+leg_h], ...
      'FaceColor', 'w', 'FaceAlpha', 0.85, 'EdgeColor', [0.5 0.5 0.5]);

plot(leg_x0+0.05, leg_y0+0.5*leg_h, 'o', 'MarkerSize', 8, 'MarkerFaceColor', gray_color, 'MarkerEdgeColor', 'k');
plot(leg_x0+0.05, leg_y0+0.5*leg_h, 'w.', 'MarkerSize', 3);
text(leg_x0+0.1, leg_y0+0.5*leg_h, '电流流出 (Out)', 'FontSize', 9, 'VerticalAlignment', 'middle');

% 坐标轴与标题配置
axis equal; axis([-1 1 -1 1]);
title('无限长直导线磁场分布 (横截面)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('x'); ylabel('y');
set(gca, 'Layer', 'top', 'TickDir', 'out');
hold off;