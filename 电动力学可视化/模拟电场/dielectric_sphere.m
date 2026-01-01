%% 均匀外电场中的介质球模拟 
clc; clear; close all;

% 物理参数设置
a = 0.4;          % 介质球半径
E0 = 2;           % 初始外电场强度
k1 = 0.05;        % 球体介电常数 (或电导率)
k0 = 0.01;        % 外部介质介电常数

% 建立全平面网格
res = 0.01; 
[X, Y] = meshgrid(-1.2:res:1.2);
[theta, r] = cart2pol(X, Y); 

% 计算电势分布（介质球物理公式）
factor = (k1 - k0) / (k1 + 2*k0);  % 极化系数
Z = zeros(size(X));
mask_out = (r >= a);
mask_in  = (r < a);

% 外部电势：均匀场 + 偶极子感应场
Z(mask_out) = -E0 .* r(mask_out) .* cos(theta(mask_out)) + ...
               factor * E0 * a^3 * cos(theta(mask_out)) ./ r(mask_out).^2;
% 内部电势：均匀场（强度被介电常数调制）
Z(mask_in) = -(3*k0 / (k1 + 2*k0)) * E0 .* r(mask_in) .* cos(theta(mask_in));

% 计算电场强度矢量（电势的负梯度）
[Ex, Ey] = gradient(-Z, res);

% 可视化绘图配置
figure('Color', 'w', 'Position', [100, 100, 850, 650]);
hold on; axis equal; grid on;

% 绘制电势等势面填充
[~, hFill] = contourf(X, Y, Z, 60);
set(hFill, 'LineStyle', 'none');
colormap(jet(256));
cb = colorbar;
ylabel(cb, '电势');

% 绘制电场流线（左侧均匀起始）
sy = -1.2:0.08:1.2;
sx = -1.2 * ones(size(sy));
hLines = streamline(X, Y, Ex, Ey, sx, sy);
set(hLines, 'Color', [0.1 0.1 0.1], 'LineWidth', 1);

% 绘制电场矢量箭头（采样降密，避免拥挤）
skip = 7; 
quiver(X(1:skip:end, 1:skip:end), Y(1:skip:end, 1:skip:end), ...
       Ex(1:skip:end, 1:skip:end), Ey(1:skip:end, 1:skip:end), ...
       1.2, 'r', 'LineWidth', 0.7);

% 绘制介质球轮廓（实线边界）
t = linspace(0, 2*pi, 100);
plot(a*cos(t), a*sin(t), 'k-', 'LineWidth', 1.5);

% 坐标轴与标签修饰
xlabel('x'); ylabel('y');
xlim([-1.2 1.2]); ylim([-1.2 1.2]);