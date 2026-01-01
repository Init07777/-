%% 5. 矩形波导 TE10 模式 (Top-Down View, H-field loops)
clc; clear; close all;

% --- 参数 ---
a = 1;              % 波导宽
lambda_g = 1.5;     % 波导波长
L_view = 2.5;       % 观察长度
res = 0.02;

% --- 网格 (x-z 平面, y=b/2 处) ---
[z, x] = meshgrid(0:res:L_view, 0:res:a);

% --- 物理计算 (TE10 Mode) ---
beta = 2*pi / lambda_g;
% Ey (电场，标量，垂直纸面)
Ey = sin(pi*x/a) .* sin(beta*z); 

% Hx, Hz (磁场，矢量，在纸面内)
Hx = -sin(pi*x/a) .* sin(beta*z); % 简化系数，仅关注形态
Hz = (lambda_g/(2*a)) * cos(pi*x/a) .* cos(beta*z);

% --- 绘图 ---
figure('Color', 'w', 'Name', 'Waveguide TE10'); hold on; axis equal;

% A. 背景：电场强度 Ey (垂直分量)
[~, hF] = contourf(z, x, Ey, 50, 'LineStyle', 'none');
colormap(jet); 
cb = colorbar; ylabel(cb, '电场 Ey (垂直纸面)');

% B. 磁场线 (流线 H) - 形成闭合环
sz = linspace(0, L_view, 20); sx = 0.5 * a * ones(size(sz));
hL = streamline(z, x, Hz, Hx, sz, sx); % 注意: plot(z,x) 所以 u=Hz, v=Hx
set(hL, 'Color', [0.2 0.2 0.2], 'LineWidth', 1.2);

% C. 磁场矢量 (箭头 H)
M = 4;
quiver(z(1:M:end, 1:M:end), x(1:M:end, 1:M:end), ...
       Hz(1:M:end, 1:M:end), Hx(1:M:end, 1:M:end), 1.0, 'r');

% D. 波导边界
plot([0 L_view], [0 0], 'k-', 'LineWidth', 3);
plot([0 L_view], [a a], 'k-', 'LineWidth', 3);

title('矩形波导 TE10 模式 (纵截面: 颜色=Ey, 流线=H场)');
xlabel('z (传播方向)'); ylabel('x (宽边)');
axis([0 L_view 0 a]);