%% 均匀外电场中的导体球模拟
clc; clear; 

% 参数设置
a = 0.4;          % 导体球半径
E0 = 2;           % 初始均匀电场强度 (沿x轴正方向)
phi0 = 0;         % 导体球电势（零势能点）

% 建立空间网格并转换极坐标
[X, Y] = meshgrid(-1:0.02:1);
[theta, r] = cart2pol(X, Y);

% 计算空间电势分布，球内电势置0
Z = -E0 .* r .* cos(theta) .* (1 - a^3 ./ r.^3);
Z(r < a) = 0;

% 计算电场强度矢量（电势的负梯度）
[Ex, Ey] = gradient(-Z, 0.02);

% 消除边界数值误差，强制球内电场为0
Ex(r < a*0.98) = 0; 
Ey(r < a*0.98) = 0;

% 可视化绘图
figure('Color', 'w', 'Name', 'Conducting Sphere in Uniform Electric Field');
hold on;

% 绘制填充等势面
[C, h] = contourf(X, Y, Z, 50);
set(h, 'LineStyle', 'none');
colormap(jet);

% 添加颜色条并标注电势
hcb = colorbar;
ylabel(hcb, '电势 ', 'FontSize', 10); 

% 绘制蓝色等势线
contour(X, Y, Z, 15, 'b', 'LineWidth', 0.5);

% 绘制电场流线
sy = -1:0.1:1;
sx = -1 * ones(size(sy));
line_handle = streamline(X, Y, Ex, Ey, sx, sy);
set(line_handle, 'Color', [0.2 0.2 0.2], 'LineWidth', 1.2);

% 绘制红色电场矢量箭头
M = 5; 
quiver(X(1:M:end, 1:M:end), Y(1:M:end, 1:M:end), ...
       Ex(1:M:end, 1:M:end), Ey(1:M:end, 1:M:end), 0.8, 'r');

% 绘制导体球（灰色填充+黑色轮廓）
theta_draw = linspace(0, 2*pi, 100);
fill(a*cos(theta_draw), a*sin(theta_draw), [0.8 0.8 0.8]);
plot(a*cos(theta_draw), a*sin(theta_draw), 'k-', 'LineWidth', 2);

% 坐标轴设置
axis equal;
axis([-1 1 -1 1]);
xlabel('x'); ylabel('y');