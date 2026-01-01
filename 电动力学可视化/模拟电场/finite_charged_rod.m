%% 有限长带电细棒电场模拟 
clc; clear;

% 参数设置
L = 0.6;           % 半棒长 (总长 1.2)
Q = 1;             % 总电量 (正电荷)
res = 0.01;        % 空间分辨率
rod_w = 0.02;      % 棒可视化宽度
V_clamping = 4.0;  % 电势显示截断值 

% 建立空间网格
[X, Y] = meshgrid(-1.3:res:1.3, -1.3:res:1.3);

% 计算电势分布（对称形式）
r1 = hypot(X, Y - L); % 网格点到棒顶端距离
r2 = hypot(X, Y + L); % 网格点到棒底端距离
V_raw = (Q / (2*L)) * log((r1 + r2 + 2*L) ./ (r1 + r2 - 2*L + eps));

% 数据预处理
[Ex, Ey] = gradient(-V_raw, res);  % 计算电场强度
V_display = V_raw;
V_display(V_display > V_clamping) = V_clamping;  % 电势截断处理

% 可视化绘图
figure('Color', 'w', 'Name', 'Symmetric Electric Field', 'Position', [100 100 750 650]);
hold on; axis equal;

% 绘制电势填充背景
[~, hF] = contourf(X, Y, V_display, 100, 'LineStyle', 'none');
colormap(jet(256)); 
caxis([0, V_clamping]);
hcb = colorbar;
ylabel(hcb, '电势');

% 绘制蓝色等势线
contour(X, Y, V_display, linspace(0.2, V_clamping-0.5, 15), 'b', 'LineWidth', 0.5);

% 绘制电场流线
n_side = 14;
sy_side = linspace(-L*0.95, L*0.95, n_side);
sx_start = [ones(1,n_side)*rod_w*1.5, -ones(1,n_side)*rod_w*1.5];
sy_start = [sy_side, sy_side];
t = linspace(0.05*pi, 0.95*pi, 10); 
sx_end = rod_w * 3 * cos(t);
sy_end_top = L + rod_w * 1.5 * sin(t);
sy_end_bot = -L - rod_w * 1.5 * sin(t);

hL = streamline(X, Y, Ex, Ey, [sx_start, sx_end, sx_end], [sy_start, sy_end_top, sy_end_bot]);
set(hL, 'Color', [0.2 0.2 0.2], 'LineWidth', 1);

% 绘制红色单位化电场矢量箭头
M = 14;
E_mag = hypot(Ex, Ey) + eps;
u_plot = Ex ./ E_mag;
v_plot = Ey ./ E_mag;

in_rod = (abs(X) < rod_w * 1.5) & (abs(Y) < L * 1.02);
u_plot(in_rod) = NaN; 
v_plot(in_rod) = NaN;

quiver(X(1:M:end, 1:M:end), Y(1:M:end, 1:M:end), ...
       u_plot(1:M:end, 1:M:end), v_plot(1:M:end, 1:M:end), ...
       0.5, 'r', 'LineWidth', 1);

% 绘制带电细棒（灰色填充+黑色轮廓）
fill([-rod_w rod_w rod_w -rod_w], [-L -L L L], [0.4 0.4 0.4], 'EdgeColor', 'k', 'LineWidth', 1.5);

% 坐标轴与界面修饰
xlabel('x'); ylabel('y');
axis([-1.3 1.3 -1.3 1.3]);
box on;
set(gca, 'Layer', 'top', 'TickDir', 'out', 'FontSize', 10);