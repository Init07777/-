clc; clear; close all;

% 1. 物理参数设置（交变场核心参数）
E0 = 2;          % 电场振幅（对应原代码E0）
B0 = E0/3e8;     % 磁场振幅（电磁波E/B=光速c）
k = 2*pi/0.5;    % 波数（波长0.5）
w = 3e8*k;       % 角频率（光速c=w/k）
t = 0;           % 初始时间
dt = 0.001;      % 时间步长

% 2. 网格生成（复刻原代码）
[X,Y] = meshgrid(-1:0.01:1);

% 3. 动态可视化循环（核心：时变场刷新）
figure('Position',[100,100,800,600]);
while t < 2  % 模拟2秒时长
    % 交变电场+磁场计算（沿X方向传播，E沿Y、B沿Z）
    Ey = E0 * sin(k*X - w*t);  % 电场Y分量
    Bz = B0 * sin(k*X - w*t);  % 磁场Z分量
    
    % 绘图刷新（完全复刻原代码框架）
    cla; hold on; box on;
    % 电场分布填色等高线
    [C,h] = contourf(X,Y,Ey,100); set(h,'LineStyle','none');
    % 电场矢量箭头（红色）
    M = 20;
    [~,dEy] = gradient(Ey,0.01);
    quiver(X(M:M:end-M,M:M:end-M),Y(M:M:end-M,M:M:end-M),...
           zeros(size(X(M:M:end-M,M:M:end-M))),dEy(M:M:end-M,M:M:end-M),1,'r');
    % 等场强线标注
    contour(X,Y,Ey,8,'b','ShowText','on');
    % 电磁波传播流线（对应原代码streamline）
    [VX,VY] = gradient(Ey,0.01);
    streamline(X,Y,VX,VY,0:0.2:1,zeros(1,6));
    
    % 图形美化
    axis equal;
    colorbar;
    title(sprintf('均匀交变电磁场分布（t=%.3fs，E⊥B，沿X传播）',t));
    xlabel('X轴'); ylabel('Y轴');
    drawnow;  % 动态刷新
    t = t + dt;
end
hold off;