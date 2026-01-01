%% 双平行长直导线磁场模拟 
clc; clear; 

% 参数配置
params.I1 = 5;          % 导线1电流
params.I2 = 5;          % 导线2电流
params.d = 0.4;         % 导线距中心距离
params.k = 1;           % 比例系数
params.gray_color = [0.4 0.4 0.4];

% 生成同向/反向电流对比图
draw_magnetic_field(params, 1);
draw_magnetic_field(params, -1);

%% 核心绘图函数
function draw_magnetic_field(p, direction_flag)
    res = 0.01;
    [X, Y] = meshgrid(-1.2:res:1.2, -1:res:1);

    % 左侧导线磁场分量计算
    r1 = sqrt((X+p.d).^2 + Y.^2) + eps;
    theta1 = atan2(Y, X+p.d);
    Bx1 = -(p.k * p.I1 ./ r1) .* sin(theta1);
    By1 = (p.k * p.I1 ./ r1) .* cos(theta1);
    
    % 右侧导线磁场分量计算
    r2 = sqrt((X-p.d).^2 + Y.^2) + eps;
    theta2 = atan2(Y, X-p.d);
    Bx2 = -(p.k * (p.I2 * direction_flag) ./ r2) .* sin(theta2);
    By2 = (p.k * (p.I2 * direction_flag) ./ r2) .* cos(theta2);

    % 总磁场与幅值计算
    Bx_tot = Bx1 + Bx2;
    By_tot = By1 + By2;
    Bmag = sqrt(Bx_tot.^2 + By_tot.^2);
    Bmag(Bmag > 10) = 10;

    % 创建绘图窗口（位置错开）
    if direction_flag == 1
        figure('Color', 'w', 'Units', 'normalized', 'Position', [0.05, 0.2, 0.45, 0.5]);
    else
        figure('Color', 'w', 'Units', 'normalized', 'Position', [0.5, 0.2, 0.45, 0.5]);
    end
    hold on; box on;

    % 绘制磁感应强度填充背景
    contourf(X, Y, Bmag, 100, 'LineStyle', 'none');
    colormap(jet);
    hcb = colorbar; ylabel(hcb, '磁感应强度 |B|');

    % 绘制磁感流线
    edge_x = linspace(-1.15, 1.15, 10);
    edge_y = linspace(-0.95, 0.95, 10);
    [sx1, sy1] = meshgrid(edge_x, [-0.95, 0.95]); 
    [sx2, sy2] = meshgrid([-1.15, 1.15], edge_y); 
    hL = streamline(X, Y, Bx_tot, By_tot, [sx1(:); sx2(:)], [sy1(:); sy2(:)]);
    set(hL, 'Color', [0.1 0.1 0.1], 'LineWidth', 0.6);

    % 绘制归一化磁场矢量箭头
    M = 18;
    u = Bx_tot(1:M:end, 1:M:end) ./ (Bmag(1:M:end, 1:M:end) + eps);
    v = By_tot(1:M:end, 1:M:end) ./ (Bmag(1:M:end, 1:M:end) + eps);
    quiver(X(1:M:end, 1:M:end), Y(1:M:end, 1:M:end), u, v, 0.4, 'r', 'LineWidth', 0.8);

    % 绘制导线实体与电流方向标识
    for pos = [-p.d, p.d]
        is_left = (pos == -p.d);
        plot(pos, 0, 'o', 'MarkerSize', 11, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', p.gray_color, 'LineWidth', 1);
        if is_left || direction_flag == 1
            plot(pos, 0, 'w.', 'MarkerSize', 4);
        else
            text(pos, 0, '×', 'Color', 'w', 'FontSize', 8, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        end
    end

    % 绘制图例框
    leg_w = 0.35; leg_h = 0.12; leg_x0 = 0.8; leg_y0 = 0.85;
    patch('XData', [leg_x0, leg_x0+leg_w, leg_x0+leg_w, leg_x0], ...
          'YData', [leg_y0, leg_y0, leg_y0+leg_h, leg_y0+leg_h], ...
          'FaceColor', 'w', 'FaceAlpha', 0.8, 'EdgeColor', [0.6 0.6 0.6]);
    plot(leg_x0+0.04, leg_y0+0.7*leg_h, 'o', 'MarkerSize', 6, 'MarkerFaceColor', p.gray_color);
    plot(leg_x0+0.04, leg_y0+0.7*leg_h, 'w.', 'MarkerSize', 2);
    text(leg_x0+0.08, leg_y0+0.7*leg_h, ' 电流流出', 'FontSize', 8.5, 'VerticalAlignment', 'middle');
    plot(leg_x0+0.04, leg_y0+0.3*leg_h, 'o', 'MarkerSize', 6, 'MarkerFaceColor', p.gray_color);
    text(leg_x0+0.04, leg_y0+0.3*leg_h, '×', 'Color', 'w', 'FontSize', 6, 'HorizontalAlignment', 'center');
    text(leg_x0+0.08, leg_y0+0.3*leg_h, ' 电流流入', 'FontSize', 8.5, 'VerticalAlignment', 'middle');

    % 标题与坐标轴配置
    if direction_flag == 1
        mode = '同向';
    else
        mode = '反向';
    end
    title(['双平行直导线磁场 (', mode, '电流)'], 'FontSize', 12);
    xlabel('x'); ylabel('y'); axis equal; axis([-1.2 1.2 -1 1]);
    set(gca, 'Layer', 'top', 'TickDir', 'out');
end