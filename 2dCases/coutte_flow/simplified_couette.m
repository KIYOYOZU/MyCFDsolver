%% simplified_couette.m
% Couette流启动过程数值模拟 - 简化教学版
%
% 功能：
% 1. 显式欧拉格式求解 ∂u/∂t = ν·∂²u/∂y²
% 2. 计算解析解（Fourier级数）
% 3. 生成可视化图形（3个静态图 + 1个动画）
%
% 作者：[自动生成]
% 日期：2025-12-04
%
% 输出文件：
%   - results/couette_results.mat (数据文件)
%   - results/velocity_evolution.png (图1：速度剖面演化)
%   - results/midplane_comparison.png (图2：8时刻对比)
%   - results/steady_state_2D.png (图3：稳态云图)
%   - results/couette_animation.mp4 (图4：动画)
%
% ⚠️ 注意：
% 本文件保留作为“单文件教学版”备份。
% 模块化版本请使用：
%   1) run_couette.m            % 计算脚本
%   2) postprocess_couette.m   % 后处理脚本

clear; clc; close all;

%% ========== Part 1: 参数设置 ==========
% 物理参数
Nx = 32;                    % X方向网格点数
Ny = 51;                    % Y方向网格点数
Lx = 2*pi;                  % X方向域长度
Ly = 2.0;                   % Y方向域高度（通道高度h）
U = 1.0;                    % 上壁速度
nu = 0.1;                   % 运动粘度

% 时间参数
dt = 0.0001;                % 时间步长
t_end = 20.0;               % 终止时间
n_snapshots = 51;           % 快照数量（包括t=0）

% 计算派生参数
dx = Lx / (Nx - 1);
dy = Ly / (Ny - 1);
alpha_y = nu * dt / (dy^2);

% 稳定性检查
fprintf('========== Couette流求解器 ==========\n');
fprintf('网格: %d × %d\n', Nx, Ny);
fprintf('时间步长: dt = %.6f s\n', dt);
fprintf('空间步长: dx = %.6f, dy = %.6f\n', dx, dy);
fprintf('扩散数 α_y = %.4f\n', alpha_y);

if alpha_y > 0.5
    error('稳定性条件违反：α_y = %.4f > 0.5', alpha_y);
end
fprintf('稳定性检查通过 (α_y ≤ 0.5)\n');
fprintf('====================================\n\n');

%% ========== Part 2: 网格初始化 ==========
x = linspace(0, Lx, Nx);
y = linspace(0, Ly, Ny);
u = zeros(Nx, Ny);          % 初始条件：u(x,y,0) = 0
u_new = zeros(Nx, Ny);

% 应用初始边界条件
u = apply_BC(u, Nx, Ny, U);

%% ========== Part 3: 预分配存储空间 ==========
output_interval = t_end / (n_snapshots - 1);
u_history = cell(1, n_snapshots);
time_history = zeros(1, n_snapshots);
L2_history = zeros(1, n_snapshots);

% 存储初始状态
snapshot_idx = 1;
u_history{snapshot_idx} = u;
time_history(snapshot_idx) = 0.0;

fprintf('存储配置:\n');
fprintf('  - 快照数量: %d\n', n_snapshots);
fprintf('  - 快照间隔: %.4f s\n', output_interval);
fprintf('  - 预计总时间步数: %d\n\n', round(t_end/dt));

%% ========== Part 4: 主时间推进循环 ==========
fprintf('开始时间推进...\n');
fprintf('------------------------------------\n');

current_time = 0.0;
step = 0;
next_output_time = output_interval;

tic;  % 计时开始

while current_time < t_end
    step = step + 1;

    % 显式欧拉推进（向量化实现）
    u_new = time_step(u, Nx, Ny, alpha_y);

    % 应用边界条件
    u = apply_BC(u_new, Nx, Ny, U);

    % 更新时间
    current_time = current_time + dt;

    % 保存快照
    if current_time >= next_output_time - 1e-10 && snapshot_idx < n_snapshots
        snapshot_idx = snapshot_idx + 1;
        u_history{snapshot_idx} = u;
        time_history(snapshot_idx) = current_time;
        next_output_time = snapshot_idx * output_interval;

        % 快照保存提示
        if mod(snapshot_idx, 10) == 0
            fprintf('快照 %d/%d: t=%.3fs\n', snapshot_idx, n_snapshots, current_time);
        end
    end

    % 进度显示
    if mod(step, 10000) == 0
        fprintf('进度: %.1f%% (t=%.3fs, 步数=%d)\n', ...
                current_time/t_end*100, current_time, step);
    end
end

% 保存最终状态（如果未保存）
if snapshot_idx < n_snapshots
    snapshot_idx = snapshot_idx + 1;
    u_history{snapshot_idx} = u;
    time_history(snapshot_idx) = current_time;
end

elapsed_time = toc;
fprintf('------------------------------------\n');
fprintf('时间推进完成！耗时: %.2f 秒\n', elapsed_time);
fprintf('实际计算步数: %d\n', step);
fprintf('平均速度: %.0f 步/秒\n\n', step/elapsed_time);

%% ========== Part 5: 计算解析解和误差 ==========
fprintf('计算解析解和误差...\n');

u_analytical = cell(1, n_snapshots);

for i = 1:n_snapshots
    u_analytical{i} = compute_analytical(y, time_history(i), U, nu, Ly);

    % 计算L2误差
    u_ana_2D = repmat(u_analytical{i}', Nx, 1);
    L2_history(i) = compute_L2_error(u_history{i}, u_ana_2D);

    % 显示进度
    if mod(i, 10) == 0 || i == n_snapshots
        fprintf('  解析解计算进度: %d/%d (L2误差=%.6e)\n', i, n_snapshots, L2_history(i));
    end
end

fprintf('\n误差分析:\n');
fprintf('  - L2误差范围: [%.4e, %.4e]\n', min(L2_history), max(L2_history));
fprintf('  - 最终L2误差: %.4e\n', L2_history(end));
fprintf('  - 相对精度: %.6f%%\n\n', L2_history(end)*100);

%% ========== Part 6: 创建输出目录 ==========
output_dir = 'results';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    fprintf('创建输出目录: %s\n\n', output_dir);
end

%% ========== Part 7: 可视化 ==========
fprintf('========== 生成可视化图形 ==========\n');

% 图1：速度剖面演化对比
fprintf('[1/4] 生成速度剖面演化图...\n');
plot_velocity_evolution(u_history, u_analytical, x, y, time_history, output_dir);

% 图2：8时刻对比图
fprintf('[2/4] 生成8时刻对比图...\n');
plot_midplane_comparison(u_history, u_analytical, x, y, time_history, output_dir);

% 图3：稳态2D流场
fprintf('[3/4] 生成稳态2D流场图...\n');
plot_steady_state_2D(u_history, x, y, time_history, output_dir);

% 图4：动画生成
fprintf('[4/4] 生成演化动画...\n');
generate_animation(u_history, u_analytical, x, y, time_history, output_dir);

fprintf('====================================\n\n');

%% ========== Part 8: 保存数据 ==========
fprintf('保存数据到MAT文件...\n');
save(fullfile(output_dir, 'couette_results.mat'), ...
     'x', 'y', 'u_history', 'time_history', 'L2_history', 'u_analytical', ...
     'Nx', 'Ny', 'Lx', 'Ly', 'U', 'nu', 'dt', 'alpha_y', 't_end', 'n_snapshots');
fprintf('  ✓ 数据已保存\n\n');

%% ========== Part 9: 最终总结 ==========
fprintf('\n========== 求解完成 ==========\n');
fprintf('总时间步数: %d\n', step);
fprintf('最终时间: %.4f s\n', current_time);
fprintf('最终L2误差: %.6e\n', L2_history(end));
fprintf('\n输出文件:\n');
fprintf('  - %s/couette_results.mat\n', output_dir);
fprintf('  - %s/velocity_evolution.png\n', output_dir);
fprintf('  - %s/midplane_comparison.png\n', output_dir);
fprintf('  - %s/steady_state_2D.png\n', output_dir);
fprintf('  - %s/couette_animation.mp4\n', output_dir);
fprintf('================================\n');

%% ========== 子函数定义 ==========

% 时间推进函数（显式欧拉格式，向量化实现）
function u_new = time_step(u, Nx, Ny, alpha_y)
    % 功能：使用显式欧拉格式推进一个时间步
    % 离散格式：u^(n+1)_{i,j} = u^n_{i,j} + α_y·(u^n_{i,j+1} - 2·u^n_{i,j} + u^n_{i,j-1})
    %
    % 输入：
    %   u      - 当前时刻速度场 [Nx×Ny]
    %   Nx, Ny - 网格点数
    %   alpha_y - 扩散数（Y方向）
    % 输出：
    %   u_new  - 下一时刻速度场 [Nx×Ny]

    u_new = u;

    % 仅更新内部节点（边界节点由边界条件函数处理）
    u_new(2:Nx-1, 2:Ny-1) = u(2:Nx-1, 2:Ny-1) + ...
        alpha_y * (u(2:Nx-1, 3:Ny) - 2*u(2:Nx-1, 2:Ny-1) + u(2:Nx-1, 1:Ny-2));
end

% 边界条件函数
function u = apply_BC(u, Nx, Ny, U)
    % 功能：应用边界条件
    %
    % 边界条件设置：
    %   Y方向：Dirichlet边界条件
    %     - 下壁（y=0）：u=0（静止壁面）
    %     - 上壁（y=h）：u=U（运动壁面）
    %   X方向：周期性边界条件
    %     - 左边界：u(0,y) = u(Lx-dx,y)
    %     - 右边界：u(Lx,y) = u(dx,y)
    %
    % 输入：
    %   u   - 速度场 [Nx×Ny]
    %   Nx, Ny - 网格点数
    %   U   - 上壁速度
    % 输出：
    %   u   - 应用边界条件后的速度场

    % Y方向：Dirichlet边界条件
    u(:, 1) = 0;        % 下壁静止
    u(:, Ny) = U;       % 上壁运动

    % X方向：周期性边界条件
    u(1, :) = u(Nx-1, :);    % 左边界 = 右侧倒数第二列
    u(Nx, :) = u(2, :);      % 右边界 = 左侧第二列
end

% 解析解计算函数（Fourier级数展开）
function u_ana = compute_analytical(y, t, U, nu, h)
    % 功能：计算Couette流启动问题的解析解
    %
    % 理论解析解（上壁运动）：
    %   u(y,t) = U·(y/h) - (2U/π)·Σ[(-1)^(n+1)/n·sin(nπy/h)·exp(-n²π²νt/h²)]
    %
    % 物理意义：
    %   - 稳态分量：U·(y/h) 线性分布
    %   - 瞬态分量：Fourier级数项随时间指数衰减
    %
    % 输入：
    %   y  - Y方向坐标向量 [Ny×1]
    %   t  - 当前时间
    %   U  - 上壁速度
    %   nu - 运动粘度
    %   h  - 通道高度
    % 输出：
    %   u_ana - 解析解速度剖面 [Ny×1]

    % 确保y是列向量
    y_col = y(:);                                % [Ny×1]

    % 稳态分量（线性分布）
    u_steady = U * (y_col / h);                  % [Ny×1]

    % 瞬态分量（Fourier级数100项）
    n = (1:100)';                                % 模态数 [100×1]
    coeff = (-1).^(n+1) .* (2*U) ./ (n*pi);     % 系数 [100×1]

    spatial = sin(n * pi * y_col' / h);          % 空间项 [100×Ny]
    temporal = exp(-n.^2 * pi^2 * nu * t / h^2); % 时间项 [100×1]

    % 瞬态分量求和
    u_transient = -sum(coeff .* temporal .* spatial, 1)';  % [Ny×1]

    % 总解析解
    u_ana = u_steady + u_transient;              % [Ny×1]
end

% L2误差计算函数
function L2_err = compute_L2_error(u_num, u_ana_2D)
    % 功能：计算相对L2范数误差
    %
    % L2误差定义：
    %   L2_err = ||u_numerical - u_analytical||_2 / ||u_analytical||_2
    %
    % 输入：
    %   u_num    - 数值解 [Nx×Ny]
    %   u_ana_2D - 解析解（2D扩展）[Nx×Ny]
    % 输出：
    %   L2_err   - 相对L2误差（标量）

    error_sq = sum((u_num - u_ana_2D).^2, 'all');  % 误差平方和
    ana_sq = sum(u_ana_2D.^2, 'all');              % 解析解平方和
    L2_err = sqrt(error_sq / ana_sq);              % 相对L2误差
end

% 可视化函数1：速度剖面演化对比
function plot_velocity_evolution(u_history, u_analytical, x, y, time_history, output_dir)
    % 功能：绘制速度剖面演化图（选择5个关键时刻）
    %
    % 输出内容：
    %   - 5条彩色曲线：t=0%, 25%, 50%, 75%, 100% 的数值解
    %   - 1条黑色虚线：稳态解析解参考线
    %
    % 输入：
    %   u_history    - 速度场历史记录（cell数组）
    %   u_analytical - 解析解历史记录（cell数组）
    %   x, y         - 网格坐标
    %   time_history - 时间历史记录
    %   output_dir   - 输出目录路径

    % 创建图形窗口
    fig = figure('Name', 'Velocity Profile Evolution', ...
                 'Position', [100, 100, 1000, 700], 'Color', 'w');

    % 获取中心线索引（x=π）
    midplane_idx = round(length(x) / 2);

    % 选择关键时刻索引（0%, 25%, 50%, 75%, 100%）
    num_frames = length(time_history);
    frame_indices = unique(round([1, 0.25*num_frames, 0.5*num_frames, ...
                                  0.75*num_frames, num_frames]));

    % 定义颜色方案
    colors = lines(length(frame_indices));

    hold on;

    % 首先绘制稳态解析解（粗黑色虚线）
    u_steady_analytical = u_analytical{end};
    plot(u_steady_analytical, y, 'k--', 'LineWidth', 2.5, ...
         'DisplayName', 'Steady State (analytical)');

    % 绘制各时刻的数值解
    for j = 1:length(frame_indices)
        idx = frame_indices(j);
        t = time_history(idx);

        % 提取中心线速度剖面 u(x=π, y)
        u_profile = u_history{idx}(midplane_idx, :);

        % 绘制数值解曲线
        plot(u_profile, y, '-', 'Color', colors(j,:), 'LineWidth', 1.8, ...
             'DisplayName', sprintf('t = %.2f (numerical)', t));
    end

    % 图形美化
    xlabel('u/U', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('y/h', 'FontSize', 13, 'FontWeight', 'bold');
    title('Velocity Profile Evolution (Numerical vs Analytical)', ...
          'FontSize', 15, 'FontWeight', 'bold');

    % 图例
    legend('Location', 'best', 'FontSize', 11);

    % 网格和坐标轴
    grid on;
    box on;
    xlim([0, 1.05]);
    ylim([0, 2.0]);
    set(gca, 'FontSize', 11);

    % 保存图形
    output_path = fullfile(output_dir, 'velocity_evolution.png');
    print(output_path, '-dpng', '-r300');
    fprintf('   ✓ 已保存: %s\n', output_path);

    close(fig);
end

% 可视化函数2：8时刻对比图（2×4子图）
function plot_midplane_comparison(u_history, u_analytical, x, y, time_history, output_dir)
    % 功能：生成8个子图对比数值解与解析解
    %
    % 输出内容：
    %   - 2×4子图布局，每个子图包含：
    %     * 蓝色实线+圆点：数值解
    %     * 红色虚线：解析解
    %
    % 输入：
    %   u_history    - 速度场历史记录（cell数组）
    %   u_analytical - 解析解历史记录（cell数组）
    %   x, y         - 网格坐标
    %   time_history - 时间历史记录
    %   output_dir   - 输出目录路径

    % 选择8个均匀分布的时刻
    n_snapshots = length(time_history);
    if n_snapshots >= 8
        time_indices = round(linspace(1, n_snapshots, 8));
    else
        time_indices = 1:n_snapshots;
        if length(time_indices) < 8
            time_indices = [time_indices, repmat(n_snapshots, 1, 8 - length(time_indices))];
        end
    end

    time_indices = unique(time_indices);
    n_plots = min(8, length(time_indices));

    % 获取中心线索引（x=π）
    midplane_idx = round(length(x) / 2);

    % 创建图形窗口
    fig = figure('Position', [100, 100, 1400, 900], 'Color', 'white');

    for i = 1:n_plots
        idx = time_indices(i);
        subplot(2, 4, i);

        % 提取中心线速度剖面 u(x=π, y)
        u_midplane = u_history{idx}(midplane_idx, :);

        % 绘制数值解（蓝色实线+圆点）
        plot(u_midplane, y, 'b-o', 'LineWidth', 1.5, ...
             'MarkerSize', 4, 'MarkerFaceColor', 'b', 'DisplayName', 'Numerical');
        hold on;

        % 绘制解析解（红色虚线）
        plot(u_analytical{idx}, y, 'r--', 'LineWidth', 2, ...
             'DisplayName', 'Analytical');

        % 图形美化
        xlabel('u/U', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('y/h', 'FontSize', 11, 'FontWeight', 'bold');
        title(sprintf('t = %.2f', time_history(idx)), 'FontSize', 12, 'FontWeight', 'bold');

        % 图例
        legend('Location', 'best', 'FontSize', 9);

        % 网格和坐标轴
        grid on;
        box on;
        xlim([0, 1]);
        ylim([0, 2]);
        set(gca, 'FontSize', 10);
    end

    % 总标题
    sgtitle('Midplane Velocity Profile vs Analytical Solution (x = \pi)', ...
            'FontSize', 16, 'FontWeight', 'bold');

    % 保存图形
    output_path = fullfile(output_dir, 'midplane_comparison.png');
    print(output_path, '-dpng', '-r300');
    fprintf('   ✓ 已保存: %s\n', output_path);

    close(fig);
end

% 可视化函数3：稳态2D流场
function plot_steady_state_2D(u_history, x, y, time_history, output_dir)
    % 功能：绘制稳态2D速度场云图
    %
    % 输出内容：
    %   - pcolor伪彩色图（shading interp）
    %   - jet色标
    %   - 中心线标记（x=π）
    %
    % 输入：
    %   u_history    - 速度场历史记录（cell数组）
    %   x, y         - 网格坐标
    %   time_history - 时间历史记录
    %   output_dir   - 输出目录路径

    % 获取最终时刻索引
    idx_final = length(time_history);

    % 创建图形窗口
    fig = figure('Position', [100, 100, 1000, 700], 'Color', 'white');

    % 创建网格（归一化坐标）
    [X, Y] = meshgrid(x / (2*pi), y / 2);

    % 绘制2D云图（转置以匹配维度）
    pcolor(X, Y, u_history{idx_final}');
    shading interp;

    % 颜色条
    h_colorbar = colorbar;
    ylabel(h_colorbar, 'u/U', 'FontSize', 13, 'FontWeight', 'bold');
    caxis([0, 1]);
    colormap(jet);

    % 叠加中心线参考线
    hold on;
    midplane_x = 0.5;  % x=π 对应 x/(2π) = 0.5
    plot([midplane_x, midplane_x], [0, 1], 'k--', 'LineWidth', 2.5);
    text(midplane_x + 0.05, 0.93, 'x = \pi', 'FontSize', 12, ...
         'Color', 'k', 'FontWeight', 'bold', 'BackgroundColor', 'white');

    % 图形美化
    xlabel('x/(2\pi)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('y/h', 'FontSize', 14, 'FontWeight', 'bold');
    title(sprintf('Steady State Velocity Field (t = %.2f)', time_history(idx_final)), ...
          'FontSize', 16, 'FontWeight', 'bold');

    % 坐标轴设置
    axis equal tight;
    set(gca, 'FontSize', 12);

    % 保存图形
    output_path = fullfile(output_dir, 'steady_state_2D.png');
    print(output_path, '-dpng', '-r300');
    fprintf('   ✓ 已保存: %s\n', output_path);

    close(fig);
end

% 可视化函数4：动画生成（双面板布局）
function generate_animation(u_history, u_analytical, x, y, time_history, output_dir)
    % 功能：生成演化动画
    %
    % 输出内容：
    %   - 左面板：2D速度场云图
    %   - 右面板：中心剖面速度曲线（数值vs解析）
    %   - 20 fps，MPEG-4格式
    %
    % 输入：
    %   u_history    - 速度场历史记录（cell数组）
    %   u_analytical - 解析解历史记录（cell数组）
    %   x, y         - 网格坐标
    %   time_history - 时间历史记录
    %   output_dir   - 输出目录路径

    fprintf('   [动画] 准备数据...\n');

    %% 1. 数据准备
    n_frames = length(time_history);
    [X, Y] = meshgrid(x, y);  % 预计算网格坐标

    % 中心线索引（x=π）
    x_mid_idx = round(length(x) / 2);

    % 预计算全局速度范围（固定色标）
    u_max = max(cellfun(@(u) max(u(:)), u_history));
    u_min = min(cellfun(@(u) min(u(:)), u_history));

    fprintf('   [动画] 数据准备完成: %d 帧, 时间范围 [%.2f, %.2f]s\n', ...
            n_frames, time_history(1), time_history(end));

    %% 2. 创建图形窗口和布局
    fig = figure('Name', 'Couette Flow Animation', ...
                 'Position', [100, 100, 1400, 600], 'Color', 'w');

    % 左面板：2D速度场云图
    ax_main = subplot(1, 2, 1);
    h_pcolor = pcolor(ax_main, X, Y, u_history{1}');
    shading(ax_main, 'interp');
    colormap(ax_main, 'jet');
    h_colorbar = colorbar(ax_main);
    h_colorbar.Label.String = 'u/U';
    h_colorbar.Label.FontSize = 11;
    caxis(ax_main, [u_min, u_max]);
    xlabel(ax_main, 'x/(2\pi)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel(ax_main, 'y/h', 'FontSize', 12, 'FontWeight', 'bold');
    title(ax_main, 'Velocity Field', 'FontSize', 13, 'FontWeight', 'bold');
    axis(ax_main, 'equal', 'tight');

    % 在云图上绘制剖面位置指示线（红色虚线）
    hold(ax_main, 'on');
    h_profile_line = plot(ax_main, [x(x_mid_idx), x(x_mid_idx)], [y(1), y(end)], ...
                          'r--', 'LineWidth', 2);

    % 右面板：速度剖面曲线
    ax_profile = subplot(1, 2, 2);

    % 提取初始速度剖面
    u_profile_numerical = u_history{1}(x_mid_idx, :)';
    u_profile_analytical = u_analytical{1};

    % 绘制数值解（蓝色实线）
    h_numerical = plot(ax_profile, u_profile_numerical, y, ...
                       'b-', 'LineWidth', 2.0, 'DisplayName', 'Numerical');
    hold(ax_profile, 'on');

    % 绘制解析解（红色虚线）
    h_analytical = plot(ax_profile, u_profile_analytical, y, ...
                        'r--', 'LineWidth', 2.0, 'DisplayName', 'Analytical');

    xlabel(ax_profile, 'u/U', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel(ax_profile, 'y/h', 'FontSize', 12, 'FontWeight', 'bold');
    title(ax_profile, sprintf('Velocity Profile at x = %.2f', x(x_mid_idx)), ...
          'FontSize', 13, 'FontWeight', 'bold');
    grid(ax_profile, 'on');
    legend(ax_profile, 'Location', 'southeast', 'FontSize', 10);
    xlim(ax_profile, [u_min - 0.05, u_max + 0.05]);
    ylim(ax_profile, [y(1), y(end)]);
    set(ax_profile, 'YDir', 'normal');

    fprintf('   [动画] 布局创建完成\n');

    %% 3. 创建视频编码器
    video_path = fullfile(output_dir, 'couette_animation.mp4');
    v = VideoWriter(video_path, 'MPEG-4');
    v.FrameRate = 20;  % 20 帧/秒
    v.Quality = 90;    % 高质量
    open(v);

    fprintf('   [动画] 视频编码器创建完成: %s\n', video_path);
    fprintf('   [动画] 帧率: %d fps, 预计时长: %.1f 秒\n', ...
            v.FrameRate, n_frames / v.FrameRate);

    %% 4. 动画循环
    fprintf('   [动画] 开始渲染帧...\n');

    for i = 1:n_frames
        % 更新左面板：2D速度场
        set(h_pcolor, 'CData', u_history{i}');

        % 更新右面板：速度剖面
        u_profile_numerical = u_history{i}(x_mid_idx, :)';
        u_profile_analytical = u_analytical{i};

        set(h_numerical, 'XData', u_profile_numerical);
        set(h_analytical, 'XData', u_profile_analytical);

        % 更新全局标题
        sgtitle(sprintf('Couette Flow Evolution | t = %.4f s | Frame %d/%d', ...
                        time_history(i), i, n_frames), ...
                'FontSize', 15, 'FontWeight', 'bold');

        % 捕获并写入视频帧
        drawnow;
        frame = getframe(fig);
        writeVideo(v, frame);

        % 进度指示
        if mod(i, 50) == 0 || i == n_frames
            fprintf('      进度: %d/%d (%.1f%%)\n', i, n_frames, 100*i/n_frames);
        end
    end

    %% 5. 清理和保存
    close(v);
    close(fig);

    fprintf('   [动画] 渲染完成！\n');
    fprintf('   ✓ 已保存: %s\n', video_path);

    % 显示文件大小
    if exist(video_path, 'file')
        file_info = dir(video_path);
        fprintf('   [动画] 文件大小: %.2f MB\n', file_info.bytes / (1024*1024));
    end
end
