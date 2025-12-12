%% postprocess_couette.m
% Couette流启动过程后处理脚本 - 模块化版本
%
% 说明：
%   - 读取 run_couette.m 生成的 results/couette_results.mat
%   - 生成 3 张静态图 + 1 个动画

clear; clc; close all;

case_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(case_dir, '..', '..', 'src')));  % 通用模块路径
output_dir = fullfile(case_dir, 'results');
data_path = fullfile(output_dir, 'couette_results.mat');

if ~exist(data_path, 'file')
    error('未找到结果文件 %s，请先运行 run_couette.m', data_path);
end

load(data_path);

fprintf('开始后处理...\n');
fprintf('输出目录: %s\n', output_dir);

%% ========== 解析解与误差（后处理阶段计算） ==========
% 兼容旧数据：可能没有 source_term 字段
G_val = 0.0;
if exist('source_term', 'var')
    G_val = source_term;
elseif exist('G', 'var')
    G_val = G;
end

if ~exist('Nx', 'var'); Nx = size(u_history{1}, 1); end
if ~exist('Ny', 'var'); Ny = size(u_history{1}, 2); end
n_snapshots = length(time_history);

% 传入 (G_val, n_terms) 两个参数，避免 G_val 为整数时被误判成 n_terms
analytical_solver = CouetteAnalytical(U, nu, Ly, G_val, 100);
u_analytical = cell(1, n_snapshots);
L2_history = zeros(1, n_snapshots);

for i = 1:n_snapshots
    u_analytical{i} = analytical_solver.compute_profile(y, time_history(i));
    u_ana_2D = repmat(u_analytical{i}', Nx, 1);
    L2_history(i) = ErrorAnalyzer.compute_L2(u_history{i}, u_ana_2D);
end

fprintf('\n误差分析:\n');
fprintf('  - L2误差范围: [%.4e, %.4e]\n', min(L2_history), max(L2_history));
fprintf('  - 最终L2误差: %.4e\n\n', L2_history(end));

% 三张静态图
plot_velocity_evolution(u_history, u_analytical, x, y, time_history, output_dir);
plot_midplane_comparison(u_history, u_analytical, x, y, time_history, output_dir);
plot_steady_state_2D(u_history, x, y, time_history, output_dir);

% 动画
generate_animation(u_history, u_analytical, x, y, time_history, output_dir);

fprintf('后处理完成！\n');

%% ========== 子函数（从 simplified_couette.m 迁移） ==========

% 可视化函数1：速度剖面演化对比
function plot_velocity_evolution(u_history, u_analytical, x, y, time_history, output_dir)
    % 功能：绘制速度剖面演化图（选择5个关键时刻）

    fig = figure('Name', 'Velocity Profile Evolution', ...
                 'Position', [100, 100, 1000, 700], 'Color', 'w');

    midplane_idx = round(length(x) / 2);

    % 全局速度范围（适配含压力驱动情况）
    u_min = min(cellfun(@(uu) min(uu(:)), u_history));
    u_max = max(cellfun(@(uu) max(uu(:)), u_history));
    margin = 0.05 * (u_max - u_min + eps);

    num_frames = length(time_history);
    frame_indices = unique(round([1, 0.25*num_frames, 0.5*num_frames, ...
                                  0.75*num_frames, num_frames]));

    colors = lines(length(frame_indices));

    hold on;

    u_steady_analytical = u_analytical{end};
    plot(u_steady_analytical, y, 'k--', 'LineWidth', 2.5, ...
         'DisplayName', 'Steady State (analytical)');

    for j = 1:length(frame_indices)
        idx = frame_indices(j);
        t = time_history(idx);
        u_profile = u_history{idx}(midplane_idx, :);
        plot(u_profile, y, '-', 'Color', colors(j,:), 'LineWidth', 1.8, ...
             'DisplayName', sprintf('t = %.2f (numerical)', t));
    end

    xlabel('u/U', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('y/h', 'FontSize', 13, 'FontWeight', 'bold');
    title('Velocity Profile Evolution (Numerical vs Analytical)', ...
          'FontSize', 15, 'FontWeight', 'bold');

    legend('Location', 'best', 'FontSize', 11);
    grid on; box on;
    xlim([u_min - margin, u_max + margin]); ylim([y(1), y(end)]);
    set(gca, 'FontSize', 11);

    output_path = fullfile(output_dir, 'velocity_evolution.png');
    print(output_path, '-dpng', '-r300');
    fprintf('   ✓ 已保存: %s\n', output_path);

    close(fig);
end

% 可视化函数2：8时刻对比图（2×4子图）
function plot_midplane_comparison(u_history, u_analytical, x, y, time_history, output_dir)
    % 功能：生成8个子图对比数值解与解析解

    u_min = min(cellfun(@(uu) min(uu(:)), u_history));
    u_max = max(cellfun(@(uu) max(uu(:)), u_history));
    margin = 0.05 * (u_max - u_min + eps);

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

    midplane_idx = round(length(x) / 2);

    fig = figure('Position', [100, 100, 1400, 900], 'Color', 'white');

    for i = 1:n_plots
        idx = time_indices(i);
        subplot(2, 4, i);

        u_midplane = u_history{idx}(midplane_idx, :);

        plot(u_midplane, y, 'b-o', 'LineWidth', 1.5, ...
             'MarkerSize', 4, 'MarkerFaceColor', 'b', 'DisplayName', 'Numerical');
        hold on;

        plot(u_analytical{idx}, y, 'r--', 'LineWidth', 2, ...
             'DisplayName', 'Analytical');

        xlabel('u/U', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('y/h', 'FontSize', 11, 'FontWeight', 'bold');
        title(sprintf('t = %.2f', time_history(idx)), 'FontSize', 12, 'FontWeight', 'bold');

        legend('Location', 'best', 'FontSize', 9);
        grid on; box on;
        xlim([u_min - margin, u_max + margin]); ylim([y(1), y(end)]);
        set(gca, 'FontSize', 10);
    end

    sgtitle('Midplane Velocity Profile vs Analytical Solution (x = \pi)', ...
            'FontSize', 16, 'FontWeight', 'bold');

    output_path = fullfile(output_dir, 'midplane_comparison.png');
    print(output_path, '-dpng', '-r300');
    fprintf('   ✓ 已保存: %s\n', output_path);

    close(fig);
end

% 可视化函数3：稳态2D流场
function plot_steady_state_2D(u_history, x, y, time_history, output_dir)
    % 功能：绘制稳态2D速度场云图

    idx_final = length(time_history);
    fig = figure('Position', [100, 100, 1000, 700], 'Color', 'white');

    Lx_local = x(end) - x(1);
    Ly_local = y(end) - y(1);
    [X, Y] = meshgrid(x / Lx_local, y / Ly_local);

    pcolor(X, Y, u_history{idx_final}');
    shading interp;

    h_colorbar = colorbar;
    ylabel(h_colorbar, 'u/U', 'FontSize', 13, 'FontWeight', 'bold');
    u_min = min(u_history{idx_final}(:));
    u_max = max(u_history{idx_final}(:));
    caxis([u_min, u_max]);
    colormap(jet);

    hold on;
    midplane_x = x(round(length(x)/2)) / Lx_local;
    plot([midplane_x, midplane_x], [0, 1], 'k--', 'LineWidth', 2.5);
    text(midplane_x + 0.05, 0.93, 'x = \pi', 'FontSize', 12, ...
         'Color', 'k', 'FontWeight', 'bold', 'BackgroundColor', 'white');

    xlabel('x/L_x', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('y/h', 'FontSize', 14, 'FontWeight', 'bold');
    title(sprintf('Steady State Velocity Field (t = %.2f)', time_history(idx_final)), ...
          'FontSize', 16, 'FontWeight', 'bold');

    axis equal tight;
    set(gca, 'FontSize', 12);

    output_path = fullfile(output_dir, 'steady_state_2D.png');
    print(output_path, '-dpng', '-r300');
    fprintf('   ✓ 已保存: %s\n', output_path);

    close(fig);
end

% 可视化函数4：动画生成（双面板布局）
function generate_animation(u_history, u_analytical, x, y, time_history, output_dir)
    % 功能：生成演化动画

    fprintf('   [动画] 准备数据...\n');

    n_frames = length(time_history);
    [X, Y] = meshgrid(x, y);
    x_mid_idx = round(length(x) / 2);

    u_max = max(cellfun(@(u) max(u(:)), u_history));
    u_min = min(cellfun(@(u) min(u(:)), u_history));

    fprintf('   [动画] 数据准备完成: %d 帧, 时间范围 [%.2f, %.2f]s\n', ...
            n_frames, time_history(1), time_history(end));

    fig = figure('Name', 'Couette Flow Animation', ...
                 'Position', [100, 100, 1400, 600], 'Color', 'w');

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

    hold(ax_main, 'on');
    plot(ax_main, [x(x_mid_idx), x(x_mid_idx)], [y(1), y(end)], ...
         'r--', 'LineWidth', 2);

    ax_profile = subplot(1, 2, 2);
    u_profile_numerical = u_history{1}(x_mid_idx, :)';
    u_profile_analytical = u_analytical{1};

    h_numerical = plot(ax_profile, u_profile_numerical, y, ...
                       'b-', 'LineWidth', 2.0, 'DisplayName', 'Numerical');
    hold(ax_profile, 'on');
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

    video_path = fullfile(output_dir, 'couette_animation.mp4');
    v = VideoWriter(video_path, 'MPEG-4');
    v.FrameRate = 20;
    v.Quality = 90;
    open(v);

    fprintf('   [动画] 视频编码器创建完成: %s\n', video_path);
    fprintf('   [动画] 帧率: %d fps, 预计时长: %.1f 秒\n', ...
            v.FrameRate, n_frames / v.FrameRate);

    fprintf('   [动画] 开始渲染帧...\n');

    for i = 1:n_frames
        set(h_pcolor, 'CData', u_history{i}');

        u_profile_numerical = u_history{i}(x_mid_idx, :)';
        u_profile_analytical = u_analytical{i};

        set(h_numerical, 'XData', u_profile_numerical);
        set(h_analytical, 'XData', u_profile_analytical);

        sgtitle(sprintf('Couette Flow Evolution | t = %.4f s | Frame %d/%d', ...
                        time_history(i), i, n_frames), ...
                'FontSize', 15, 'FontWeight', 'bold');

        drawnow;
        frame = getframe(fig);
        writeVideo(v, frame);

        if mod(i, 50) == 0 || i == n_frames
            fprintf('      进度: %d/%d (%.1f%%)\n', i, n_frames, 100*i/n_frames);
        end
    end

    close(v);
    close(fig);

    fprintf('   [动画] 渲染完成！\n');
    fprintf('   ✓ 已保存: %s\n', video_path);

    if exist(video_path, 'file')
        file_info = dir(video_path);
        fprintf('   [动画] 文件大小: %.2f MB\n', file_info.bytes / (1024*1024));
    end
end
