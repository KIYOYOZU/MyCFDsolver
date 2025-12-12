%% postprocess_poiseuille.m
% Poiseuille 启动过程后处理脚本
%
% 说明：
%   - 读取 run_poiseuille.m 生成的 results/poiseuille_results.mat
%   - 计算 Poiseuille 稳态解析解并叠加到数值结果上（不画解析瞬态演化）
%   - 生成 3 张静态图 + 1 个动画

clear; clc; close all;

case_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(case_dir, '..', '..', 'src')));  % 通用模块路径

% CouetteAnalytical 在 Couette 算例目录中，Poiseuille 复用其解析解实现
analytical_dir = fullfile(case_dir, '..', 'coutte_flow');
if exist(fullfile(analytical_dir, 'CouetteAnalytical.m'), 'file')
    addpath(analytical_dir);
end

output_dir = fullfile(case_dir, 'results');
data_path = fullfile(output_dir, 'poiseuille_results.mat');

if ~exist(data_path, 'file')
    error('未找到结果文件 %s，请先运行 run_poiseuille.m', data_path);
end

load(data_path);

fprintf('开始后处理...\n');
fprintf('输出目录: %s\n', output_dir);

%% ========== 解析解与误差 ==========
% 兼容旧数据：可能没有 source_term 字段
G_val = 0.0;
if exist('source_term', 'var')
    G_val = source_term;
elseif exist('G', 'var')
    G_val = G;
end

if ~exist('u_bulk', 'var')
    u_bulk = 1.0;  % 兜底，用于归一化显示
end

U_val = 0.0;
if exist('U', 'var')
    U_val = U;
end

if ~exist('Nx', 'var'); Nx = size(u_history{1}, 1); end
if ~exist('Ny', 'var'); Ny = size(u_history{1}, 2); end
n_snapshots = length(time_history);

% 稳态解析解（Poiseuille-Couette 稳态分量）
% u_s(y) = U*(y/h) + (G/(2*nu))*(h*y - y^2)
y_col = y(:);
u_steady_analytical = U_val * (y_col / Ly) + (G_val / (2 * nu)) * (Ly * y_col - y_col.^2);

L2_history = zeros(1, n_snapshots);

for i = 1:n_snapshots
    u_ana_2D = repmat(u_steady_analytical', Nx, 1);
    L2_history(i) = ErrorAnalyzer.compute_L2(u_history{i}, u_ana_2D);
end

fprintf('\n误差分析:\n');
fprintf('  - L2误差范围: [%.4e, %.4e]\n', min(L2_history), max(L2_history));
fprintf('  - 最终L2误差: %.4e\n\n', L2_history(end));

% 三张静态图（解析线只取稳态）
plot_velocity_evolution(u_history, u_steady_analytical, x, y, time_history, output_dir, u_bulk);
plot_midplane_comparison(u_history, u_steady_analytical, x, y, time_history, output_dir, u_bulk);
plot_steady_state_2D(u_history, x, y, time_history, output_dir, u_bulk);

% 动画
generate_animation(u_history, u_steady_analytical, x, y, time_history, output_dir, u_bulk);

fprintf('后处理完成！\n');

%% ========== 子函数 ==========

% 可视化函数1：速度剖面演化对比
function plot_velocity_evolution(u_history, u_steady_analytical, x, y, time_history, output_dir, u_bulk)
    fig = figure('Name', 'Velocity Profile Evolution', ...
                 'Position', [100, 100, 1000, 700], 'Color', 'w');

    midplane_idx = round(length(x) / 2);

    u_min = min(cellfun(@(uu) min(uu(:)), u_history));
    u_max = max(cellfun(@(uu) max(uu(:)), u_history));
    margin = 0.05 * (u_max - u_min + eps);
    margin_n = margin / u_bulk;

    num_frames = length(time_history);
    frame_indices = unique(round([1, 0.25*num_frames, 0.5*num_frames, ...
                                  0.75*num_frames, num_frames]));

    colors = lines(length(frame_indices));
    hold on;

    plot(u_steady_analytical / u_bulk, y, 'k--', 'LineWidth', 2.5, ...
         'DisplayName', 'Steady State (analytical)');

    for j = 1:length(frame_indices)
        idx = frame_indices(j);
        t = time_history(idx);
        u_profile = u_history{idx}(midplane_idx, :);
        plot(u_profile / u_bulk, y, '-', 'Color', colors(j,:), 'LineWidth', 1.8, ...
             'DisplayName', sprintf('t = %.2f (numerical)', t));
    end

    xlabel('u/u_{bulk}', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('y/h', 'FontSize', 13, 'FontWeight', 'bold');
    title('Velocity Profile Evolution (Numerical vs Analytical)', ...
          'FontSize', 15, 'FontWeight', 'bold');

    legend('Location', 'best', 'FontSize', 11);
    grid on; box on;
    xlim([u_min/u_bulk - margin_n, u_max/u_bulk + margin_n]); ylim([y(1), y(end)]);
    set(gca, 'FontSize', 11);

    output_path = fullfile(output_dir, 'velocity_evolution.png');
    print(output_path, '-dpng', '-r300');
    fprintf('   ✓ 已保存: %s\n', output_path);

    close(fig);
end

% 可视化函数2：8时刻对比图
function plot_midplane_comparison(u_history, u_steady_analytical, x, y, time_history, output_dir, u_bulk)
    u_min = min(cellfun(@(uu) min(uu(:)), u_history));
    u_max = max(cellfun(@(uu) max(uu(:)), u_history));
    margin = 0.05 * (u_max - u_min + eps);
    margin_n = margin / u_bulk;

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

        plot(u_midplane / u_bulk, y, 'b-o', 'LineWidth', 1.5, ...
             'MarkerSize', 4, 'MarkerFaceColor', 'b', 'DisplayName', 'Numerical');
        hold on;

        plot(u_steady_analytical / u_bulk, y, 'r--', 'LineWidth', 2, ...
             'DisplayName', 'Analytical');

        xlabel('u/u_{bulk}', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('y/h', 'FontSize', 11, 'FontWeight', 'bold');
        title(sprintf('t = %.2f', time_history(idx)), 'FontSize', 12, 'FontWeight', 'bold');

        legend('Location', 'best', 'FontSize', 9);
        grid on; box on;
        xlim([u_min/u_bulk - margin_n, u_max/u_bulk + margin_n]); ylim([y(1), y(end)]);
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
function plot_steady_state_2D(u_history, x, y, time_history, output_dir, u_bulk)
    idx_final = length(time_history);
    fig = figure('Position', [100, 100, 1000, 700], 'Color', 'white');

    Lx_local = x(end) - x(1);
    Ly_local = y(end) - y(1);
    [X, Y] = meshgrid(x / Lx_local, y / Ly_local);

    pcolor(X, Y, (u_history{idx_final}'/u_bulk));
    shading interp;

    h_colorbar = colorbar;
    ylabel(h_colorbar, 'u/u_{bulk}', 'FontSize', 13, 'FontWeight', 'bold');
    u_min = min(u_history{idx_final}(:));
    u_max = max(u_history{idx_final}(:));
    caxis([u_min/u_bulk, u_max/u_bulk]);
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

% 可视化函数4：动画生成
function generate_animation(u_history, u_steady_analytical, x, y, time_history, output_dir, u_bulk)
    fprintf('   [动画] 准备数据...\n');

    n_frames = length(time_history);
    [X, Y] = meshgrid(x, y);
    x_mid_idx = round(length(x) / 2);

    u_max = max(cellfun(@(u) max(u(:)), u_history));
    u_min = min(cellfun(@(u) min(u(:)), u_history));

    fprintf('   [动画] 数据准备完成: %d 帧, 时间范围 [%.2f, %.2f]s\n', ...
            n_frames, time_history(1), time_history(end));

    fig = figure('Name', 'Poiseuille Flow Animation', ...
                 'Position', [100, 100, 1400, 600], 'Color', 'w');

    ax_main = subplot(1, 2, 1);
    h_pcolor = pcolor(ax_main, X, Y, (u_history{1}'/u_bulk));
    shading(ax_main, 'interp');
    colormap(ax_main, 'jet');
    h_colorbar = colorbar(ax_main);
    h_colorbar.Label.String = 'u/u_{bulk}';
    h_colorbar.Label.FontSize = 11;
    caxis(ax_main, [u_min/u_bulk, u_max/u_bulk]);
    xlabel(ax_main, 'x/(2\pi)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel(ax_main, 'y/h', 'FontSize', 12, 'FontWeight', 'bold');
    title(ax_main, 'Velocity Field', 'FontSize', 13, 'FontWeight', 'bold');
    axis(ax_main, 'equal', 'tight');

    hold(ax_main, 'on');
    plot(ax_main, [x(x_mid_idx), x(x_mid_idx)], [y(1), y(end)], ...
         'r--', 'LineWidth', 2);

    ax_profile = subplot(1, 2, 2);
    u_profile_numerical = u_history{1}(x_mid_idx, :)' / u_bulk;
    u_profile_analytical = u_steady_analytical / u_bulk;

    h_numerical = plot(ax_profile, u_profile_numerical, y, ...
                       'b-', 'LineWidth', 2.0, 'DisplayName', 'Numerical');
    hold(ax_profile, 'on');
    h_analytical = plot(ax_profile, u_profile_analytical, y, ...
                        'r--', 'LineWidth', 2.0, 'DisplayName', 'Analytical');

    xlabel(ax_profile, 'u/u_{bulk}', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel(ax_profile, 'y/h', 'FontSize', 12, 'FontWeight', 'bold');
    title(ax_profile, 'Midplane Profile', 'FontSize', 13, 'FontWeight', 'bold');
    legend(ax_profile, 'Location', 'best');
    grid(ax_profile, 'on'); box(ax_profile, 'on');
    xlim(ax_profile, [u_min/u_bulk, u_max/u_bulk]);
    ylim(ax_profile, [y(1), y(end)]);

    video_path = fullfile(output_dir, 'poiseuille_animation.mp4');
    writerObj = VideoWriter(video_path, 'MPEG-4');
    writerObj.FrameRate = 20;
    open(writerObj);

    for k = 1:n_frames
        set(h_pcolor, 'CData', u_history{k}'/u_bulk);
        set(h_numerical, 'XData', u_history{k}(x_mid_idx, :)'/u_bulk);
        % 解析线为稳态，不随时间变化

        title(ax_main, sprintf('Velocity Field (t = %.2f s)', time_history(k)), ...
              'FontSize', 13, 'FontWeight', 'bold');
        title(ax_profile, sprintf('Midplane Profile (t = %.2f s)', time_history(k)), ...
              'FontSize', 13, 'FontWeight', 'bold');

        drawnow;
        frame = getframe(fig);
        writeVideo(writerObj, frame);
    end

    close(writerObj);
    fprintf('   ✓ 动画已保存: %s\n', video_path);

    close(fig);
end
