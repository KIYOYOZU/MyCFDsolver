%% run_couette.m
% Couette流启动过程数值模拟 - 模块化计算脚本
%
% 说明：
%   - 本脚本只负责“计算 + 保存数据”
%   - 可视化与动画请运行 postprocess_couette.m

clear; clc; close all;

%% 添加 src 路径
case_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(case_dir, '..', '..', 'src')));

%% ========== Part 1: 参数设置 ==========
% 物理参数
Nx = 32;
Ny = 51;
Lx = 2*pi;
Ly = 2.0;        % 通道高度 h
U = 1.0;         % 上壁速度
nu = 0.1;        % 运动粘度
source_term = 1; % 常数源项 S（等效 x方向压力梯度体力），S = -(1/ρ)·dp/dx

% 时间参数
t_end = 20.0;
n_snapshots = 51;
cfl_diff = 0.45;   % 扩散 CFL(=alpha_y) 系数
cfl_adv = 0.5;     % 对流 CFL 系数

% 网格生成（通用 mesh 模块）
mesh = Mesh2D(Nx, Ny, Lx, Ly);
grid_info = mesh.to_struct();
x = grid_info.x;
y = grid_info.y;

% 估计稳态最大速度用于对流 CFL（Poiseuille-Couette 稳态解析上界）
U_val = U;
G_val = source_term;
h = Ly;
u_max_est = abs(U_val);
if abs(G_val) > eps
    y_star = h/2 + (U_val * nu) / (G_val * h);
    y_star = min(max(y_star, 0.0), h);
    u_star = U_val * (y_star / h) + (G_val / (2 * nu)) * (h * y_star - y_star^2);
    u_max_est = max(u_max_est, abs(u_star));
end
v_max_est = 0.0;

[dt, dt_details] = TimeStepCalculator.compute_explicit_advection_diffusion_dt( ...
    grid_info, nu, u_max_est, v_max_est, cfl_adv, cfl_diff);

fprintf('自动时间步长: dt = %.6f\n', dt);
fprintf('  dt_adv = %.6f, dt_diff = %.6f\n', dt_details.dt_adv, dt_details.dt_diff);
fprintf('  alpha_x=%.4f, alpha_y=%.4f (sum=%.4f)\n', ...
        dt_details.alpha_x, dt_details.alpha_y, dt_details.alpha_x + dt_details.alpha_y);
fprintf('  cfl_x=%.4f, cfl_y=%.4f\n\n', dt_details.cfl_x, dt_details.cfl_y);

physical_params = struct('U', U, 'nu', nu, 'rho', 1.0, ...
                         'source_term_x', source_term);
time_params = struct('dt', dt, 't_end', t_end, 'n_snapshots', n_snapshots, ...
                     'cfl_diff', cfl_diff, 'cfl_adv', cfl_adv);

%% ========== Part 2: 构建求解器与边界条件 ==========
solver = FractionalStepNSSolver2D(grid_info, physical_params, time_params);

bc_u = CompositeBoundaryCondition();
bc_u.add_bc(DirichletBC('y', 'min', 0));   % 下壁无滑移 u=0
bc_u.add_bc(DirichletBC('y', 'max', U));   % 上壁 u=U
bc_u.add_bc(PeriodicBC('x'));             % x 方向周期

bc_v = CompositeBoundaryCondition();
bc_v.add_bc(DirichletBC('y', 'min', 0));   % 下壁 v=0
bc_v.add_bc(DirichletBC('y', 'max', 0));   % 上壁 v=0
bc_v.add_bc(PeriodicBC('x'));

%% ========== Part 3: 初始条件与时间推进 ==========
u0 = zeros(Nx, Ny);
v0 = zeros(Nx, Ny);
u0 = bc_u.apply_all(u0, grid_info, physical_params);
v0 = bc_v.apply_all(v0, grid_info, physical_params);

[u_history, v_history, p_history, time_history] = solver.run_simulation(u0, v0, bc_u, bc_v);

%% ========== Part 4: 保存数据 ==========
output_dir = fullfile(case_dir, 'results');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

alpha_y = nu * dt / (grid_info.dy^2);

save(fullfile(output_dir, 'couette_results.mat'), ...
     'x', 'y', 'u_history', 'v_history', 'p_history', 'time_history', ...
     'Nx', 'Ny', 'Lx', 'Ly', 'U', 'nu', 'dt', 'alpha_y', 't_end', 'n_snapshots', ...
     'source_term');

fprintf('✓ 数据已保存: %s\n', fullfile(output_dir, 'couette_results.mat'));
