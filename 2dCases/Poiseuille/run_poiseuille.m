%% run_poiseuille.m
% Poiseuille 平板通道压力驱动二维不可压 NS 启动算例 - 计算脚本（非充分发展）
%
% 控制方程（二维不可压 Navier–Stokes，ρ=常数，默认ρ=1）:
%   ∂u/∂t + u∂u/∂x + v∂u/∂y = -(1/ρ)∂p/∂x + ν(∂²u/∂x²+∂²u/∂y²) + f_x
%   ∂v/∂t + u∂v/∂x + v∂v/∂y = -(1/ρ)∂p/∂y + ν(∂²v/∂x²+∂²v/∂y²) + f_y
%   ∂u/∂x + ∂v/∂y = 0
%
% 数值方法:
%   - 显式欧拉推进对流/二维粘性/体力得到中间速度 (u*, v*)
%   - 投影(分步)法解压力泊松方程并修正速度
%
% 边界条件:
%   - y=0 与 y=h: 无滑移/无渗透 u=v=0
%   - x 方向周期
%
% 本算例参数(用户指定):
%   Lx = 2*pi
%   Ly = h = 2.0
%   Re = 100 (基于 u_bulk 与 h)
%   u_bulk = 1
%
% 由 Re = u_bulk * h / nu 得:
%   nu = u_bulk * h / Re
% 为保证稳态平均速度为 u_bulk, 压力梯度等效体力取:
%   f_x = 12 * nu * u_bulk / h^2
%
% 初始条件:
%   u(x,y,0)=u_bulk (壁面强制为0), v=0, p=0
%   只计算到用户设定的 t_end（非充分发展阶段）

clear; clc; close all;

%% 添加 src 路径
case_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(case_dir, '..', '..', 'src')));

%% ========== Part 1: 参数设置 ==========
% 网格参数
Nx = 32;
Ny = 51;
Lx = 2*pi;
Ly = 2.0;           % 通道高度 h

% 物理参数（无量纲约束）
Re = 100;
u_bulk = 1.0;       % 体积平均速度
nu = u_bulk * Ly / Re;

% 常数源项（等效压力梯度体力）
source_term = 12.0 * nu * u_bulk / (Ly^2);

% 时间参数
t_end = 100;    % 只计算非充分发展阶段(未到稳态)
n_snapshots = 51;
cfl_diff = 0.45;            % 扩散 CFL(=alpha_y) 系数
cfl_adv = 0.5;              % 对流 CFL 系数

fprintf('========== Poiseuille 算例 ==========\n');
fprintf('Lx = %.6f, Ly = %.6f\n', Lx, Ly);
fprintf('Re = %.2f, u_bulk = %.2f\n', Re, u_bulk);
fprintf('nu = %.6f\n', nu);
fprintf('source_term S = %.6f\n', source_term);
fprintf('t_end = %.2f, snapshots = %d\n', t_end, n_snapshots);
fprintf('cfl_diff = %.2f (alpha_y target)\n', cfl_diff);
fprintf('cfl_adv  = %.2f (advection CFL)\n', cfl_adv);
fprintf('====================================\n\n');

% 网格生成
mesh = Mesh2D(Nx, Ny, Lx, Ly);
grid_info = mesh.to_struct();
x = grid_info.x;
y = grid_info.y;

physical_params = struct('U', 0.0, 'nu', nu, 'rho', 1.0, ...
                         'source_term_x', source_term, ...
                         'Re', Re, 'u_bulk', u_bulk);

% 估计速度尺度用于对流 CFL（取稳态解析最大值与初始/壁面速度的保守上界）
U_val = physical_params.U;
G_val = source_term;   % 体力/压力梯度等效项
h = Ly;
u_max_est = max(abs([u_bulk, U_val]));
if abs(G_val) > eps
    y_star = h/2 + (U_val * nu) / (G_val * h);  % 稳态极值位置
    y_star = min(max(y_star, 0.0), h);
    u_star = U_val * (y_star / h) + (G_val / (2 * nu)) * (h * y_star - y_star^2);
    u_max_est = max(u_max_est, abs(u_star));
end
v_max_est = 0.0;

% 根据对流-扩散 CFL 自动计算 dt（显式对流 + 二维 Laplacian 粘性项）
[dt, dt_details] = TimeStepCalculator.compute_explicit_advection_diffusion_dt( ...
    grid_info, nu, u_max_est, v_max_est, cfl_adv, cfl_diff);

fprintf('自动时间步长: dt = %.6f\n', dt);
fprintf('  dt_adv = %.6f, dt_diff = %.6f\n', dt_details.dt_adv, dt_details.dt_diff);
fprintf('  alpha_x=%.4f, alpha_y=%.4f (sum=%.4f)\n', ...
        dt_details.alpha_x, dt_details.alpha_y, dt_details.alpha_x + dt_details.alpha_y);
fprintf('  cfl_x=%.4f, cfl_y=%.4f\n\n', dt_details.cfl_x, dt_details.cfl_y);

time_params = struct('dt', dt, 't_end', t_end, 'n_snapshots', n_snapshots, ...
                     'cfl_diff', cfl_diff, 'cfl_adv', cfl_adv);

%% ========== Part 2: 构建求解器与边界条件 ==========
solver = FractionalStepNSSolver2D(grid_info, physical_params, time_params);

bc_u = CompositeBoundaryCondition();
bc_u.add_bc(DirichletBC('y', 'min', 0.0));  % 下壁无滑移 u=0
bc_u.add_bc(DirichletBC('y', 'max', 0.0));  % 上壁无滑移 u=0
bc_u.add_bc(PeriodicBC('x'));              % x 方向周期

bc_v = CompositeBoundaryCondition();
bc_v.add_bc(DirichletBC('y', 'min', 0.0));  % 下壁无渗透 v=0
bc_v.add_bc(DirichletBC('y', 'max', 0.0));  % 上壁无渗透 v=0
bc_v.add_bc(PeriodicBC('x'));              % x 方向周期

%% ========== Part 3: 初始条件与时间推进 ==========
% 初始条件：从 u=1 的均匀场开始（随后在无滑移壁面/压力梯度作用下演化）
u0 = u_bulk * ones(Nx, Ny);
v0 = zeros(Nx, Ny);
u0 = bc_u.apply_all(u0, grid_info, physical_params);
v0 = bc_v.apply_all(v0, grid_info, physical_params);

[u_history, v_history, p_history, time_history] = solver.run_simulation(u0, v0, bc_u, bc_v);

alpha_y = nu * dt / (grid_info.dy^2);

%% ========== Part 4: 保存数据 ==========
output_dir = fullfile(case_dir, 'results');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

save(fullfile(output_dir, 'poiseuille_results.mat'), ...
     'x', 'y', 'u_history', 'v_history', 'p_history', 'time_history', ...
     'Nx', 'Ny', 'Lx', 'Ly', 'Re', 'u_bulk', 'nu', 'dt', 'alpha_y', ...
     't_end', 'n_snapshots', 'source_term');

fprintf('✓ 数据已保存: %s\n', fullfile(output_dir, 'poiseuille_results.mat'));
