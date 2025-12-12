classdef FractionalStepNSSolver2D < handle
    % FractionalStepNSSolver2D 二维不可压 Navier-Stokes 分步(投影)求解器
    %
    % 求解方程(ρ=常数，默认ρ=1):
    %   ∂u/∂t + u∂u/∂x + v∂u/∂y = -1/ρ ∂p/∂x + ν(∂²u/∂x²+∂²u/∂y²) + f_x
    %   ∂v/∂t + u∂v/∂x + v∂v/∂y = -1/ρ ∂p/∂y + ν(∂²v/∂x²+∂²v/∂y²) + f_y
    %   ∂u/∂x + ∂v/∂y = 0
    %
    % 数值方法:
    %   - 时间: 显式欧拉
    %   - 空间: 二阶中心差分(对流/扩散)
    %   - 压力: 分步法(Chorin) + SOR 解 Poisson
    %
    % 适用网格:
    %   Mesh2D 结构化均匀网格，周期方向使用 ghost cell (PeriodicBC 同约定)

    properties (Access = private)
        grid_info
        physical_params
        time_params

        nu
        rho
        fx_const
        fy_const

        poisson_omega = 1.7
        poisson_max_iter = 5000
        poisson_tol = 1e-6
    end

    methods
        function obj = FractionalStepNSSolver2D(grid_info, physical_params, time_params)
            obj.grid_info = grid_info;
            obj.physical_params = physical_params;
            obj.time_params = time_params;

            obj.nu = physical_params.nu;
            if isfield(physical_params, 'rho')
                obj.rho = physical_params.rho;
            else
                obj.rho = 1.0;
            end

            % 兼容不同命名的体力/源项：默认作用在 x 方向
            obj.fx_const = 0.0;
            obj.fy_const = 0.0;
            if isfield(physical_params, 'source_term_x')
                obj.fx_const = physical_params.source_term_x;
            elseif isfield(physical_params, 'source_term')
                obj.fx_const = physical_params.source_term;
            elseif isfield(physical_params, 'G')
                obj.fx_const = physical_params.G;
            end
            if isfield(physical_params, 'source_term_y')
                obj.fy_const = physical_params.source_term_y;
            end

            if isfield(time_params, 'poisson_omega'); obj.poisson_omega = time_params.poisson_omega; end
            if isfield(time_params, 'poisson_max_iter'); obj.poisson_max_iter = time_params.poisson_max_iter; end
            if isfield(time_params, 'poisson_tol'); obj.poisson_tol = time_params.poisson_tol; end

            fprintf('========== 2D NS 分步求解器 ==========\n');
            fprintf('网格: %d × %d\n', grid_info.Nx, grid_info.Ny);
            fprintf('dt = %.6f, t_end = %.3f\n', time_params.dt, time_params.t_end);
            fprintf('dx = %.6f, dy = %.6f\n', grid_info.dx, grid_info.dy);
            fprintf('nu = %.6f, rho = %.3f\n', obj.nu, obj.rho);
            fprintf('body force: fx = %.6f, fy = %.6f\n', obj.fx_const, obj.fy_const);
            fprintf('Poisson: omega=%.2f, max_iter=%d, tol=%.1e\n', ...
                obj.poisson_omega, obj.poisson_max_iter, obj.poisson_tol);
            fprintf('====================================\n\n');
        end

        function [u_hist, v_hist, p_hist, t_hist] = run_simulation(obj, u0, v0, bc_u, bc_v, p0)
            if nargin < 6 || isempty(p0)
                p0 = zeros(obj.grid_info.Nx, obj.grid_info.Ny);
            end

            Nx = obj.grid_info.Nx;
            Ny = obj.grid_info.Ny;
            dt = obj.time_params.dt;
            t_end = obj.time_params.t_end;
            n_snap = obj.time_params.n_snapshots;

            output_interval = t_end / (n_snap - 1);

            u_hist = cell(1, n_snap);
            v_hist = cell(1, n_snap);
            p_hist = cell(1, n_snap);
            t_hist = zeros(1, n_snap);

            % 初始边界条件
            u = bc_u.apply_all(u0, obj.grid_info, obj.physical_params);
            v = bc_v.apply_all(v0, obj.grid_info, obj.physical_params);
            p = obj.apply_pressure_bc(p0);

            snap_idx = 1;
            u_hist{snap_idx} = u;
            v_hist{snap_idx} = v;
            p_hist{snap_idx} = p;
            t_hist(snap_idx) = 0.0;
            next_out = output_interval;

            current_time = 0.0;
            step = 0;

            fprintf('开始时间推进...\n');
            fprintf('------------------------------------\n');
            tic;

            while current_time < t_end - 1e-12
                step = step + 1;

                [u, v, p] = obj.step_once(u, v, p, bc_u, bc_v);

                current_time = current_time + dt;

                if current_time >= next_out - 1e-10 && snap_idx < n_snap
                    snap_idx = snap_idx + 1;
                    u_hist{snap_idx} = u;
                    v_hist{snap_idx} = v;
                    p_hist{snap_idx} = p;
                    t_hist(snap_idx) = current_time;
                    next_out = snap_idx * output_interval;
                end

                if mod(step, 5000) == 0
                    fprintf('进度: %.1f%% (t=%.3fs, 步数=%d)\n', ...
                        current_time / t_end * 100, current_time, step);
                end
            end

            if snap_idx < n_snap
                snap_idx = snap_idx + 1;
                u_hist{snap_idx} = u;
                v_hist{snap_idx} = v;
                p_hist{snap_idx} = p;
                t_hist(snap_idx) = current_time;
            end

            elapsed_time = toc;
            fprintf('------------------------------------\n');
            fprintf('时间推进完成！耗时: %.2f 秒\n', elapsed_time);
            fprintf('实际计算步数: %d\n', step);
            fprintf('平均速度: %.0f 步/秒\n\n', step / elapsed_time);
        end
    end

    methods (Access = private)
        function [u_new, v_new, p_new] = step_once(obj, u, v, p, bc_u, bc_v)
            % 1) 应用速度边界条件(含 periodic ghost)
            u = bc_u.apply_all(u, obj.grid_info, obj.physical_params);
            v = bc_v.apply_all(v, obj.grid_info, obj.physical_params);

            Nx = obj.grid_info.Nx;
            Ny = obj.grid_info.Ny;
            dx = obj.grid_info.dx;
            dy = obj.grid_info.dy;
            dt = obj.time_params.dt;
            nu = obj.nu;

            % 2) 计算对流/扩散项（只更新内部点 2:Nx-1,2:Ny-1）
            du_dx = (u(3:Nx, 2:Ny-1) - u(1:Nx-2, 2:Ny-1)) / (2*dx);
            du_dy = (u(2:Nx-1, 3:Ny) - u(2:Nx-1, 1:Ny-2)) / (2*dy);
            dv_dx = (v(3:Nx, 2:Ny-1) - v(1:Nx-2, 2:Ny-1)) / (2*dx);
            dv_dy = (v(2:Nx-1, 3:Ny) - v(2:Nx-1, 1:Ny-2)) / (2*dy);

            conv_u = u(2:Nx-1, 2:Ny-1) .* du_dx + v(2:Nx-1, 2:Ny-1) .* du_dy;
            conv_v = u(2:Nx-1, 2:Ny-1) .* dv_dx + v(2:Nx-1, 2:Ny-1) .* dv_dy;

            lap_u = (u(3:Nx, 2:Ny-1) - 2*u(2:Nx-1, 2:Ny-1) + u(1:Nx-2, 2:Ny-1)) / dx^2 + ...
                    (u(2:Nx-1, 3:Ny) - 2*u(2:Nx-1, 2:Ny-1) + u(2:Nx-1, 1:Ny-2)) / dy^2;

            lap_v = (v(3:Nx, 2:Ny-1) - 2*v(2:Nx-1, 2:Ny-1) + v(1:Nx-2, 2:Ny-1)) / dx^2 + ...
                    (v(2:Nx-1, 3:Ny) - 2*v(2:Nx-1, 2:Ny-1) + v(2:Nx-1, 1:Ny-2)) / dy^2;

            u_star = u;
            v_star = v;
            u_star(2:Nx-1, 2:Ny-1) = u(2:Nx-1, 2:Ny-1) + dt * (-conv_u + nu*lap_u + obj.fx_const);
            v_star(2:Nx-1, 2:Ny-1) = v(2:Nx-1, 2:Ny-1) + dt * (-conv_v + nu*lap_v + obj.fy_const);

            % 3) 施加速度边界到中间速度
            u_star = bc_u.apply_all(u_star, obj.grid_info, obj.physical_params);
            v_star = bc_v.apply_all(v_star, obj.grid_info, obj.physical_params);

            % 4) 组装压力 Poisson 右端项 rhs = rho/dt * div(u_star)
            du_star_dx = (u_star(3:Nx, 2:Ny-1) - u_star(1:Nx-2, 2:Ny-1)) / (2*dx);
            dv_star_dy = (v_star(2:Nx-1, 3:Ny) - v_star(2:Nx-1, 1:Ny-2)) / (2*dy);
            rhs = (obj.rho / dt) * (du_star_dx + dv_star_dy);  % 尺寸 (Nx-2)×(Ny-2)

            % 5) 解 Poisson 得到 p^{n+1}
            p_new = obj.solve_poisson(p, rhs);

            % 6) 速度修正
            dp_dx = (p_new(3:Nx, 2:Ny-1) - p_new(1:Nx-2, 2:Ny-1)) / (2*dx);
            dp_dy = (p_new(2:Nx-1, 3:Ny) - p_new(2:Nx-1, 1:Ny-2)) / (2*dy);

            u_new = u_star;
            v_new = v_star;
            u_new(2:Nx-1, 2:Ny-1) = u_star(2:Nx-1, 2:Ny-1) - dt/obj.rho * dp_dx;
            v_new(2:Nx-1, 2:Ny-1) = v_star(2:Nx-1, 2:Ny-1) - dt/obj.rho * dp_dy;

            % 7) 最终速度边界
            u_new = bc_u.apply_all(u_new, obj.grid_info, obj.physical_params);
            v_new = bc_v.apply_all(v_new, obj.grid_info, obj.physical_params);
        end

        function p = solve_poisson(obj, p, rhs)
            % solve_poisson SOR 解压力 Poisson:
            %   ∇²p = rhs
            Nx = obj.grid_info.Nx;
            Ny = obj.grid_info.Ny;
            dx2 = obj.grid_info.dx^2;
            dy2 = obj.grid_info.dy^2;
            coef = 1.0 / (2*(dx2 + dy2));

            omega = obj.poisson_omega;
            tol = obj.poisson_tol;

            for iter = 1:obj.poisson_max_iter
                for i = 2:Nx-1
                    for j = 2:Ny-1
                        rhs_ij = rhs(i-1, j-1);
                        p_new = ((p(i+1,j) + p(i-1,j))*dy2 + (p(i,j+1) + p(i,j-1))*dx2 ...
                                 - rhs_ij*dx2*dy2) * coef;
                        p(i,j) = (1-omega)*p(i,j) + omega*p_new;
                    end
                end

                p = obj.apply_pressure_bc(p);
                p(2,2) = 0.0; % 固定参考压力防止奇异

                if mod(iter, 50) == 0
                    lap_p = (p(3:Nx, 2:Ny-1) - 2*p(2:Nx-1, 2:Ny-1) + p(1:Nx-2, 2:Ny-1)) / dx2 + ...
                            (p(2:Nx-1, 3:Ny) - 2*p(2:Nx-1, 2:Ny-1) + p(2:Nx-1, 1:Ny-2)) / dy2;
                    res = max(abs(lap_p - rhs), [], 'all');
                    if res < tol
                        break;
                    end
                end
            end
        end

        function p = apply_pressure_bc(obj, p)
            % apply_pressure_bc 压力边界:
            %   x 方向周期(ghost)
            %   y 方向 Neumann dp/dy=0
            Nx = obj.grid_info.Nx;
            Ny = obj.grid_info.Ny;
            % 周期 x
            p(1, :) = p(Nx-1, :);
            p(Nx, :) = p(2, :);
            % Neumann y
            p(:, 1) = p(:, 2);
            p(:, Ny) = p(:, Ny-1);
        end
    end
end

