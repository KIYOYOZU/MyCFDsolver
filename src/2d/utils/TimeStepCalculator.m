classdef TimeStepCalculator
    % TimeStepCalculator 时间步长计算工具
    %
    % 用于显式时间推进格式的稳定时间步长估算。
    %
    % 支持三类常见显式稳定时间步长:
    % 1) 1D/单方向扩散(当前 ExplicitEulerSolver 的形式):
    %    ∂u/∂t = ν ∂²u/∂y² + S
    %    alpha_y = ν dt / dy^2 <= 0.5
    %    dt = cfl * dy^2 / ν, 0 < cfl <= 0.5
    %
    % 2) 2D 扩散(五点 Laplacian):
    %    ∂u/∂t = ν(∂²u/∂x² + ∂²u/∂y²)
    %    alpha_x + alpha_y <= 0.5
    %    dt = cfl / (ν(1/dx^2 + 1/dy^2)), 0 < cfl <= 0.5
    %
    % 3) 2D 对流 CFL:
    %    dt <= cfl_adv * min(dx/|u|max, dy/|v|max), 0 < cfl_adv <= 1
    %
    % 组合对流-扩散:
    %    dt = min(dt_adv, dt_diff)

    methods (Static)
        function [dt, alpha_y] = compute_explicit_diffusion_dt(grid_info, nu, cfl)
            % compute_explicit_diffusion_dt 计算显式扩散稳定时间步长
            %
            % 输入:
            %   grid_info - Mesh2D.to_struct() 返回的网格结构体，需包含 dy
            %   nu        - 运动粘度 ν
            %   cfl       - 扩散 CFL(=alpha_y) 系数，默认 0.45
            %
            % 输出:
            %   dt        - 稳定时间步长
            %   alpha_y   - 对应扩散数 ν dt / dy^2

            if nargin < 3 || isempty(cfl)
                cfl = 0.45;
            end

            if cfl <= 0 || cfl > 0.5
                error('TimeStepCalculator:cfl', ...
                    'cfl must be in (0, 0.5] for explicit diffusion stability.');
            end
            if nu <= 0
                error('TimeStepCalculator:nu', 'nu must be positive.');
            end
            if ~isfield(grid_info, 'dy')
                error('TimeStepCalculator:grid', 'grid_info must contain dy.');
            end

            dy = grid_info.dy;
            dt = cfl * dy^2 / nu;
            alpha_y = nu * dt / dy^2;
        end

        function [dt, alpha_x, alpha_y] = compute_explicit_diffusion_dt2D(grid_info, nu, cfl)
            % compute_explicit_diffusion_dt2D 计算 2D 显式扩散稳定时间步长
            %
            % 输入:
            %   grid_info - 需包含 dx, dy
            %   nu        - 运动粘度 ν
            %   cfl       - 扩散 CFL(=alpha_x+alpha_y) 系数，默认 0.45
            %
            % 输出:
            %   dt        - 稳定时间步长
            %   alpha_x   - ν dt / dx^2
            %   alpha_y   - ν dt / dy^2

            if nargin < 3 || isempty(cfl)
                cfl = 0.45;
            end

            if cfl <= 0 || cfl > 0.5
                error('TimeStepCalculator:cfl', ...
                    'cfl must be in (0, 0.5] for explicit diffusion stability.');
            end
            if nu <= 0
                error('TimeStepCalculator:nu', 'nu must be positive.');
            end
            if ~isfield(grid_info, 'dx') || ~isfield(grid_info, 'dy')
                error('TimeStepCalculator:grid', 'grid_info must contain dx and dy.');
            end

            dx = grid_info.dx;
            dy = grid_info.dy;

            dt = cfl / (nu * (1/dx^2 + 1/dy^2));
            alpha_x = nu * dt / dx^2;
            alpha_y = nu * dt / dy^2;
        end

        function [dt, cfl_x, cfl_y] = compute_explicit_advection_dt(grid_info, u_max, v_max, cfl_adv)
            % compute_explicit_advection_dt 计算显式对流 CFL 时间步长
            %
            % 输入:
            %   grid_info - 需包含 dx, dy
            %   u_max     - |u| 最大值(标量)
            %   v_max     - |v| 最大值(标量)
            %   cfl_adv   - 对流 CFL 系数，默认 0.5, 0 < cfl_adv <= 1
            %
            % 输出:
            %   dt        - 稳定时间步长（若 u_max=v_max=0, 返回 Inf）
            %   cfl_x     - 实际 x 方向 CFL = |u|max dt/dx
            %   cfl_y     - 实际 y 方向 CFL = |v|max dt/dy

            if nargin < 4 || isempty(cfl_adv)
                cfl_adv = 0.5;
            end

            if cfl_adv <= 0 || cfl_adv > 1.0
                error('TimeStepCalculator:cfl_adv', ...
                    'cfl_adv must be in (0, 1] for explicit advection CFL.');
            end
            if ~isfield(grid_info, 'dx') || ~isfield(grid_info, 'dy')
                error('TimeStepCalculator:grid', 'grid_info must contain dx and dy.');
            end

            dx = grid_info.dx;
            dy = grid_info.dy;
            umax = abs(u_max);
            vmax = abs(v_max);

            dt_x = Inf;
            dt_y = Inf;
            if umax > eps
                dt_x = dx / umax;
            end
            if vmax > eps
                dt_y = dy / vmax;
            end

            dt = cfl_adv * min(dt_x, dt_y);

            if isinf(dt)
                cfl_x = 0; cfl_y = 0;
            else
                cfl_x = umax * dt / dx;
                cfl_y = vmax * dt / dy;
            end
        end

        function [dt, details] = compute_explicit_advection_diffusion_dt(grid_info, nu, u_max, v_max, cfl_adv, cfl_diff)
            % compute_explicit_advection_diffusion_dt 对流-扩散显式稳定 dt
            %
            % dt = min(dt_adv, dt_diff)
            % 返回 dt 及各分量信息，便于诊断。

            [dt_adv, cfl_x, cfl_y] = TimeStepCalculator.compute_explicit_advection_dt(grid_info, u_max, v_max, cfl_adv);
            [dt_diff, alpha_x, alpha_y] = TimeStepCalculator.compute_explicit_diffusion_dt2D(grid_info, nu, cfl_diff);

            dt = min(dt_adv, dt_diff);
            details = struct( ...
                'dt_adv', dt_adv, 'dt_diff', dt_diff, ...
                'cfl_x', cfl_x, 'cfl_y', cfl_y, ...
                'alpha_x', alpha_x, 'alpha_y', alpha_y);
        end
    end
end
