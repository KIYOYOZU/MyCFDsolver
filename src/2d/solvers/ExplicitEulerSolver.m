classdef ExplicitEulerSolver < BaseSolver
    % ExplicitEulerSolver 显式欧拉求解器
    %
    % 求解扩散方程（可含常数源项）:
    %   ∂u/∂t = ν · ∂²u/∂y² + S
    % 其中 S 为等效体力/源项（常数或函数句柄返回常数）。
    % 为向后兼容：
    %   - physical_params.source_term 优先生效
    %   - 若未提供 source_term，则读取 physical_params.G
    %   - 两者都没有则 S = 0（退化为 Couette 启动问题）
    % 离散:
    %   u^(n+1) = u^n + α_y (u_{j+1} - 2u_j + u_{j-1})

    properties (Access = private)
        alpha_y  % Y方向扩散数 α_y = ν dt / dy^2
        source_term_value  % 常数源项 S
        use_parallel = false  % 是否使用 parfor 并行加速（可选）
    end

    methods
        function obj = ExplicitEulerSolver(grid_info, physical_params, time_params)
            obj@BaseSolver(grid_info, physical_params, time_params);

            obj.alpha_y = physical_params.nu * time_params.dt / (grid_info.dy^2);

            % 可选并行配置：time_params.use_parallel = true
            if isfield(time_params, 'use_parallel')
                obj.use_parallel = logical(time_params.use_parallel);
            end

            % 若启用并行但无并行工具箱，则回退为串行
            if obj.use_parallel
                try
                    if ~license('test', 'Distrib_Computing_Toolbox')
                        warning('ExplicitEulerSolver:ParallelUnavailable', ...
                            'use_parallel=true 但未检测到 Parallel Computing Toolbox，回退为串行。');
                        obj.use_parallel = false;
                    else
                        % 启动并行池（若尚未启动）
                        if isempty(gcp('nocreate'))
                            % 默认启动本地并行池（兼容旧版本 MATLAB）
                            parpool;
                        end
                    end
                catch
                    warning('ExplicitEulerSolver:ParallelInitFailed', ...
                        '并行初始化失败，回退为串行。');
                    obj.use_parallel = false;
                end
            end

            % 解析常数源项（仅支持“常数/常数函数”以保持 time_step 接口不变）
            obj.source_term_value = 0.0;
            if isfield(physical_params, 'source_term')
                st = physical_params.source_term;
                if isa(st, 'function_handle')
                    obj.source_term_value = st(grid_info, physical_params);
                else
                    obj.source_term_value = st;
                end
            elseif isfield(physical_params, 'G')
                obj.source_term_value = physical_params.G;
            end

            fprintf('========== Couette流求解器 ==========\n');
            fprintf('网格: %d × %d\n', grid_info.Nx, grid_info.Ny);
            fprintf('时间步长: dt = %.6f s\n', time_params.dt);
            fprintf('空间步长: dx = %.6f, dy = %.6f\n', grid_info.dx, grid_info.dy);
            fprintf('扩散数 α_y = %.4f\n', obj.alpha_y);
            if isfield(physical_params, 'source_term') || isfield(physical_params, 'G')
                fprintf('常数源项 S = %.6f\n', obj.source_term_value);
            end
            if obj.use_parallel
                fprintf('并行: parfor 开启\n');
            end
            if obj.alpha_y <= 0.5
                fprintf('稳定性检查通过 (α_y ≤ 0.5)\n');
            end
            fprintf('====================================\n\n');
        end

        function u_new = time_step(obj, u)
            % time_step 显式欧拉推进一个时间步
            % 默认向量化；当 use_parallel=true 时，按 x 方向切片使用 parfor 并行。
            Nx = obj.grid_info.Nx;
            Ny = obj.grid_info.Ny;
            dt = obj.time_params.dt;

            u_new = u;

            if obj.use_parallel
                % parfor 在 Nx 很大或 time_step 较重时更有意义；
                % 对当前简单扩散问题可能增益有限，但保留接口便于扩展。
                %
                % 这里使用独立的 sliced 输出数组避免 parfor 变量分类问题。
                u_slice = zeros(Nx-2, Ny-2);
                parfor ii = 1:(Nx-2)
                    i = ii + 1;
                    u_slice(ii, :) = u(i, 2:Ny-1) + ...
                        obj.alpha_y * (u(i, 3:Ny) - 2*u(i, 2:Ny-1) + u(i, 1:Ny-2)) + ...
                        dt * obj.source_term_value;
                end
                u_new(2:Nx-1, 2:Ny-1) = u_slice;
            else
                % 向量化串行更新
                u_new(2:Nx-1, 2:Ny-1) = u(2:Nx-1, 2:Ny-1) + ...
                    obj.alpha_y * (u(2:Nx-1, 3:Ny) - 2*u(2:Nx-1, 2:Ny-1) + u(2:Nx-1, 1:Ny-2)) + ...
                    dt * obj.source_term_value;
            end
        end

        function is_stable = check_stability(obj)
            % check_stability 扩散格式稳定性条件
            is_stable = obj.alpha_y <= 0.5;
        end

        function a = get_alpha_y(obj)
            % 便于外部读取扩散数
            a = obj.alpha_y;
        end
    end
end
