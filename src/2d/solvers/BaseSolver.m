classdef (Abstract) BaseSolver < handle
    % BaseSolver 抽象基类
    %
    % 统一求解器接口：
    %   - time_step: 推进一个时间步
    %   - check_stability: 稳定性检查
    % 提供通用 run_simulation 时间推进循环。

    properties (Access = protected)
        grid_info          % 网格信息结构体: Nx, Ny, dx, dy, x, y, Lx, Ly
        physical_params    % 物理参数结构体: U, nu, ...
        time_params        % 时间参数结构体: dt, t_end, n_snapshots
    end

    methods
        function obj = BaseSolver(grid_info, physical_params, time_params)
            obj.grid_info = grid_info;
            obj.physical_params = physical_params;
            obj.time_params = time_params;
        end
    end

    methods (Abstract)
        u_new = time_step(obj, u)
        is_stable = check_stability(obj)
    end

    methods
        function [u_history, time_history] = run_simulation(obj, u0, bc_manager, n_snapshots)
            % run_simulation 执行完整时间推进并按快照输出
            %
            % 输入:
            %   u0          - 初始速度场 [Nx×Ny]
            %   bc_manager  - CompositeBoundaryCondition 对象
            %   n_snapshots - 输出快照数量（可选，默认使用 time_params 中值）
            %
            % 输出:
            %   u_history   - 速度场历史记录（cell数组）
            %   time_history- 时间历史记录（向量）

            if nargin < 4 || isempty(n_snapshots)
                n_snapshots = obj.time_params.n_snapshots;
            end

            dt = obj.time_params.dt;
            t_end = obj.time_params.t_end;

            % 稳定性检查
            if ~obj.check_stability()
                error('稳定性条件违反：求解器不稳定，请检查 dt/网格/物理参数。');
            end

            output_interval = t_end / (n_snapshots - 1);
            u_history = cell(1, n_snapshots);
            time_history = zeros(1, n_snapshots);

            % 初始边界条件
            u = bc_manager.apply_all(u0, obj.grid_info, obj.physical_params);

            snapshot_idx = 1;
            u_history{snapshot_idx} = u;
            time_history(snapshot_idx) = 0.0;
            next_output_time = output_interval;

            current_time = 0.0;
            step = 0;

            fprintf('开始时间推进...\n');
            fprintf('------------------------------------\n');
            tic;

            while current_time < t_end - 1e-12
                step = step + 1;

                % 一个时间步
                u_new = obj.time_step(u);

                % 应用边界条件
                u = bc_manager.apply_all(u_new, obj.grid_info, obj.physical_params);

                current_time = current_time + dt;

                % 保存快照
                if current_time >= next_output_time - 1e-10 && snapshot_idx < n_snapshots
                    snapshot_idx = snapshot_idx + 1;
                    u_history{snapshot_idx} = u;
                    time_history(snapshot_idx) = current_time;
                    next_output_time = snapshot_idx * output_interval;
                end

                % 进度显示
                if mod(step, 10000) == 0
                    fprintf('进度: %.1f%% (t=%.3fs, 步数=%d)\n', ...
                            current_time / t_end * 100, current_time, step);
                end
            end

            % 确保保存最终状态
            if snapshot_idx < n_snapshots
                snapshot_idx = snapshot_idx + 1;
                u_history{snapshot_idx} = u;
                time_history(snapshot_idx) = current_time;
            end

            elapsed_time = toc;
            fprintf('------------------------------------\n');
            fprintf('时间推进完成！耗时: %.2f 秒\n', elapsed_time);
            fprintf('实际计算步数: %d\n', step);
            fprintf('平均速度: %.0f 步/秒\n\n', step / elapsed_time);
        end
    end
end

