classdef PeriodicBC < BoundaryCondition
    % PeriodicBC 周期性边界条件
    %
    % 使用方式:
    %   PeriodicBC('x') -> x方向周期
    %   PeriodicBC('y') -> y方向周期

    methods
        function obj = PeriodicBC(direction)
            obj@BoundaryCondition(direction, '');
        end

        function u = apply(obj, u, grid_info, ~)
            Nx = grid_info.Nx;
            Ny = grid_info.Ny;

            switch obj.direction
                case 'x'
                    % 左右边界周期
                    u(1, :) = u(Nx-1, :);
                    u(Nx, :) = u(2, :);
                case 'y'
                    % 上下边界周期
                    u(:, 1) = u(:, Ny-1);
                    u(:, Ny) = u(:, 2);
                otherwise
                    error('PeriodicBC:direction', '未知周期方向: %s', obj.direction);
            end
        end
    end
end

