classdef DirichletBC < BoundaryCondition
    % DirichletBC 固定值边界条件
    %
    % 支持任意方向/位置:
    %   DirichletBC('y','min',0)  -> y=0 边界 u=0
    %   DirichletBC('y','max',U)  -> y=Ly 边界 u=U

    properties
        value   % 数值或函数句柄
    end

    methods
        function obj = DirichletBC(direction, position, value)
            obj@BoundaryCondition(direction, position);
            obj.value = value;
        end

        function u = apply(obj, u, grid_info, physical_params)
            Nx = grid_info.Nx;
            Ny = grid_info.Ny;

            val = obj.value;
            if isa(val, 'function_handle')
                val = val(grid_info, physical_params);
            end

            switch obj.direction
                case 'y'
                    if strcmp(obj.position, 'min')
                        u(:, 1) = val;
                    elseif strcmp(obj.position, 'max')
                        u(:, Ny) = val;
                    else
                        error('DirichletBC:position', '未知边界位置: %s', obj.position);
                    end
                case 'x'
                    if strcmp(obj.position, 'min')
                        u(1, :) = val;
                    elseif strcmp(obj.position, 'max')
                        u(Nx, :) = val;
                    else
                        error('DirichletBC:position', '未知边界位置: %s', obj.position);
                    end
                otherwise
                    error('DirichletBC:direction', '未知边界方向: %s', obj.direction);
            end
        end
    end
end

