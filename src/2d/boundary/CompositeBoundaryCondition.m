classdef CompositeBoundaryCondition < handle
    % CompositeBoundaryCondition 边界条件组合管理器
    %
    % 负责按顺序应用多个 BoundaryCondition。

    properties (Access = private)
        bc_list   % cell数组，存储边界条件对象
    end

    methods
        function obj = CompositeBoundaryCondition()
            obj.bc_list = {};
        end

        function add_bc(obj, bc)
            if ~isa(bc, 'BoundaryCondition')
                error('CompositeBoundaryCondition:add_bc', ...
                      '输入必须是 BoundaryCondition 子类对象。');
            end
            obj.bc_list{end+1} = bc;
        end

        function u = apply_all(obj, u, grid_info, physical_params)
            for k = 1:numel(obj.bc_list)
                u = obj.bc_list{k}.apply(u, grid_info, physical_params);
            end
        end

        function n = count(obj)
            n = numel(obj.bc_list);
        end
    end
end

