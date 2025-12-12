classdef (Abstract) BoundaryCondition < handle
    % BoundaryCondition 边界条件抽象基类
    %
    % 子类需实现:
    %   u = apply(obj, u, grid_info, physical_params)

    properties (Access = protected)
        direction   % 'x' 或 'y'
        position    % 'min' 或 'max'（周期边界可为空）
    end

    methods
        function obj = BoundaryCondition(direction, position)
            if nargin < 1
                direction = '';
            end
            if nargin < 2
                position = '';
            end
            obj.direction = lower(direction);
            obj.position = lower(position);
        end
    end

    methods (Abstract)
        u = apply(obj, u, grid_info, physical_params)
    end
end

