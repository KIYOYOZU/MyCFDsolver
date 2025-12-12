classdef Mesh2D < handle
    % Mesh2D 2D 结构化均匀网格模块（通用）
    %
    % 负责生成并保存常用网格信息：
    %   Nx, Ny, Lx, Ly, dx, dy, x, y
    %
    % 使用:
    %   mesh = Mesh2D(Nx, Ny, Lx, Ly);
    %   grid_info = mesh.to_struct();

    properties (SetAccess = private)
        Nx
        Ny
        Lx
        Ly
        dx
        dy
        x
        y
    end

    methods
        function obj = Mesh2D(Nx, Ny, Lx, Ly)
            obj.Nx = Nx;
            obj.Ny = Ny;
            obj.Lx = Lx;
            obj.Ly = Ly;

            obj.dx = Lx / (Nx - 1);
            obj.dy = Ly / (Ny - 1);
            obj.x = linspace(0, Lx, Nx);
            obj.y = linspace(0, Ly, Ny);
        end

        function s = to_struct(obj)
            % to_struct 返回与求解器兼容的 grid_info 结构体
            s = struct('Nx', obj.Nx, 'Ny', obj.Ny, ...
                       'Lx', obj.Lx, 'Ly', obj.Ly, ...
                       'dx', obj.dx, 'dy', obj.dy, ...
                       'x', obj.x, 'y', obj.y);
        end
    end
end

