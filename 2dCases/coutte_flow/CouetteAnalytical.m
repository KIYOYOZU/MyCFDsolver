classdef CouetteAnalytical < handle
    % CouetteAnalytical Couette/Poiseuille-Couette 启动问题解析解（算例专用）
    %
    % 控制方程:
    %   ∂u/∂t = ν · ∂²u/∂y² + G
    % 其中 G = -(1/ρ)·dp/dx 为常数压力梯度等效体力。
    %
    % 理论解分解:
    %   u = u_s(y) + Σ B_n sin(nπy/h) exp(-n²π²νt/h²)
    %   u_s(y) = U·(y/h) + (G/(2ν))·(h y - y²)

    properties
        U
        nu
        h
        G = 0.0
        n_terms = 100
    end

    methods
        function obj = CouetteAnalytical(U, nu, h, varargin)
            obj.U = U;
            obj.nu = nu;
            obj.h = h;
            % 兼容两种调用：
            %   CouetteAnalytical(U,nu,h)                     -> G=0, n_terms=100
            %   CouetteAnalytical(U,nu,h,n_terms)            -> G=0, n_terms=指定
            %   CouetteAnalytical(U,nu,h,G,n_terms)          -> 压力驱动 + 指定项数
            if numel(varargin) == 1
                v = varargin{1};
                if isscalar(v) && v == round(v) && v >= 1
                    obj.n_terms = v;
                else
                    obj.G = v;
                end
            elseif numel(varargin) >= 2
                obj.G = varargin{1};
                obj.n_terms = varargin{2};
            end
        end

        function u_ana = compute_profile(obj, y, t)
            % compute_profile 返回给定时间的速度剖面（列向量）
            y_col = y(:);

            h = obj.h;
            nu = obj.nu;
            U = obj.U;
            G = obj.G;

            % 稳态分量（Poiseuille-Couette）
            u_steady = U * (y_col / h) + (G / (2 * nu)) * (h * y_col - y_col.^2);

            % 瞬态分量 Fourier 级数（齐次边界）
            n = (1:obj.n_terms)';
            temporal = exp(-n.^2 * pi^2 * nu * t / h^2);

            % 系数 B_n
            B_couette = (2 * U ./ (n * pi)) .* (-1).^n;
            B_poiseuille = -(2 * G * h^2 / nu) .* (1 - (-1).^n) ./ (n.^3 * pi^3);
            B = B_couette + B_poiseuille;

            spatial = sin(n * pi * y_col' / h); % [n_terms×Ny]
            u_transient = sum((B .* temporal) .* spatial, 1)';

            u_ana = u_steady + u_transient;
        end
    end
end
