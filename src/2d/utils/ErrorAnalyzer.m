classdef ErrorAnalyzer
    % ErrorAnalyzer 误差分析工具
    %
    % 提供常用误差指标的静态方法。

    methods (Static)
        function L2_err = compute_L2(u_num, u_ref)
            % compute_L2 计算相对 L2 范数误差
            %
            % 输入:
            %   u_num - 数值解
            %   u_ref - 参考解（同维度）
            %
            % 输出:
            %   L2_err - ||u_num-u_ref||_2 / ||u_ref||_2

            error_sq = sum((u_num - u_ref).^2, 'all');
            ref_sq = sum(u_ref.^2, 'all');
            L2_err = sqrt(error_sq / ref_sq);
        end
    end
end

