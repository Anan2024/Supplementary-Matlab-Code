function [theta, A] = computeAEMMedian(X)
      % 输入: 
    % X - 样本数据矩阵 (m x p)
    
    % 初始化 theta 和 A
    [m, p] = size(X);
    theta = mean(X)'; % 初始值为均值
    A = eye(p);
    tol = 1e-6;
    max_iter = 1000;
    
    for iter = 1:max_iter
        theta_old = theta;
        A_old = A;
        
        % 更新 theta
        diffs = X' - repmat(theta_old, 1, m);
        norms = sqrt(sum((A_old * diffs).^2, 1));
        weights = 1 ./ norms;
        theta = sum(diffs .* repmat(weights, p, 1), 2) / sum(weights);
        
        % 更新 A
        diffs = X' - repmat(theta, 1, m);
        S = (diffs * diffs') / m;
        A = (S^(-0.5)) / mean(norms);
        
        % 检查收敛性
        if norm(theta - theta_old) < tol && norm(A - A_old, 'fro') < tol
            break;
        end
    end
    
%     if iter == max_iter
%         warning('迭代达到最大次数，未收敛');
%     end
end