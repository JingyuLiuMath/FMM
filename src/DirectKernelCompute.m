function potential = DirectKernelCompute(k_fun, source_points, source_charges, target_points)
% DirectCompute Compute potential from source points to target points via kernel function directly.

% Jingyu Liu, November 24, 2022.

arguments
    k_fun;
    source_points(:, :);
    source_charges(:, 1);
    target_points(:, :);
end

% n = size(source_points, 1);
% m = size(target_points, 1);
% potential = zeros(m, 1);
% for i = 1 : m
%     for j = 1 : n
%         potential(i) = potential(i) + source_charges(j) * k_fun(target_points(i, :), source_points(j, :));
%     end
% end

potential = k_fun(target_points, source_points') * source_charges;

end