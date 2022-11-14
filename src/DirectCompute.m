function potential = DirectCompute(source_points, source_charges, target_points)
% DirectCompute Compute potential from source points to target points
% directly.

% Jingyu Liu, November 14, 2022.

arguments
    source_points(:, 2);
    source_charges(:, 1);
    target_points(:, 2);
end

n = size(source_points, 1);
m = size(target_points, 1);
source_points_complex = complex(source_points(:, 1), source_points(:, 2));
target_points_complex = complex(target_points(:, 1), target_points(:, 2));
potential = zeros(m, 1);
for i = 1 : m
    for j = 1 : n
        potential(i) = potential(i) - source_charges(j) * log(abs(target_points_complex(i) - source_points_complex(j)));
    end
end

end