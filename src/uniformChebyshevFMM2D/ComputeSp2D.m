function weight = ComputeSp2D(left, right, p, points1, points2)
% ComputeSp Compute Sp(point1, point2).

% Jingyu Liu, November 27, 2022.

d = size(points1, 2);  % Dimension.
m = size(points1, 1);
n = size(points2, 1);

% Convert to [-1, 1].
points1 = (points1 - left) * 2 ./ (right - left) - 1;
points2 = (points2 - left) * 2 ./ (right - left) - 1;

% Compute Sp.
weight_tensor = ones(m, n, d);
for t = 1 : d
    for k = 1 : p - 1
        weight_tensor(:, :, t) =  weight_tensor(:, :, t) + ...
            2 * kron(cos(k * acos(points1(:, t))), cos(k * acos(points2(:, t))'));
    end
end

weight_tensor = weight_tensor / p;

weight = weight_tensor(:, :, 1);
for t = 2 : d
    weight = weight .* weight_tensor(:, :, t);
end

end