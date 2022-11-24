function weight = ComputeSp(a, b, p, points1, points2)
% ComputeSp Compute Sp(point1, point2).

% Jingyu Liu, November 24, 2022.

m = length(points1);
n = length(points2);

% Convert to [-1, 1].
points1 = (points1 - a) * 2 / (b - a) - 1;
points2 = (points2 - a) * 2 / (b - a) - 1;

% Compute Sp.
weight = ones(m, n);
for k = 1 : p - 1
    weight = weight + 2 * kron(cos(k * acos(points1)), cos(k * acos(points2))');
end

weight = weight / p;

end