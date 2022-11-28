function points = ChebyPoints2D(left, right, n)
% Chebypoints Compute n^d Chebyshev points of the box [left, right] where d is the dimension.

% Jingyu Liu, November 27, 2022.

% Chebyshev points on [-1, 1].
theta = 1 : 2 : 2 * n - 1;
theta = theta * pi / 2 / n;
points = cos(theta);

% Linear transform to [left, right]. 
points1 = (points' + 1) * (right(1) - left(1)) / 2 + left(1);
points2 = (points' + 1) * (right(2) - left(2)) / 2 + left(2);

points = [kron(points1, ones(n, 1)), kron(ones(n, 1), points2)];

end