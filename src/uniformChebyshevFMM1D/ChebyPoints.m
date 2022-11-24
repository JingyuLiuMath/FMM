function points = ChebyPoints(a, b, n)
% Chebypoints Compute n Chebyshev points of the interval [a, b].

% Jingyu Liu, November 22, 2022.

% Chebyshev points on [-1, 1].
theta = 1 : 2 : 2 * n - 1;
theta = theta * pi / 2 / n;
points = cos(theta);

% Linear transform to [a, b].
points = (b - a) / 2 * (points + 1) + a;
points = points';

end