function points = GeneratePoints(n, type)
% GeneratePoints Generate points in [0, 1]^2.

% type 1: uniform distributed points in [0, 1]^2.
% type 2: uniformly distributed points in a circle of radius 0.5 about the center of [0, 1]^2.

% Jingyu Liu, December 30, 2022.

arguments
    n(1, 1);
    type(1, 1) = 1;
end

switch type
    case 1
        points = rand(n, 2);
    case 2
        angle = rand(n, 1) * 2 * pi;
        r = 0.5;
        points = zeros(n, 1);
        points(:, 1) = r .* cos(angle);
        points(:, 2) = r .* sin(angle);
        points = points + 0.5;
end

end