function points = GeneratePoints(n, type)
% GeneratePoints Generate points in [0, 1]^2.

% type 1: uniform distributed points in [0, 1]^2.
% type 2: uniformly distributed points in a circle of radius 0.5 about the center of [0, 1]^2.
% type 3: A fifth are uniformly assigned in [0, 1]^2, two fifths are randomly distributed about the 
% center of the square in a circle of radius 0.003, the rest are uniformly inside a circle of radius 
% 0.5.

% Jingyu Liu, November 17, 2022.

arguments
    n(1, 1);
    type(1, 1) = 1;
end

switch type
    case 1
        points = rand(n, 2);
    case 2
        angle = rand(n, 1) * 2 * pi;
        r = sqrt(rand(n, 1)) * 0.5;
        points = zeros(n, 1);
        points(:, 1) = r .* cos(angle);
        points(:, 2) = r .* sin(angle);
        points = points + 0.5;
    case 3
        n1 = ceil(n / 5);
        points1 = rand(n1, 2);
        n2 = ceil(n * 2 / 5);
        angle = rand(n2, 1) * 2 * pi;
        r = sqrt(rand(n2, 1)) * 0.003;
        points2 = zeros(n2, 1);
        points2(:, 1) = r .* cos(angle);
        points2(:, 2) = r .* sin(angle);
        points2 = points2 + 0.5;
        n3 = max(n - n1 - n2, 0);
        angle = rand(n3, 1) * 2 * pi;
        r = sqrt(rand(n3, 1)) * 0.5;
        points3 = zeros(n3, 2);
        points3(:, 1) = r .* cos(angle);
        points3(:, 2) = r .* sin(angle);
        points3 = points3 + 0.5;
        points = [points1; points2; points3];
end

end