function points = GeneratePoints(n, type)
% GeneratePoints Generate points in [0, 1]^2.

% type 1: uniform points in [0, 1]^2.
% type 2: unform points in [0.1, 0.9]^2.
% type 3: almost half uniform in [0, 0.5]^2, almost half uniform in [0.5, 1]^2, several in the rest.

% Jingyu Liu, November 16, 2022.

arguments
    n(1, 1);
    type(1, 1) = 0;
end

switch type
    case 0    
        points = rand(n, 2);
    case 1
        points = rand(n, 2) * 0.8 + 0.1;
    case 2
        n1 = floor(n / 1000);
        n4 = n1;
        n2 = floor((n - n1 - n4) / 2);
        n3 = n - n1 - n2 - n4;
        points1 = [0.5 * rand(n1, 1), 0.5 * 0.5 * rand(n1, 1) + 0.5];
        points2 = 0.5 * rand(n2, 2) + 0.5;
        points3 = 0.5 * rand(n3, 2);
        points4 = [0.5 * 0.5 * rand(n4, 1) + 0.5, 0.5 * rand(n4, 1)];
        points = [points1; points2; points3; points4];
end

end