%% Basic settings.
n = 20000;
m = 20000;

min_points = 512;
tol = 1e-6;

source_points = GeneratePoints(n, 2);
source_charges = rand(n, 1);
target_points = rand(m, 2);

uniformFMM_flag = 0;  % Whether need to be compared with uniform FMM. 

%% Compute directly.
tic
potential_direct = DirectCompute(source_points, source_charges, target_points);
time_direct = toc;

%% Compute via uniformFMM
if uniformFMM_flag == 1
    tic;
    uniformFMM = uniformFMM_Tree(source_points, source_charges);
    uniformFMM.FMM_Alg(tol);
    potential_uniformFMM = FMM_Compute(uniformFMM, target_points);
    time_uniformFMM = toc;
end

%% Compute via aFMM
tic;
aFMM = aFMM_Tree(source_points, source_charges, min_points);
aFMM.FMM_Alg(tol);
potential_aFMM = FMM_Compute(aFMM, target_points);
time_aFMM = toc;

%% Show results.
disp("Number of source points: " + size(source_points, 1));
disp("Number of target points: " + size(target_points, 1));
disp("Minimum points: " + min_points)
disp("Tolerance: " + tol);
if uniformFMM_flag == 1
    disp("Error of uniformFMM: " + max(abs(potential_direct - potential_uniformFMM)));
end
disp("Error of aFMM: " + max(abs(potential_direct - potential_aFMM)));
disp("Direct time: " + time_direct);
if uniformFMM_flag == 1
    disp("uniformFMM time: " + time_uniformFMM);
end
disp("aFMM time: " + time_aFMM);
