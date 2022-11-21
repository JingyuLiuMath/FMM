%% Basic settings.
n = 10000;
m = 10000;

tol = 1e-6;
source_points = GeneratePoints(n, 1);
source_charges = rand(n, 1);
target_points = rand(m, 2);

%% Compute directly.
tic
potential_direct = DirectCompute(source_points, source_charges, target_points);
time_direct = toc;

%% Compute via uniformFMM
tic;
uniformFMM = uniformFMM_Tree(source_points, source_charges);
uniformFMM.FMM_Alg(tol);
potential_uniformFMM = FMM_Compute(uniformFMM, target_points);
time_uniformFMM = toc;

%% Show results.
disp("Number of source points: " + size(source_points, 1));
disp("Number of target points: " + size(target_points, 1));
disp("Tolerance: " + tol);
disp("Error of uniformFMM: " + max(abs(potential_direct - potential_uniformFMM)));
disp("Direct time: " + time_direct);
disp("uniformFMM time: " + time_uniformFMM);
