%% Basic settings.
n = 10000;
m = 10000;

tol = 1e-6;
source_points = rand(n, 1);
source_charges = rand(n, 1);
target_points = rand(m, 1);

k_fun = @(x, y) 1 ./ abs(x - y);

%% Compute directly.
tic
potential_direct = DirectKernelCompute(k_fun, source_points, source_charges, target_points);
time_direct = toc;

%% Compute via uniformChebyFMM
tic;
uniformChebyFMM = uniformChebyFMM1D_Tree(source_points, source_charges);
uniformChebyFMM.FMM_Alg(k_fun, tol);
potential_uniformChebyFMM = uniformChebyFMM.FMM_Compute(target_points);
time_uniformChebyFMM = toc;

%% Show results.
disp("Number of source points: " + size(source_points, 1));
disp("Number of target points: " + size(target_points, 1));
disp("Tolerance: " + tol);
disp("Error of uniformChebyFMM: " + max(abs(potential_direct - potential_uniformChebyFMM)));
disp("Relative error of uniformChebyFMM: " + max(abs(potential_direct - potential_uniformChebyFMM) ./ abs(potential_direct)));
disp("Direct time: " + time_direct);
disp("uniformChebyFMM time: " + time_uniformChebyFMM);
