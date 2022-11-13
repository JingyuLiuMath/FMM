%% Basic settings.
n = 200;
m = 200;

tol = 1e-3;

source_points = rand(n, 2);
source_charges = rand(n, 1);
target_points = rand(m , 2);

%% Compute directly.
[potential_direct, time_direct] = DirectlyCompute(source_points, source_charges, target_points);

%% Compute via FMM
tic;
FMM = uniformFMM_Tree(source_points, source_charges);
FMM = FMM_Alg(FMM, tol);
potential_FMM = FMM_Compute(FMM, target_points);
time_FMM = toc;

%% Results.
disp("error: " + string(max(abs(potential_direct - potential_FMM))));
disp("Direct time: " + time_direct);
disp("FMM time: " + time_FMM);
