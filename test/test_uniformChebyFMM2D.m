%% Basic settings.
n = 20000;
m = 20000;

tol = 1e-6;
source_points = rand(n, 2);
source_charges = rand(n, 1);
target_points = rand(m, 2);

k_fun = @(x, y) log(1 ./ pdist2(x, y));

uniformFMM_flag = 1;  % Whether need to be compared with uniform FMM. 

%% Compute directly.
tic
potential_direct = DirectKernelCompute2D(k_fun, source_points, source_charges, target_points);
time_direct = toc;

%% Compute via uniformFMM
if uniformFMM_flag == 1
    tic;
    uniformFMM = uniformFMM_Tree(source_points, source_charges);
    uniformFMM.FMM_Alg(tol);
    potential_uniformFMM = uniformFMM.FMM_Compute(target_points);
    time_uniformFMM = toc;
end

%% Compute via uniformChebyFMM
tic;
uniformChebyFMM = uniformChebyFMM2D_Tree(source_points, source_charges);
uniformChebyFMM.FMM_Alg(k_fun, tol);
potential_uniformChebyFMM = uniformChebyFMM.FMM_Compute(target_points);
time_uniformChebyFMM = toc;
 
%% Show results.
disp("Number of source points: " + size(source_points, 1));
disp("Number of target points: " + size(target_points, 1));
disp("Tolerance: " + tol);
if uniformFMM_flag == 1
    disp("Error of uniformFMM: " + max(abs(potential_direct - potential_uniformFMM)));
    disp("Relative error of uniformFMM: " + max(abs(potential_direct - potential_uniformFMM) ./ abs(potential_direct)));
end
disp("Error of uniformChebyFMM: " + max(abs(potential_direct - potential_uniformChebyFMM)));
disp("Relative error of uniformChebyFMM: " + max(abs(potential_direct - potential_uniformChebyFMM) ./ abs(potential_direct)));
disp("Direct time: " + time_direct);
if uniformFMM_flag == 1
    disp("uniformFMM time: " + time_uniformFMM);
end
disp("uniformChebyFMM time: " + time_uniformChebyFMM);
