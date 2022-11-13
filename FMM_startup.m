function FMM_startup()
% FMM_startup Startup function for FMM.

file_path = mfilename('fullpath');
tmp = strfind(file_path,'FMM_startup');
file_path = file_path(1:(tmp(end)-1));
addpath(genpath([file_path 'src']));
addpath(genpath([file_path 'test']));

end