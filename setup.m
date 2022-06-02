function setup()
%SETUP sets up base and results path for current session, and adds
%necessary directories to MATLAB path.

base_path = pwd;
results_path = fullfile(base_path, "Results");

addpath(genpath(fullfile(base_path, "Aerodynamics")))
addpath(genpath(fullfile(base_path, "External")))
addpath(genpath(fullfile(base_path, "Functions")))
addpath(genpath(fullfile(base_path, "Geometry")))
addpath(genpath(fullfile(base_path, "Optimisation")))

setenv("BASE_PATH", base_path)
setenv("RESULTS_PATH", results_path)