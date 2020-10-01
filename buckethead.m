addpath(genpath('Functions'))
addpath(genpath('Optimisation'))
addpath(genpath('Geometry'))
addpath(genpath('Aerodynamics'))

% read number of procs
fid = fopen('matlab-parallel.sge');
% set linenum to the desired line number that you want to import
linenum = 15;
C = textscan(fid,'%s',4,'headerlines',linenum-1);
C = C{:};

nProc = str2double(C{4});
parpool(nProc)
fclose(fid);

foilopt