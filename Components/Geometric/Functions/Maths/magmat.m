function magnitudes = magmat(mat, p)
%MAGMAT Wrapper for vecnorm function when finding magnitude of x, y, z
%coordinates.
%   Inputs:
%   mat - x, y or x, y, z matrix where final dimension signifies is where
%       magnitude is to be found.
%   p - p norm, see norm function for more details

dims = size(mat);

% If p not specified, assume Euclidean norm
if nargin < 2 || isempty(p), p = 2; end

% Allowing x, y and x, y, z inputs
if ~any(dims(end) == [2 3])
    
    error("Incompatible format, final dimension should equal 3 for an x, y, z matrix")
end

% Specifies dimension to calculate norm based on total dimensions of input
% matrix
magnitudes = vecnorm(mat, p, numel(dims));