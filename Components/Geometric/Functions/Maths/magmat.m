function magnitudes = magmat(mat, p)
%MAGMAT Wrapper for vecnorm function when finding magnitude of x, y, z
%coordinates.
%   Inputs:
%   mat - x, y, z matrix specified as nx3 matrix or n x m x 3 matrix
%   p - p norm, see norm function for more details

dims = size(mat);

% If p not specified, assume Euclidean norm
if nargin < 2 || isempty(p), p = 2; end

if dims(end) ~= 3
    
    error("Incompatible format, final dimension should equal 3 for an x, y, z matrix")
end

% Specifies dimension to calculate norm based on total dimensions of input
% matrix
magnitudes = vecnorm(mat, p, numel(dims));