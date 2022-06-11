function output = crossmat(mat1, mat2)
%CROSSMAT Calculates elemental cross product x, y, z vectors and matrices.
%   Input:
%   mat1, mat2 - x, y, z vectors or matrices specified as 1x3 vector, 
%       nx3 matrices or n x m x 3. 
%
%   Matrices can have different dimensions e.g. mat1 = 1x3 vector crossed 
%   5x5x3 matrix will return a 5x5x3 matrix where each matrix element has 
%   been crossed by the vector.

[x1, y1, z1] = get_xyz(mat1);
[x2, y2, z2] = get_xyz(mat2);

% Dimension to concatenate cross product values on, defined by highest 
% dimension input
highest_dim = max(numel(size(mat1)), numel(size(mat2)));

x_out = y1 .* z2 - z1 .* y2;
y_out = z1 .* x2 - x1 .* z2;
z_out = x1 .* y2 - y1 .* x2;

output = cat(highest_dim, x_out, y_out, z_out);

end