function output = dotmat(mat1, mat2)
%DOTMAT Calculates dot product in vectorised fashion for vectors and
%matrices.
%   Input:
%   mat1, mat2 - x, y, z vectors or matrices specified as 1x3 vector, 
%       nx3 matrices or n x m x 3. 
%
%   Matrices can have different dimensions e.g. mat1 = 1x3 vector and 
%   5x5x3 matrix will return a 5x5 matrix where each matrix has been
%   multiplied by the vector counterpart

[x1, y1, z1] = get_xyz(mat1);
[x2, y2, z2] = get_xyz(mat2);

output = x1 .* x2 + y1 .* y2 + z1 .* z2;