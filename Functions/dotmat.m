function vectors = dotmat(mat1, mat2)
%% EXPIRED: USE DOT(MAT1, MAT2) OR SUM(MAT1 .* MAT2, 2) INSTEAD

% error("EXPIRED")

%% Calculates elemental cross product for various combinations of x,y,z matrices

[x1, y1, z1] = format(mat1);
[x2, y2, z2] = format(mat2);

vectors = x1 .* x2 + y1 .* y2 + z1 .* z2;

function [x, y, z, form] = format(mat)

[~,cols,threeD] = size(mat);

if threeD > 1
    
    x = mat(:,:,1);
    y = mat(:,:,2);
    z = mat(:,:,3);
    form = 3;
    
elseif cols == 3
    
    % Assume matrices are set up as [x y z]
    x = mat(:,1);
    y = mat(:,2);
    z = mat(:,3);
    form = 2;
    
elseif numel(mat) == 3
    
    x = mat(1);
    y = mat(2);
    z = mat(3);
    form = 1;
else
    error('Setup requried for this type of matrice input')
end