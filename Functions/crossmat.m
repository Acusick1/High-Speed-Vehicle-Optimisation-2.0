function vectors = crossmat(mat1, mat2, dim)
%% Calculates elemental cross product x,y,z vectors of up to 3D matrices

if nargin < 3
    
    [x1, y1, z1] = format(mat1);
    [x2, y2, z2] = format(mat2);

    vectors(:,:,1) = y1 .* z2 - z1 .* y2;
    vectors(:,:,2) = z1 .* x2 - x1 .* z2;
    vectors(:,:,3) = x1 .* y2 - y1 .* x2;
    
    %% TODO: Something more robust to ensure output if matched
    if ~isequal(size(vectors), size(mat1))
        
        vectors = squeeze(vectors);
    end
    
    if ~isequal(size(vectors), size(mat1))
        
        vectors = vectors';
    end

elseif dim == 1
    
    vectors = one(mat1, mat2);
    
elseif dim == 2
    
    vectors = two(mat1, mat2);
    
elseif dim == 3
    
    vectors = three(mat1, mat2);
end

end

function vectors = one(mat1, mat2)
    
    x1 = mat1(1);
    y1 = mat1(2);
    z1 = mat1(3);
    
    x2 = mat2(1);
    y2 = mat2(2);
    z2 = mat2(3);
    
    vectors(3) = x1 * y2 - y1 * x2;
    vectors(2) = z1 * x2 - x1 * z2;
    vectors(1) = y1 * z2 - z1 * y2;
end

function vectors = two(mat1, mat2)
    
    x1 = mat1(:,1);
    y1 = mat1(:,2);
    z1 = mat1(:,3);
    
    x2 = mat2(:,1);
    y2 = mat2(:,2);
    z2 = mat2(:,3);
    
    vectors(:,3) = x1 * y2 - y1 * x2;
    vectors(:,2) = z1 * x2 - x1 * z2;
    vectors(:,1) = y1 * z2 - z1 * y2;

end

function vectors = three(mat1, mat2)
    
    x1 = mat1(:,:,1);
    y1 = mat1(:,:,2);
    z1 = mat1(:,:,3);
    
    x2 = mat2(:,:,1);
    y2 = mat2(:,:,2);
    z2 = mat2(:,:,3);
    
    vectors(:,:,3) = x1 * y2 - y1 * x2;
    vectors(:,:,2) = z1 * x2 - x1 * z2;
    vectors(:,:,1) = y1 * z2 - z1 * y2;
end

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

end