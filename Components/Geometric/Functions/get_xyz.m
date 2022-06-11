function [x, y, z] = get_xyz(mat)
%GET_XYZ split matrix into x, y, z components.
%   Input: 
%   mat - vector or matrix to extract x, y, z components from

[~, columns, matrices] = size(mat);

if matrices > 1

    x = mat(:,:,1);
    y = mat(:,:,2);
    z = mat(:,:,3);
    
elseif columns == 3
    
    x = mat(:,1);
    y = mat(:,2);
    z = mat(:,3);
    
elseif n == 3
    
    x = mat(1);
    y = mat(2);
    z = mat(3);
else
    error("Setup requried for this type of matrice input")
end

end

