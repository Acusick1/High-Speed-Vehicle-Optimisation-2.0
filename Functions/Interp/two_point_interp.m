function yq = two_point_interp(x, y ,xq)
%TWO_POINT_INTERP linearly inter/extrapolates between two points
%   Vectorised to handle arbitrary number of two point sets
%   x, y, specified as 2 column vectors
%   xq specified as column vector

if any([size(y, 2) size(x, 2)] ~= 2)
    
    error("x and y must be specified as matrices with exactly two columns")
end

yq = y(:,1) + (xq - x(:,1)) .* (y(:,2) - y(:,1))./(x(:,2) - x(:,1));