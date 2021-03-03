function yq = twoPointInterp(x, y ,xq)
%TWOPOINTINTERP linearly inter/extrapolates between two points
%   Vectorised to handle arbitrary number of two point sets
%   x, y, specified as 2 column vectors
%   xq specified as single column vector

yq = y(:,1) + (xq - x(:,1)) .* (y(:,2) - y(:,1))./(x(:,2) - x(:,1));