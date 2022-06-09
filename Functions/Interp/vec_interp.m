function yq = vec_interp(x0, x1, y0, y1 ,xq)
%VEC_INTERP linearly inter/extrapolates between two point sets
%   Vectorised to handle arbitrary and dissimilar point sets

yq = y0 + (xq - x0) .* (y1 - y0)./(x1 - x0);