function inc_array(x, dx)

% x: array to be increasing
% dx: minimum spacing required between array

if nargin < 2 || isempty(dx), dx = 0; end

for i = 2:numel(x)
    
    if x(i) <= x(i-1), x(i) = x(i-1) + dx; end
end