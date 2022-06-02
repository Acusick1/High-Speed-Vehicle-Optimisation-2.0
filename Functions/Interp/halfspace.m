function [id, y_out] = halfspace(target, x, y, extrap)
%HALFSPACE vectorised search within lookup table data for input values.
%   Find closest x to target and return id, interpolate for y if required.
%   Inputs:
%       target - vector of values to search for
%       x - column vector to search for target within, or matrix containing 
%           y data
%       y - column vector/matrix with rows corresponding to x, number of 
%           rows must be equal to number of x rows
%       extrap - boolean to allow extrapolation if target lies at or beyond
%           x boundaries (default = true)
%
%   Outputs:
%       id - indices at which target value is closest to in x vector,
%           always given as the lower of the two indices where the value
%           resides
%       y_out - corresponding y values (if provided) found through 
%           interpolation 

if nargin < 3 || isempty(y)
    
    % If y data not provided, check if x is a matrix and split into x
    % (first column) and y (remaining column)
    if ~isvector(x) && size(x, 2) > 1
        
        y = x(:, 2:end);
        x = x(:, 1);
    else
        y = [];
        y_out = [];
    end
end

if nargin < 4 || isempty(extrap), extrap = true; end

% Ensure x is a column vector
x = x(:);

if ~isempty(y) && size(x, 1) ~= size(y, 1)
    
    error("Number of y rows does not match number of x rows")
end

target = target(:);
dim = size(target);

left = ones(dim);
right = left * numel(x);
previous_id = zeros(dim);
found = false(dim);

% Main halfspace search loop
% Always outputs id as last value lower than target
while any(~found)
    
    id = max(floor((left + right)/2), 1);
    
    less_than = target < x(id);
    right(less_than) = id(less_than) - 1;
    
    greater_than = target > x(id);
    left(greater_than) = id(greater_than) + 1;
    
    found(previous_id == id) = true;
    previous_id = id;
end

% If y is provided, do interpolation to find corresponding values
if ~isempty(y)
    
    % id is last value below target, therefore must inter/extrapolate
    % upwards
    lower_id = id;
    reduce = lower_id == numel(x);
    lower_id(reduce) = lower_id(reduce) - 1;
    upper_id = lower_id + 1;
    
    x_interp = [x(lower_id,:), x(upper_id, :)];
    y_interp = [y(lower_id,:), y(upper_id, :)];
    
    y_out = two_point_interp(x_interp, y_interp, target);
    
    % If extrapolation not allowed, set values to x boundaries
    if ~extrap
        
        con = any(id == [1 numel(x)], 2);
        y_out(con,:) = y(id(con),:);
    end
end