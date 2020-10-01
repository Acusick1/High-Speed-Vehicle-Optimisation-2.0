function [id, y_out] = halfspace(target, x, y)
%% Halfspace search for lookup table data
% Find closest x to target and return id, interpolate for y if required

if nargin < 3
    
    if ~isvector(x) && size(x, 2) > 1
        
        y = x(:, 2:end);
        x = x(:, 1);
    else
        y = [];
    end
end

target = target(:);
x = x(:);
dim = size(target);

L = ones(dim);
R = L * numel(x);
mPrev = zeros(dim);
found = false(dim);

% Always outputs m as last value lower than target
while any(~found)
    
    m = max(floor((L + R)/2), 1);
    
    lt = target < x(m);
    R(lt) = m(lt) - 1;
    
    gt = target > x(m);
    L(gt) = m(gt) + 1;
    
    found(mPrev == m) = true;
    mPrev = m;
end

% Save id here, as will be changed below if m = numel(x)
id = m;

if isempty(y)
    
    y_out = [];
else
    % Find corresponding y values
    maxCon = m == numel(x);
    m(maxCon) = m(maxCon) - 1;
    
    y0 = y(m, :);
    y1 = y(m + 1, :);
    x0 = x(m, :);
    x1 = x(m + 1, :);
    
    y_out = y0 + (target - x0).*((y1 - y0)./(x1 - x0));
end