function [bool, xcode, ycode] = is_within_parallelogram(x, y, xt, yt)
%IS_WITHIN_POLYGON test if points lie within polygon bound by x and y,
%vectorised to test series of parallelograms at once
%   Inputs:
%   x - parallel component, must be monotonically increasing column vector 
%       which contains both upper and lower y bounds
%   y - matrix of points at given x locations, with first and last column
%       giving boundaries 

[id, ybound] = halfspace(xt, x, y(:,[1 end]), false);

xcode = zeros(size(xt));
ycode = xcode;

xcode(id == 1) = 1;
xcode(id == numel(x)) = 2;

ycode(yt < ybound(:,end)) = 3;
ycode(yt > ybound(:,1)) = 4;

% If a point lies outside xbound, ybound will not be accurate (or relevant)
% ycode(yt < ybound(:,2) & ~xcode) = 3;
% ycode(yt > ybound(:,1) & ~xcode) = 4;

bool = ~xcode & ~ycode;