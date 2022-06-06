function [nx, ny] = normal2D(x, y, dim, invert)
%NORMAL2D Calculate two dimensional unit normal vector
%   Input:
%   x, y - x and y coordinates of a single line
% 	dim - which dimension normal is to be calculated (default = 1)
%   invert - compute inverted normal

if nargin < 3 || isempty(dim)
    
    dim = 1;
    
    if all(size(x) > 1)
    
        error('Matrix input without dimension in which to calulcate normal')
    end
end

if nargin < 4 || isempty(invert), invert = false; end 

% Calculate vector components of each line segment
diffx = diff(x, [], dim);
diffy = diff(y, [], dim);
    
% Calculate normal components
if invert
    
    nx = diffy;
    ny = -diffx;
else
    nx = -diffy;
    ny = diffx;
end

mag = (nx.^2 + ny.^2).^0.5;

nx = nx./mag;
ny = ny./mag;

%% Proof
% xc = (x(1:end-1) + x(2:end))/2;
% yc = (y(1:end-1) + y(2:end))/2;
% 
% figure
% hold on
% plot(x,y)
% plot(xc + nx,yc + ny);