function [bool, xcode, ycode] = is_within_polygon(x, y, xt, yt)

[id, ybound] = halfspace(xt, x, y(:,[1 end]));

xcode = zeros(size(xt));
ycode = xcode;

xcode(id == 1) = 1;
xcode(id == numel(x)) = 2;

ycode(yt < ybound(:,2)) = 3;
ycode(yt > ybound(:,1)) = 4;

bool = ~xcode & ~ycode;