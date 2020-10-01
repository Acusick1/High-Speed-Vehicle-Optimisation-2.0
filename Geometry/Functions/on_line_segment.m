function [bool, code] = on_line_segment(p1, p2, p)

[dim, ~] = size(p);
bool = true(dim, 1);
code = zeros(dim, 1);

% https://stackoverflow.com/questions/328107/how-can-you-determine-a-point-is-between-two-other-points-on-a-line-segment
u = p1 - p2;
v = p - p2;

test = sum(u .* v, 2);

before = test < 0;
after = test > sum(u.^2);

code(before) = 1;
code(after) = 2;
bool(before | after) = false;

end