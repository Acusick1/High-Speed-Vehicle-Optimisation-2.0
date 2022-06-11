function I = line_plane_intersection_simple(u, N, n, M)
%LINE_PLANE_INTERSECTION_SIMPLE simplified version of
%line_plan_intersection, see documentation there for usage info

% Plane offset parameter
d = -sum(n .* M, 2);

% Parametric line parameter t
t = -(d + sum(n .* N, 2)) ./ sum(n .* u, 2);

% Intersection coordinates
I = N + u .* t;