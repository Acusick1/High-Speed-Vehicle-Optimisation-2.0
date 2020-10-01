function I = new_line_plane_intersection(u, N, n, M)

% Plane offset parameter
d = -sum(n .* M, 2);

% Parametric line parameter t
t = -(d + sum(n .* N, 2)) ./ sum(n .* u, 2);

% Intersection coordinates
I = N + u .* t;