function within = is_within_triangle(P, P1, P2, P3)
%IS_WITHIN_TRIANGLE whether point P lies within triangle

% P: point
% c: corners of triangle, each column containing x, y, z coords
% c updated to P1, P2, P3

% Barycentric

% x = pt(1,:);
% y = pt(2,:);
% 
% x1 = p(1,1);
% x2 = p(1,2);
% x3 = p(1,3);
% 
% y1 = p(2,1);
% y2 = p(2,2);
% y3 = p(2,3);
% 
% W1 = ((y2 - y3).*(x - x3) + (x3 - x2).*(y - y3))./((y2 - y3).*(x1 - x3) + (x3 - x2).*(y1 - y3));
% W2 = ((y3 - y1).*(x - x3) + (x1 - x3).*(y - y3))./((y2 - y3).*(x1 - x3) + (x3 - x2).*(y1 - y3));
% 
% W3 = 1 - W1 - W2;
% 
% W = [W1 W2 W3];
% 
% within = all(W > 0 & W < 1);

% Rewrite of above
% https://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates

v0 = P2 - P1; 
v1 = P3 - P1;
v2 = P - P1;
% Replacing dot(vi, vj, 2) for speed
d00 = sum(v0 .* v0, 2);
d01 = sum(v0 .* v1, 2);
d11 = sum(v1 .* v1, 2);
d20 = sum(v2 .* v0, 2);
d21 = sum(v2 .* v1, 2);

denom = d00 .* d11 - d01 .* d01;

v = (d11 .* d20 - d01 .* d21) ./ denom;
w = (d00 .* d21 - d01 .* d20) ./ denom;
u = 1 - v - w;

within = u >= 0 & u <= 1 & v >= 0 & v <= 1 & w >= 0 & w <= 1;

% W. Heidrich, Journal of Graphics, GPU, and Game Tools,Volume 10, Issue 3, 2005

% if nargin == 2
%     
%     P3 = P1(3,:);
%     P2 = P1(2,:);
%     P1 = P1(1,:);
% end
% 
% u = P2 - P1;
% v = P3 - P1;
% w = P - P1;
% 
% n = crossmat(u,v);
% cuw = crossmat(u,w);
% cwv = crossmat(w,v);
% 
% for i = size(n, 1):-1:1
%     
%     if any(n(i,:))
%         
%         gamma = (cuw(i,:) .* n(i,:)) / (n(i,:).^2);
%         beta = (cwv(i,:) .* n(i,:)) / (n(i,:).^2);
%         
%         alpha = 1 - gamma - beta;
%         
%         abg = [alpha beta gamma];
%         
%         within(i) = all(abg >= 0 & abg <= 1);
%     else
%         within(i) = false;
%     end
% end

% https://math.stackexchange.com/questions/4322/check-whether-a-point-is-within-a-3d-triangle

% A = crossmat(P2 - P1, P3 - P1);
% C = crossmat(P3 - P, P1 - P);
% D = crossmat(P1 - P, P2 - P);
% 
% k1 = C/A;
% k2 = D/A;
% 
% within = all([k1 k2] >= -0.1);

