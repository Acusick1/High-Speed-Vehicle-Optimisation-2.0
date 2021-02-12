function [wing, lbody, ubody] = vec_wingbody(wing, body)

cont        = true;
attempt     = 1;
anchor      = 0.1; % Alterations by 10%

if sum(wing.span) > 1
    % Renormalise so total wingspan < below dimensionaliser
    wing.span = wing.span/sum(wing.span);
end

wing.chord  = wing.chord  .* body.length;
wing.span   = wing.span   .* 0.5 * (body.length + body.width/2);
wing.offset = wing.offset .* [body.length, body.height];

while cont
    
    wing = wing.dogenerate();
    
    [ujoint, usurf_id, uline_id, uexit] = linesurface(wing.upper, body);
    [ljoint, lsurf_id, lline_id, lexit] = linesurface(wing.lower, body);
    
    exit = [uexit, lexit];
    
    if any(exit)
        
        % Do not use elseifs, alterations should be stand-alone
        % Thus one error cannot result in two alterations also
        if any(exit == 1)
            
            wing.offset(1) = wing.offset(1) * (1 + anchor);
        end
        
        zOut = [any(exit == 3), any(exit == 4)];
        
        if sum(zOut) == 1
            
            % Will tend towards midline
            wing.offset(2) = wing.offset(2) * (1 - anchor);
        end
        
        if any(exit == 2) || all(zOut) || attempt > 10
            
            wing.chord(1) = wing.chord(1) * (1 - anchor);
        end
        
        attempt = attempt + 1;
        
        if wing.chord(1) < 0.1 * body.length
            % geometry.plotter(body, wing)
            wing = [];
            % Only needed if handle
            % lbody = copy(body);
            lbody = body;
            ubody = [];
            return
        end
    else
        cont = false;
    end
end

%% Redefine wing
% Only needed if handle
% wing = copy(wing);

x = wing.x;
y = wing.y;
z = wing.z;

[~, dim] = size(x);
remove = false(1, dim);

for i = 1:dim/2
    
    replace = i <= uline_id;
        
    if any(replace)
        
        x(replace, i) = ujoint(replace, 1);
        y(replace, i) = ujoint(replace, 2);
        z(replace, i) = ujoint(replace, 3);
        
        if i > 1 && all(replace)
        
            remove(i) = true;
        end
    else
        break
    end
end

% Same as above for lower surface, have to use separate increment counter
% since i starts at last point on wing wrap
j = 1;
for i = dim:-1:dim/2 + 1
    
    replace = j <= lline_id;
        
    if any(replace)
        
        x(replace, i) = ljoint(replace, 1);
        y(replace, i) = ljoint(replace, 2);
        z(replace, i) = ljoint(replace, 3);
        
        if j > 1 && all(replace)
        
            remove(i) = true;
        end
    else
        break
    end
    
    j = j + 1;
end

% Commented out to maintain consistency ie ensuring number of columns same
% for upper and lower surfaces
% Removing duplicate point sets
% x(:,remove) = [];
% y(:,remove) = [];
% z(:,remove) = [];

wing.x = x;
wing.y = y;
wing.z = z;

%% Redefine body

x = body.x(:,1);
y = body.y;
z = body.z;

%% UPPER AND LOWER
f_split_id = usurf_id(1, 2);
b_split_id = usurf_id(end, 2);

% Last panel upper joint interferes with (radially)
uradbound = min(max(usurf_id(:,2)), size(usurf_id, 1) - 1);
% First panel lower joint interferes with
lradbound = max(min(lsurf_id(:,2), 2));

xu_int = ujoint(:,1);
xl_int = ljoint(:,1);

% Including one extra point for interpolation below
yu_int = interp1(x(:,1), y(:, 1:uradbound+1), xu_int);
zu_int = interp1(x(:,1), z(:, 1:uradbound+1), xu_int);
yl_int = interp1(x(:,1), y(:, lradbound-1:end), xl_int);
zl_int = interp1(x(:,1), z(:, lradbound-1:end), xl_int);

yz_le = ujoint(1, [2 3]);
yzu_te = ujoint(end, [2 3]);
yzl_te = ljoint(end, [2 3]);

ids = [0 1] + f_split_id;

%% Leading edge line
% Same for upper and lower
line = [yu_int(1, ids); zu_int(1, ids)]';
t = positionOnLine(line(1,:), line(2,:), yz_le);

y_le = sum([1-t t] .* y(:, ids), 2);
z_le = sum([1-t t] .* z(:, ids), 2);

%% Trailing edge line
% Can be different

ids = [0 1] + b_split_id;

line = [yu_int(end, ids); zu_int(end, ids)]';
t = positionOnLine(line(1,:), line(2,:), ujoint(end, [2 3]));

yu_te = sum([1-t t] .* y(:, ids), 2);
zu_te = sum([1-t t] .* z(:, ids), 2);

if isequal(yzu_te, yzl_te)
    
    yl_te = yu_te;
    zl_te = zu_te;
else
    line = [yl_int(end, [1 2]); zl_int(end, [1 2])]';
    t = positionOnLine(line(1,:), line(2,:), ljoint(end, [2 3]));
    
    yl_te = sum([1-t t] .* y(:, ids), 2);
    zl_te = sum([1-t t] .* z(:, ids), 2);
end

xuedges = [min(xu_int), max(xu_int)];
xledges = [min(xl_int), max(xl_int)];

le = [x y_le z_le];
u_te = [x yu_te zu_te];
l_te = [x yl_te zl_te];

ujoint = [le(le(:,1) < xuedges(1), :);
    ujoint;
    u_te(u_te(:,1) > xuedges(2), :)];

ljoint = [le(le(:,1) < xledges(1), :);
    ljoint;
    l_te(l_te(:,1) > xledges(2), :)];

con = x < xuedges(1) | x > xuedges(2);
xu = [x(con, :); xu_int];
yu = [y(con, 1:uradbound+1); yu_int];
zu = [z(con, 1:uradbound+1); zu_int];

con = x < xledges(1) | x > xledges(2);
xl = [x(con, :); xl_int];
yl = [y(con, lradbound-1:end); yl_int];
zl = [z(con, lradbound-1:end); zl_int];

[xu, idu] = sort(xu);
[xl, idl] = sort(xl);

yu = yu(idu,:);
zu = zu(idu,:);
yl = yl(idl,:);
zl = zl(idl,:);

rad = abs(atan2d(yu, zu));
rad_joint = abs(atan2d(ujoint(:,2), ujoint(:,3)));
[~, udim] = size(rad);

for i = 1:udim
    
    replace = rad(:,i) > rad_joint;
    yu(replace, i) = ujoint(replace, 2);
    zu(replace, i) = ujoint(replace, 3);
end

rad = abs(atan2d(yl, zl));
rad_joint = abs(atan2d(ljoint(:,2), ljoint(:,3)));
[~, ldim] = size(rad);

for i = 1:ldim
    
    replace = rad(:,i) < rad_joint;
    yl(replace, i) = ljoint(replace, 2);
    zl(replace, i) = ljoint(replace, 3);
end

% Replacing extra point with joint
yu(:, end) = ujoint(:,2);
zu(:, end) = ujoint(:,3);
yl(:, 1) = ljoint(:,2);
zl(:, 1) = ljoint(:,3);
%% TODO: Only needed if handle
% ubody = copy(body);
% lbody = copy(body);
ubody = body;
lbody = body;

ubody.x = repmat(xu, 1, udim);
ubody.y = yu;
ubody.z = zu;

lbody.x = repmat(xl, 1, ldim);
lbody.y = yl;
lbody.z = zl;

% geometry.plotter(lbody, ubody, wing)

% figure(gcf)
% hold on
% plot3(ujoint(:,1), ujoint(:,2), ujoint(:,3), 'b*')
% plot3(ljoint(:,1), ljoint(:,2), ljoint(:,3), 'b*')
% hold off
end

function [joint, surf_id, seg_id, exit] = linesurface(lines, surface)
% Some code optimisation specific to expected input. Assumes each line
% proceeds the last in terms of x-coordinate

[nLines, nPoints, ~] = size(lines);

joint = nan(nLines, 3);
surf_id = zeros(nLines, 2);
seg_id = zeros(nLines, 1);

data = surface.quad_data;
[s1, s2, ~] = size(data.centre);
tri_id = surface.tri_data.id;

%% GPU TEST

% x = gpuArray(surface.x);
% y = gpuArray(surface.y);
% z = gpuArray(surface.z);

x = surface.x;
y = surface.y;
z = surface.z;

%% Initial check to see if either first or second point lie in 2D polygon
% If lines lies within body at midline then they have to pass through the
% surface (assuming infinite length)
% Any lines are completely ahead or behind surface
% Any lines that are completely above or below surface

for i = 2:-1:1
    
    [~, xexit(:,i), yexit(:,i)] = is_within_polygon(x(:,1), z(:,[1 end]), lines(:,i,1), lines(:,i,3));
end

xout = diff(xexit, [], 2) == 0;
yout = diff(yexit, [], 2) == 0 & ~xout;

flag = [xexit(:,1) yexit(:,1)];
exit = unique([flag(xout, 1); flag(yout, 2)])';
if any(exit), return; end

%% Main init

xtri = x(tri_id);
ytri = y(tri_id);
ztri = z(tri_id);

centre = data.centre;
mag = (centre(:,:,2).^2 + centre(:,:,3).^2).^0.5;

% r defined at mean z for each x ie. origin, y=0
r = abs(mag - mean(centre(:,:,3), 2));
r = r(:);

[~, idq] = sort(r, 'descend');
idt = reshape([idq*2-1 idq*2]', [], 1);

P1 = [xtri(idt, 1) ytri(idt, 1) ztri(idt, 1)];
P2 = [xtri(idt, 2) ytri(idt, 2) ztri(idt, 2)];
P3 = [xtri(idt, 3) ytri(idt, 3) ztri(idt, 3)];

centre = reshape(centre, [], 3, 1);
norm = reshape(data.norm, [], 3, 1);

M = centre(idq,:);
n = norm(idq,:);

lead = lines(1,:,:) - lines(1,1,:);
lead_rad = (lead(:,:,2).^2 + lead(:,:,3).^2).^0.5;
c = find(any(r >= lead_rad), 1, 'last');

% Order to go through line segments
nSegs = nPoints - 1;
seg_order = getorder(nSegs, min(c, nSegs), 2);

nq = 1:length(idq);
hack = reshape([nq; nq], [], 1);

%% Within triangle init

v0 = P2 - P1; 
v1 = P3 - P1;
d00 = sum(v0 .* v0, 2);
d01 = sum(v0 .* v1, 2);
d11 = sum(v1 .* v1, 2);
denom = d00 .* d11 - d01 .* d01;

%% Main loop
for i = [nLines, 1:nLines-1]
    
    if i > 1 && i < nLines
        % Reorder to have closest panels first
        seg_order = getorder(nSegs, c(i-1), 2);
    end
    
    [joint(i,:), j, k, c(i,:)] = main(lines(i,:,:));
    surf_id(i,:) = [j k];
    
    exit = exitcode(joint(i,:), j, k);
    if any(exit), return; end
end

seg_id = c;

    function [joint, j, k, c, flag] = main(line)
        
        joint = nan(1, 3);
        flag = 0;
        j = 0;
        k = 0;
        c = seg_order(1);
        seg_count = 1;
        
        while true
            
            p1 = squeeze(line(:, c, :));
            p2 = squeeze(line(:, c+1, :));
            u = p1 - p2;
            
            % Get intersection point
            I = new_line_plane_intersection(u', p1', n, M);
            I = I(hack, :);
            
            inTri = iswithintriangle(I);
            
            if any(inTri)
                
                I_in = I(inTri, :);
                ids = idt(inTri, :);
                
                % Does point lie within line segment bounds
                [onLine, ~] = on_line_segment(p1', p2', I_in);
                
                if any(onLine)
                    
                    onLine = find(onLine, 1);
                    
                    joint = I_in(onLine, :);
                    id = ids(onLine);
                    
                    if mod(id, 2)
                        
                        id = id + 1;
                    end
                    
                    [j, k] = ind2sub([s1, s2], id/2);
                    return
                end
            end
            
            % Go to next segment on line
            seg_count = seg_count+ 1;
            
            if seg_count > nSegs
                
                flag = true;
                return
            else
                c = seg_order(seg_count);
            end
        end
    end

    function within = iswithintriangle(P)
        
        v2 = P - P1;
        % Replacing dot(vi, vj, 2) for speed
        d20 = sum(v2 .* v0, 2);
        d21 = sum(v2 .* v1, 2);

        v = (d11 .* d20 - d01 .* d21) ./ denom;
        w = (d00 .* d21 - d01 .* d20) ./ denom;
        u = 1 - v - w;

        within = u >= 0 & u <= 1 & v >= 0 & v <= 1 & w >= 0 & w <= 1;
    end

    function exit = exitcode(joint, j, k)
        
        exit = 0;
        
        if any(isnan(joint))
            
            if any(flag(i,:))
                
                exit = flag(i,:);
            else
                exit = 5;
            end
        else
            con = [j == [1 s1] k == [1 s2]];
            
            if any(con)
                
                arr = 1:numel(con);
                exit = arr(con);
            end
        end
    end
end

function order = getorder(dim, begin, method)
%% Order indices from defined start point
% Function determining which order to search indices based on last
% successful index

if nargin < 3 || isempty(method)
    
    method = 1;
end

one = begin:dim;
two = begin-1:-1:1;

if method == 1
    % Increasing decreasing sequence [i, i+1, ... dim, i-1, i-2, ... 1]
    order = [one, two];
    
elseif method == 2
    % Interleave to get sequence [i, i-1, i+1, i-2...]
    N = min(numel(one),numel(two));
    order = [reshape([one(1:N); two(1:N)], 1, []), one(N+1:end), ...
        two(N+1:end)];
end
end

function t = positionOnLine(p1, p2, p)

shortline = [p1; p];
longline = [p1; p2];

% Magnitude (line length)
d = sum(diff(longline).^2).^0.5;
dt = sum(diff(shortline).^2).^0.5;

t = dt/d;
end