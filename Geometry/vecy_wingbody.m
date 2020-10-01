function [wing, lbody, ubody] = vecy_wingbody(wing, body)

cont = true;
attempt = 1;
% Alterations by 10%
anchor = 0.1;
body_dim = [body.Length, body.Height];
prev_exit = [];

cut_chord = [1, 2; 3, 4];

while cont
    
    [up.joint, up.surf_id, up.line_id, up.exit] = line_surface(wing.upper, body);
    [lo.joint, lo.surf_id, lo.line_id, lo.exit] = line_surface(wing.lower, body);
    
    ucon = any(isnan(up.joint), 'all');
    lcon = any(isnan(lo.joint), 'all');
    
    exit = [upExit, loExit];
    
    if ucon || lcon
        %% Fixer stuff
        % 1/2: Increase/Decrease x-offset
        % 3/4: Increase/Decrease z-offset
        
        % Reduce chord if line_surface exited because of 1 & 2 or 3 & 4
        cut = all(ismember(cut_chord, [exit, prev_exit]), 2);
        con = [isempty(exit), any(cut), attempt > 5];
        
        if any(con)
            
            wing.chord(1) = wing.chord(1) * (1 - anchor);
        end
        
        if ~con(1)
            % Have to do it here, unique([]) = 0x1 not 0x0
            exit = unique(exit);
        end
        
        for i = exit
            
            half = ceil(i/2);
            
            if i == 5
                
                lbody = body;
                %                 geometry.plotter(lbody, wing)
                wing = [];
                ubody = [];
                return
                
            elseif i == 0
                
                continue
                
            elseif ~cut(half) || i == 1
                
                % Only alter offsets if one of [1, 2] or one of [3, 4] is
                % true. Always increase offset if i == 1
                if mod(i, 2)
                    
                    fix = anchor;
                else
                    fix = -anchor;
                end
                
                % Automatically selects offset and direction
                wing.offset(half) = wing.offset(half) + ...
                    fix * body_dim(half);
            end
        end
        
        if isequal(prev_exit, exit)
            
            attempt = attempt + 1;
        else
            attempt = 1;
        end
        
        prev_exit = exit;
        
        if wing.chord(1) < 0.1 * body.Length
            %% TODO: Make sure only coming in here in worst case scenarios
            %             geometry.plotter(body, wing)
            wing = [];
            lbody = [];
            ubody = [];
            return
        end
    else
        cont = false;
    end
end

%% Redefine wing

wing = copy(wing);

x = wing.x;
y = wing.y;
z = wing.z;

[~, dim] = size(x);
remove = false(1, dim);

for i = 1:dim/2
    
    replace = i <= upLine_id;
    
    if all(replace) && i > 1
        
        remove(i) = true;
        
    elseif any(replace)
        
        x(replace, i) = up.joint(replace, 1);
        y(replace, i) = up.joint(replace, 2);
        z(replace, i) = up.joint(replace, 3);
    else
        break
    end
end

% Same as above for lower surface, have to use separate increment counter
% since i starts at last point on wing wrap
j = 1;
for i = dim:-1:dim/2 + 1
    
    replace = j <= lo.line_id;
    
    if all(replace) && j > 1
        
        remove(i) = true;
        
    elseif any(replace)
        
        x(replace, i) = lo.joint(replace, 1);
        y(replace, i) = lo.joint(replace, 2);
        z(replace, i) = lo.joint(replace, 3);
    else
        break
    end
    
    j = j + 1;
end

x(:,remove) = [];
y(:,remove) = [];
z(:,remove) = [];

wing.x = x;
wing.y = y;
wing.z = z;

%% Redefine body

x = body.x(:,1);
y = body.y;
z = body.z;

%% UPPER AND LOWER
f_split_id = up.surf_id(1, 2);
b_split_id = up.surf_id(end, 2);

% Last panel upper joint interferes with (radially)
uradbound = max(up.surf_id(:,2));
% First panel lower joint interferes with
lradbound = min(lo.surf_id(:,2));

xu_int = up.joint(:,1);
xl_int = lo.joint(:,1);

% Including one extra point for interpolation below
yu_int = interp1(x(:,1), y(:, 1:uradbound+1), xu_int);
zu_int = interp1(x(:,1), z(:, 1:uradbound+1), xu_int);
yl_int = interp1(x(:,1), y(:, lradbound-1:end), xl_int);
zl_int = interp1(x(:,1), z(:, lradbound-1:end), xl_int);

yz_le = up.joint(1, [2 3]);
yzu_te = up.joint(end, [2 3]);
yzl_te = lo.joint(end, [2 3]);

ids = [0 1] + f_split_id;

%% Leading edge line
% Same for upper and lower
line = [yu_int(1, ids); zu_int(1, ids)]';
t = position_on_line(line(1,:), line(2,:), yz_le);

y_le = sum([1-t t] .* y(:, ids), 2);
z_le = sum([1-t t] .* z(:, ids), 2);

%% Trailing edge line
% Can be different

ids = [0 1] + b_split_id;

line = [yu_int(end, ids); zu_int(end, ids)]';
t = position_on_line(line(1,:), line(2,:), up.joint(end, [2 3]));

yu_te = sum([1-t t] .* y(:, ids), 2);
zu_te = sum([1-t t] .* z(:, ids), 2);

if isequal(yzu_te, yzl_te)
    
    yl_te = yu_te;
    zl_te = zu_te;
else
    line = [yl_int(end, [1 2]); zl_int(end, [1 2])]';
    t = positionOnLine(line(1,:), line(2,:), lo.joint(end, [2 3]));
    
    yl_te = sum([1-t t] .* y(:, ids), 2);
    zl_te = sum([1-t t] .* z(:, ids), 2);
end

xuedges = [min(xu_int), max(xu_int)];
xledges = [min(xl_int), max(xl_int)];

le = [x y_le z_le];
u_te = [x yu_te zu_te];
l_te = [x yl_te zl_te];

up.joint = [le(le(:,1) < xuedges(1), :);
    up.joint;
    u_te(u_te(:,1) > xuedges(2), :)];

lo.joint = [le(le(:,1) < xledges(1), :);
    lo.joint;
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
rad_joint = abs(atan2d(up.joint(:,2), up.joint(:,3)));
[~, udim] = size(rad);

for i = 1:udim
    
    replace = rad(:,i) > rad_joint;
    yu(replace, i) = up.joint(replace, 2);
    zu(replace, i) = up.joint(replace, 3);
end

rad = abs(atan2d(yl, zl));
rad_joint = abs(atan2d(lo.joint(:,2), lo.joint(:,3)));
[~, ldim] = size(rad);

for i = 1:ldim
    
    replace = rad(:,i) < rad_joint;
    yl(replace, i) = lo.joint(replace, 2);
    zl(replace, i) = lo.joint(replace, 3);
end

% Replacing extra point with joint
yu(:, end) = up.joint(:,2);
zu(:, end) = up.joint(:,3);
yl(:, 1) = lo.joint(:,2);
zl(:, 1) = lo.joint(:,3);

ubody = copy(body);
lbody = copy(body);

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

function [joint, surf_id, seg_id, exit] = line_surface(lines, surface)
% Some code optimisation specific to expected input. Assumes each line
% proceeds the last in terms of x-coordinate

[nLines, nPoints, ~] = size(lines);

joint = nan(nLines, 3);
surf_id = zeros(nLines, 2);
seg_id = zeros(nLines, 1);

data = surface.quad_data;
[s1, s2, ~] = size(data.centre);
tri_id = surface.tri_id();

x = surface.x;
y = surface.y;
z = surface.z;

xTri = x(tri_id);
yTri = y(tri_id);
zTri = z(tri_id);

centre = data.centre;
mag = (centre(:,:,2).^2 + centre(:,:,3).^2).^0.5;

% r defined at mean z for each x ie. origin, y=0
r = abs(mag - mean(centre(:,:,3), 2));

centre = reshape(centre, [], 3, 1);
norm = reshape(data.norm, [], 3, 1);

body_rad = mean(r, 1);

[~, y_order] = sort(body_rad, 'descend');

lead = lines(1,:,:) - lines(1,1,:);
lead_rad = (lead(:,:,2).^2 + lead(:,:,3).^2).^0.5;
c = find(any(body_rad' >= lead_rad), 1, 'last');

% Order to go through line segments
nSegs = nPoints - 1;
seg_order = get_order(nSegs, min(c, nSegs), 2);

%% Initial check to see if either first or second point lie in 2D polygon
% If lines lies within body at midline then they have to pass through the
% surface (assuming infinite length)
% Any lines are completely ahead or behind surface
% Any lines that are completely above or below surface

for i = 2:-1:1
    
    [~, xExit(:,i), yExit(:,i)] = is_within_polygon(x(:,1), z(:,[1 end]), lines(:,i,1), lines(:,i,3));
end

xOut = diff(xExit, [], 2) == 0;
yOut = diff(xExit, [], 2) == 0 & ~xOut;

flag = [xExit(:,1) yExit(:,1)];
exit = unique([flag(xOut, 1); flag(yOut, 2)])';
if any(exit), return; end

%% Failures at front/end can still slip through above exit conditions, so
% best to check directly that they work before trying all lines

%% Main loop
for i = [nLines, 1:nLines-1]
    
    if i > 1 && i < nLines
        % Reorder to have closest panels first
        % y_order = getorder(s2, k(i-1));
        seg_order = get_order(nSegs, c(i-1), 2);
    end
    
    [joint(i,:), j(i,:), k(i,:), c(i,:)] = main(lines(i,:,:));
    
    exit = exit_code(joint(i,:), j(i), k(i));
    if any(exit), return; end
end

surf_id = [j k];
seg_id = c;

    function [joint, j, k, c, flag] = main(line)
        
        joint = nan(1, 3);
        flag = 0;
        j = 0;
        k = 0;
        c = seg_order(1);
        seg_count = 1;
        rows = 1:s1;
        
        while true
            
            p1 = squeeze(line(:, c, :));
            p2 = squeeze(line(:, c+1, :));
            u = p1 - p2;

            % y-dimension loop
            for k = y_order
                
                % Radial dimension loop
                idq = rows + (k-1)*s1;
                
                M = centre(idq,:);
                n = norm(idq,:);
                
                % Get intersection point
                inter = new_line_plane_intersection(u', p1', n, M);
                
                for ii = 0:1
                    
                    P1 = [xTri(idq*2+ii-1, 1) yTri(idq*2+ii-1, 1) zTri(idq*2+ii-1, 1)];
                    P2 = [xTri(idq*2+ii-1, 2) yTri(idq*2+ii-1, 2) zTri(idq*2+ii-1, 2)];
                    P3 = [xTri(idq*2+ii-1, 3) yTri(idq*2+ii-1, 3) zTri(idq*2+ii-1, 3)];
                    
                    inTri = is_within_triangle(inter, P1, P2, P3);
                    
                    if any(inTri)
                        
                        inInter = inter(inTri, :);
                        j = find(inTri);
                        
                        for s = 1:size(inInter, 1)
                            % Does point lie within line segment bounds
                            [onLine, ~] = on_line_segment(p1, p2, inInter(s, :)');
                            
                            if onLine
                                %% TODO: May need to loop all the way
                                % instead of return, save intersection
                                % and continue, compare based on polar
                                % r
                                joint = inInter(s, :);
                                j = j(s);
                                return
                            end
                        end
                    end
                end
            end
            
            % Go to next segment on line
            seg_count = seg_count+ 1;
            
            if seg_count > nSegs
                
                flag = true;
                j = nan;
                k = nan;
                c = nan;
                return
            else
                c = seg_order(seg_count);
            end
        end
    end

    function exit = exit_code(joint, j, k)
        
        exit = 0;
        
        if any(isnan(joint))
            
            if any(flag(i,:))
                
                exit = flag(i,:);
            else
                exit = 5;
            end
        else
            con = [j == [1 s1], k == [1 s2]];
            
            if any(con)
                
                arr = 1:numel(con);
                exit = arr(con);
            end
        end
    end
end

function order = get_order(dim, begin, method)
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
    n = min(numel(one), numel(two));
    order = [reshape([one(1:n); two(1:n)], 1, []), one(n+1:end), ...
        two(n+1:end)];
end
end

function t = position_on_line(p1, p2, p)

short_line = [p1; p];
long_line = [p1; p2];

% Magnitude (line length)
d = sum(diff(long_line).^2).^0.5;
dt = sum(diff(short_line).^2).^0.5;

t = dt/d;
end