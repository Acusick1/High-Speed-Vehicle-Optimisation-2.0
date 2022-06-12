function [wing, lbody, ubody, extWing] = vec_wingbody(wing, body)

body = body.generate();

% Removing body boundary from interpolant so that wing cannot alter it
bx = body.x(2:end-1, 2:end-1);
by = body.y(2:end-1, 2:end-1);
bz = body.z(2:end-1, 2:end-1);

Fy = scatteredInterpolant(bx(:), bz(:), by(:));
Fy.ExtrapolationMethod = 'none';

cont        = true;
attempt     = 1;
anchor      = 0.1; % Alterations by 10%

withinLength = wing.chord(1) + wing.offset(1);

if withinLength > 1

    wing.chord(1) = wing.chord(1)/withinLength;
    wing.offset(1) = wing.offset(1)/withinLength;
end

if sum(wing.span) > 1
    % Renormalise so total wingspan < below dimensionaliser
    wing.span = wing.span/sum(wing.span);
end

maxSpan = 0.5 * body.length;

wing.chord  = wing.chord  .* body.length;
% Body width halved again to give symmetry axis width
wing.span   = wing.span .* maxSpan;
wing.offset(1) = wing.offset(1) .* body.length;

% Ratio defining minimum local width wing can be placed at
yMinOffsetRatio = 0.3;

zOffset = get_zOffsets(body, wing.offset(1));
yOffset = Fy(repmat(wing.offset(1), size(zOffset)), zOffset);
offCount = find(yOffset > yMinOffsetRatio * max(yOffset), 1);

wing.offset(2) = zOffset(offCount);
quit = false;

while cont
    
    midDiff = body.length/2 - (wing.offset(1) + wing.chord(1)/2);
    if midDiff > 0 && attempt < 20
        
    	wing.offset(1) = wing.offset(1) + midDiff;
        zOffset = get_zOffsets(body, wing.offset(1));
        yOffset = Fy(repmat(wing.offset(1), size(zOffset)), zOffset);
        
        offCount = find(yOffset > yMinOffsetRatio * max(yOffset), 1);
    end
    
    wing = wing.generate();
        
    extWing = wing;
    extWing.y = extWing.y + yOffset(offCount);
    % Geometry.plotter(body, temp)
    extWing = extWing.extendWithin();
    
    [ujoint, ~, uline_id, uexit] = linesurface(extWing.upper, body);
    [ljoint, ~, lline_id, lexit] = linesurface(extWing.lower, body);
    
    exit = unique([uexit, lexit]);
    % Prioritise zOffset first, then deal with x
    prioritise = any(exit == [3, 4]', 1);
    if any(prioritise)
        
        exit(~prioritise) = [];
    end
    
    if any(exit)
        
        % Do not use elseifs, alterations should be stand-alone
        % Thus one error cannot result in two alterations also
        if any(exit == 1)
            
            wing.offset(1) = wing.offset(1) * (1 + anchor);
            zOffset = get_zOffsets(body, wing.offset(1));
            yOffset = Fy(repmat(wing.offset(1), size(zOffset)), zOffset);
            
            offCount = find(yOffset > yMinOffsetRatio * max(yOffset), 1);
        end
        
        if any(exit == 3)
            
            if offCount < numel(zOffset)
                
                offCount = offCount + 1;
                wing.offset(2) = zOffset(offCount);
            
%             elseif wing.dihedral(1) > deg2rad(20)
%                 
%                 wing.dihedral(1) = wing.dihedral(1) - deg2rad(5);
            else
                wing.chord(1) = wing.chord(1) * (1 - anchor);
            end
            % Will tend towards midline
            % wing.offset(2) = wing.offset(2) * (1 - anchor);
        end
        
        % Wing should never be above body since its placed below mid-height
        % So if exit == 4, treat same as chord beyond body
        if any(exit == 2) || any(exit == 4) %||attempt > 10
            
            lead_sweep = wing.get_lead_sweep();
            
            if lead_sweep(1) < 0
                
                wing.trail_sweep(1) = wing.trail_sweep(1) * (1 - anchor);
            
            elseif wing.chord(1) > 0.5 * body.length || midDiff > 0 || all(any(exit' == [3, 4]))
                
                wing.chord(1) = wing.chord(1) * (1 - anchor);
            else
                wing.offset(1) = wing.offset(1) * (1 - anchor);
                zOffset = get_zOffsets(body, wing.offset(1));
                yOffset = Fy(repmat(wing.offset(1), size(zOffset)), zOffset);
                
                offCount = find(yOffset > yMinOffsetRatio * max(yOffset), 1);
            end
        end
        
        if any(exit == 5)
            
            wing.span(1) = wing.span(1) * (1 + anchor);
            
            if sum(wing.span) > maxSpan
                
                quit = true;
            end
        end
        
        attempt = attempt + 1;
        
        if mod(attempt, 500) == 0
            
            quit = true;
            % error('%s failed, cannot merge components', mfilename) 
        end
        
        if wing.chord(1) < 0.1 * body.length || quit
            % Geometry.plotter(body, wing)
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
wing = extWing;

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

% Only keep body points outwith wing chord
% Merge body, wing points and interpolate to find new discretisation
% Cannot use sort here, sometimes joint will not be monotonically
% increasing in x
% xunew = [x(x < min(wing.x(:,1))); wing.x(:,1); x(x > max(wing.x(:,1)))];
% xlnew = [x(x < min(wing.x(:,end))); wing.x(:,end); x(x > max(wing.x(:,end)))];
xunew = [x(x < wing.x(1,1)); wing.x(:,1); x(x > wing.x(end,1))];
xlnew = [x(x < wing.x(1,end)); wing.x(:,end); x(x > wing.x(end,end))];
yunew = interp1(x, y, xunew);
zunew = interp1(x, z, xunew);
ylnew = interp1(x, y, xlnew);
zlnew = interp1(x, z, xlnew);

%% UPPER AND LOWER
%% Leading edge (upper and lower)
le_id = find(wing.x(1,1) == xunew, 1);

lastAbove = find(zunew(le_id,:) - wing.z(1,1) < 0, 1) - 1;
lastAbove = max(min(lastAbove, size(zunew, 2) - 1), 1);
ids = [0 1] + lastAbove;

yz_le = [wing.y(1), wing.z(1)];

line = [yunew(le_id, ids); zunew(le_id, ids)]';
t = positionOnLine(line(1,:), line(2,:), yz_le);

y_le = sum([1-t t] .* yunew(:, ids), 2);
z_le = sum([1-t t] .* zunew(:, ids), 2);

yu_le = y_le(1:le_id-1);
zu_le = z_le(1:le_id-1);

%% Leading edge (lower)
le_id = find(wing.x(1,end) == xlnew, 1);

lastAbove = find(zlnew(le_id,:) - wing.z(1,1) < 0, 1) - 1;
lastAbove = max(min(lastAbove, size(zlnew, 2) - 1), 1);
ids = [0 1] + lastAbove;

yz_le = [wing.y(1), wing.z(1)];

line = [ylnew(le_id, ids); zlnew(le_id, ids)]';
t = positionOnLine(line(1,:), line(2,:), yz_le);

y_le = sum([1-t t] .* ylnew(:, ids), 2);
z_le = sum([1-t t] .* zlnew(:, ids), 2);

yl_le = y_le(1:le_id-1);
zl_le = z_le(1:le_id-1);

%% Trailing edge (upper)
te_id = find(wing.x(end, 1) == xunew, 1, 'last');

lastAbove = find(zunew(te_id,:) - wing.z(end, 1) < 0, 1) - 1;
lastAbove = max(min(lastAbove, size(zunew, 2) - 1), 1);
ids = [0 1] + lastAbove;

yz_te = [wing.y(end, 1), wing.z(end, 1)];

line = [yunew(te_id, ids); zunew(te_id, ids)]';
t = positionOnLine(line(1,:), line(2,:), yz_te);

yu_te = sum([1-t t] .* yunew(:, ids), 2);
zu_te = sum([1-t t] .* zunew(:, ids), 2);

yu_te = yu_te(te_id+1:end);
zu_te = zu_te(te_id+1:end);

%% Trailing edge (lower)
te_id = find(wing.x(end, end) == xlnew, 1, 'last');

lastAbove = find(zlnew(te_id,:) - wing.z(end, end) < 0, 1) - 1;
lastAbove = max(min(lastAbove, size(zlnew, 2) - 1), 1);
ids = [0 1] + lastAbove;

yz_te = [wing.y(end, end), wing.z(end, end)];

line = [ylnew(te_id, ids); zlnew(te_id, ids)]';
t = positionOnLine(line(1,:), line(2,:), yz_te);

yl_te = sum([1-t t] .* ylnew(:, ids), 2);
zl_te = sum([1-t t] .* zlnew(:, ids), 2);

yl_te = yl_te(te_id+1:end);
zl_te = zl_te(te_id+1:end);

yuJoint = [yu_le; wing.y(:,1); yu_te];
zuJoint = [zu_le; wing.z(:,1); zu_te];

ylJoint = [yl_le; wing.y(:,end); yl_te];
zlJoint = [zl_le; wing.z(:,end); zl_te];

rad = abs(atan2d(yunew, zunew));
rad_joint = abs(atan2d(yuJoint, zuJoint));

%% Upper and lower seperately defined by orginal body points and joints
% 2:end-1 to exclude collapsing edges (0 degrees)
udim = find(min(rad(2:end-1,:)) - max(rad_joint(2:end-1)) > 0, 1) - 1;

yub = yunew(:,1:udim);
zub = zunew(:,1:udim);
rad = rad(:,1:udim);

for i = 1:udim
    
    replace = rad(:,i) > rad_joint;
    yub(replace, i) = yuJoint(replace);
    zub(replace, i) = zuJoint(replace);
end

yub = [yub yuJoint];
zub = [zub zuJoint];

rad = abs(atan2d(ylnew, zlnew));
rad_joint = abs(atan2d(ylJoint, zlJoint));

% Lower
% 2:end-1 to exclude collapsing edges (0 degrees)
ldim = find(max(rad(2:end-1,:)) - min(rad_joint(2:end-1)) > 0, 1);

ylb = ylnew(:,ldim:end);
zlb = zlnew(:,ldim:end);
rad = rad(:,ldim:end);

for i = 1:size(ylb, 2)
    
    replace = rad(:,i) < rad_joint;
    ylb(replace, i) = ylJoint(replace);
    zlb(replace, i) = zlJoint(replace);
end

ylb = [ylJoint ylb];
zlb = [zlJoint zlb];

%% TODO: Only needed if handle
% ubody = copy(body);
% lbody = copy(body);
ubody = body;
lbody = body;

ubody.x = repmat(xunew, 1, size(yub, 2));
ubody.y = yub;
ubody.z = zub;

lbody.x = repmat(xlnew, 1, size(ylb, 2));
lbody.y = ylb;
lbody.z = zlb;

% Geometry.plotter(lbody, ubody, wing)

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
    
    [~, xexit(:,i), yexit(:,i)] = is_within_parallelogram(x(:,1), z(:,[1 end]), lines(:,i,1), lines(:,i,3));
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
    if any(exit)
        
        return
    end
end

% Last resort to say increase span, would be better to have condition
% triggering this within exitcode
if any(isnan(joint(:,1)))
    
    exit = 5;
    return
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
            I = line_plane_intersection_simple(u', p1', n, M);
            I = I(hack, :);
            
            inTri = iswithintriangle(I);
            
            if any(inTri)
                
                I_in = I(inTri, :);
                ids = idt(inTri, :);
                
                % Does point lie within line segment bounds
                [onLine, ~] = on_line_segment(p1', p2', I_in);
                
                if any(onLine)
                    
                    % If more than one, pick one with largest y
                    if sum(onLine) > 1
                        
                        onLine = I_in(:,2) == max(I_in(:,2));
                        
                        % Still get more than one point if intersection
                        % lies on panel line
                        if sum(onLine) > 1
                            
                            del = find(onLine);
                            onLine(del(2:end):end) = 0;
                        end
                    end
                    
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
%             else
%                 exit = 5;
            end
        else
%             con = [j == [1 s1] k == [1 s2]];
%             
%             if any(con)
%                 
%                 arr = 1:numel(con);
%                 exit = arr(con);
%             end
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

function offsets = get_zOffsets(body, xOffset)

offsetsNd = 0.1:0.05:0.9;

id = find(body.x(:,1) - xOffset > 0, 1) - 1;

x0 = body.x(id,:)';
x1 = body.x(id + 1,:)';
z0 = body.z(id,:)';
z1 = body.z(id + 1,:)';
zint = vec_interp(x0, x1, z0, z1, xOffset);

offsets = vec_interp(0, 1, min(zint), max(zint), offsetsNd);
end