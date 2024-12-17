%% Find the closest point to a given point x over the boundary of the CYL
% WARNING: This is only implemented to search for lines and arc segments
function [boundary_xy,normal_angle,boundary_id] = find_boundary(xy, r, geo_data)
% Inputs:
% - xy: 2D test point coordinates
% - r: search radius
% - geo_data: geometry data file (see test_problem json example)

% Outputs:
% - boundary_xy: nearest boundary points to input coordinates (max. 1 per
% geometric entity [n x 2]
% - normal_angle: boundary outwards normal vector at boundary_xy location
% [n x 2]
% - boundary_id: geometric_entity ID corresponding to each entry in
% boundary_xy and normal_angle [n]

n_bound = numel(geo_data);
boundary_xy = zeros(n_bound, 2);
normal_angle = zeros(n_bound, 2);
boundary_id = 1:n_bound;
d = zeros(n_bound, 1);

for i_bound = 1:n_bound
    i_type = geo_data(i_bound).type;
    i_data = geo_data(i_bound).properties;
    switch i_type
        case 'line'
            [close_xy, dist, normal] = closestPointLineSegment(...
                i_data.p1', i_data.p2', xy);
        case 'arc'
            [close_xy, dist, normal] = closestPointArc(...
                i_data.p1', i_data.p2', i_data.radius, i_data.normal_dir, xy);
        otherwise
            error('Geometry type %s is not implemented', i_type)
    end
    boundary_xy(i_bound, :) = close_xy;
    normal_angle(i_bound, :) = normal;
    d(i_bound) = dist;
end

% Filter with patch support radius
index = find(d<=r);
boundary_xy = boundary_xy(index, :);
normal_angle = normal_angle(index, :);
boundary_id = boundary_id(index);
end

function [closestPoint, dist, normal] = closestPointLineSegment(p1, p2, x)
    v = p2 - p1;
    w = x - p1;
    t = (w * v') / (v * v');
    t = max(0, min(1, t));
    closestPoint = p1 + t * v;
    dist = sqrt(sum((x - closestPoint).^2));
    normal = [v(2), -v(1)];
    normal = normal / sqrt(sum(normal.^2));
end

function [closestPoint, distance, normal_p] = closestPointArc(p1, p2, r, normal_direction, x)
    % Arc line defined by ini-end points p1, p2, radius r and normal vector
    % direction ('in', 'out'). 
    % p1 and p2 are defined such that the arc is counter-clockwise
    midPoint = (p1 + p2) / 2;
    chord = p2 - p1;
    chordLength = sqrt(sum(chord.^2));
    centerOffset = sqrt(r^2 - (chordLength / 2)^2);
    normal = [-chord(2), chord(1)] / chordLength; % Perpendicular direction
    center = midPoint + centerOffset * normal;
    % Angles of p1, p2, and x relative to the circle center
    theta1 = atan2(p1(2) - center(2), p1(1) - center(1));
    theta2 = atan2(p2(2) - center(2), p2(1) - center(1));
    thetaX = atan2(x(2) - center(2), x(1) - center(1));
    % Ensure theta1 < theta2 for a counter-clockwise arc
    if theta2 < theta1
        theta2 = theta2 + 2 * pi;
    end
    if thetaX < theta1
        thetaX = thetaX + 2 * pi;
    end
    % Clamp thetaX to the range [theta1, theta2]
    clampedTheta = max(theta1, min(thetaX, theta2));
    normal_p = [cos(clampedTheta), sin(clampedTheta)];
    % Closest point on the arc
    closestPoint = center + r * normal_p;
    % Distance from x to the closest point
    distance = norm(x - closestPoint);
    % output normal vector
    if strcmp(normal_direction, 'in')
        normal_p = -normal_p;
    end
end