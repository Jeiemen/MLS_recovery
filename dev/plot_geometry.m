function plot_geometry(geometry)
figure;
hold on;
for i = 1:numel(geometry)
    i_type = geometry(i).type;
    i_data = geometry(i).properties;
    switch i_type
        case "line"
            plot([i_data.p1(1) i_data.p2(1)], [i_data.p1(2) i_data.p2(2)],'-k', 'LineWidth', 3)
        case "arc"
            plot_arc(i_data);
        otherwise
            error('Geometry type %s is not implemented', i_type)
    end
end
hold off;
end

function arcPoints = plot_arc(data, n_points)
    p1 = data.p1';
    p2 = data.p2';
    r = data.radius;
    if nargin == 1
        n_points = 50;
    end
    % Compute the center of the circle
    midPoint = (p1 + p2) / 2;
    chord = p2 - p1;
    chordLength = norm(chord);
    centerOffset = sqrt(r^2 - (chordLength / 2)^2);
    normal = [-chord(2), chord(1)] / chordLength; % Perpendicular direction
    center = midPoint + centerOffset * normal;
    % Angles of p1 and p2 relative to the circle center
    theta1 = atan2(p1(2) - center(2), p1(1) - center(1));
    theta2 = atan2(p2(2) - center(2), p2(1) - center(1));
    % Ensure counter-clockwise order
    if theta2 < theta1
        theta2 = theta2 + 2 * pi;
    end
    % Interpolate angles between theta1 and theta2
    angles = linspace(theta1, theta2, n_points);
    % Compute points on the arc
    arcPoints = center + r * [cos(angles(:)), sin(angles(:))];
    plot(arcPoints(:, 1), arcPoints(:, 2), '-k', 'LineWidth', 3);
end