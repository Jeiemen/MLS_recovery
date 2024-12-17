function plot_patch_support(xy_patch, xy_out, r_patch, geometry, bound_id)
figure;
hold on;
plot_geometry(geometry);
hold on;
rectangle('Position', [xy_out(1) - r_patch, xy_out(2) - r_patch,...
    2*r_patch, 2*r_patch],'Curvature',[1,1]);
plot(xy_out(1), xy_out(2), 'bo');
plot(xy_patch(:, 1), xy_patch(:, 2), 'r*');
if ~isempty(bound_id)
    plot(bound_id(:, 1), bound_id(:, 2), 'kx');
end
hold off
axis equal