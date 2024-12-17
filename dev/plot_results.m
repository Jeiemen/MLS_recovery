function plot_results(folder, results_name, xy, magnitude, format, geometry)
arguments
    folder
    results_name
    xy
    magnitude
    format = '-dpng'
    geometry = []
end

mkdir(folder);

if isempty(geometry)
    h = figure;
else
    plot_geometry(geometry);
    h = gcf();
end

hold on;
scatter(xy(:,1), xy(:,2), [], magnitude, 'filled'); 
colorbar;
hold off
axis equal
title(results_name)
set(h,'NumberTitle','off','Name',results_name,...
    'MenuBar','none','ToolBar','figure','DockControls','off');
hgsave(h,[folder '\' results_name])
print(h,[folder '\' results_name], format);
% close(h);