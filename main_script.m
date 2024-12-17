restoredefaultpath
clearvars

%% Problem Data
problem_folder = './test_problem_cylinder';
% problem folder must contain the following functions (see test_problem_cylinder for examples):
% - find_boundary : find the closest point from boundary entities to input coordinates
% - neumann_conditions : evaluate neumann boundary conditions for a given
% coordinate and boundary segment
% - volume_load : body forces of the problem for a given coordinate

% Data files
geo_file = [problem_folder '/cylinder_geo.json'];
mesh_sample = 2;
data_file = [problem_folder '/data/FE_stress_mesh_' num2str(mesh_sample) '.mat'];
% mesh_sample = mesh iteration during a uniform mesh refinement process
% data_file contains, for each column, the following data referred to FE 
% quadrature points (9 per element): [coordinates_xy; integration_weight; stress_fe (x, y, xy)];

% Load data and geometry
path(problem_folder, path);
path('./dev', path);
geometry = read_geometry(geo_file);
data = load(data_file);
xy_data = data.data(1:2, :)';
w_data = data.data(3, :)';
stress_data = data.data(4:6, :)';

%% Query points
% output coordinates can be different from input data

% % suitable for cylinder problem
% n_rad = 30; 
% n_angle = 30;
% rad_coord = linspace(5, 20, n_rad);
% angle_coord = linspace(0, pi/2, n_angle);
% [rad_out, angle_out] = meshgrid(rad_coord, angle_coord);
% cil_output = [rad_out(:), angle_out(:)];
% xy_output = cil_output(:,1) .* [cos(cil_output(:,2)), sin(cil_output(:,2))];

% generic output can be the same input points
xy_output = xy_data; 

%% Run MLS
stress_output = MLS_book(stress_data, xy_data, xy_output, geometry);

%% Plot results
% If you want to save/plot magnitudes you can use this for example
% plot_results(problem_folder, 'test_plot', xy_output, stress_output(:,1),'-dpng', geometry)

%% Evaluate exact solution
% You can query the analytical stress solution of the test problem with:
% stress_exact = exact_solu(xy_output);