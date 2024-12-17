function stress_out = MLS_book(stress_data, xy_data, xy_out, geometry)
% Inputs:
% - stress_data: input data with 2D stress field as a matrix [n_input x 3]
% - xy_data: input point cloud (n_input points) coordinates as a matrix [n_input x 2]
% - xy_out: point cloud (n_output points) coordinates as a matrix [n_output x 2]

% Outputs:
% - stress_out: recovered stress at output points as a matrix [n_output x 3]

%% Settings
% Constraint data
c_inter = true; % Enforce iternal equilibrium (bool)
c_bound = true; % Enforce boundary equilibrium (bool)
c_inter_full = 1; % Internal equilibrium scaling factor, value between [0, 1]

% check inputs
n_stress_comp = 3; % implemented only for 2D problems
assert(size(stress_data, 2) == n_stress_comp, 'Stress input must be a 3-column matrix')
assert(size(xy_data, 2) == 2, 'Input coordinates must be a 2-column matrix')
assert(size(xy_out, 2) == 2, 'Output coordinates must be a 2-column matrix')
n_input = size(xy_data, 1); % number of data points
n_output = size(xy_out, 1); % number of output points

% General configuration
recovery_degree = 3; % Degree of the recovered field
max_patch_points = 151; % Number of influence points for the MLS.
bound_tol = 1e-12; %Tolerance for boundary conditions 

% Evaluate the maximum influence radius considering a uniform GP
% distribution
x_min = min(xy_data(:,1)); x_max = max(xy_data(:,1));
y_min = min(xy_data(:,2)); y_max = max(xy_data(:,2));
r_domain = sqrt((x_max-x_min)^2 + (y_max-y_min)^2)/2;
A_domain = pi() * r_domain^2;
point_density = n_input / A_domain;
r_max_patch = sqrt(4 * max_patch_points / (point_density * pi()));

%% MSL loop
%Initilize output matrix
stress_out = zeros(n_output, n_stress_comp);
r_patch_all = zeros(n_output, 1);
patch_degree = inf(n_output, 1);
for i_out = 1:n_output
    xy = xy_out(i_out,:);
    
    % search closest data points to build MLS patch
    distance = sqrt(sum((xy_data - xy).^2, 2));
    if n_input <= max_patch_points
        target_points = 1:n_input;
    else
        inside_points = find(distance <= r_max_patch);
        if numel(inside_points) > max_patch_points % choose the closest points
            [~, aux_sort] = sort(distance(inside_points),'ascend');
            target_points = inside_points(aux_sort(1:max_patch_points));
        else
            target_points = inside_points;
        end
    end
    n_patch_pts = numel(target_points);
    r_patch = max(distance(target_points));
    r_patch_all(i_out) = r_patch;
    
    % Local coordinates of the patch points
    local_x = (xy_data(target_points,1) - xy(1))/r_patch;
    local_y = (xy_data(target_points,2) - xy(2))/r_patch;
    S = sqrt(local_x.^2 + local_y.^2);
    
    % Evaluate weighting function
    W = (1 + (-6*S + 8*S.^2 - 3*S.^3).*S);
    W(W<0) = 0; % needed due to round-off errors
    
    % Adjust degree of interpolation in patch for small data situation
    poli_len = [3,6,10,15,21,28]; % length of the 2D-polynomial basis
    patch_deg = recovery_degree;
    while n_patch_pts < poli_len(patch_deg) % at least one point per term
        patch_deg = patch_deg - 1;
        if patch_deg < 1
            error(['Not enough data points inside patch.' ...
                ' Consider increasing data resolution or patch size.']);
        end
    end
    patch_degree(i_out) = patch_deg;
    recovery_size = poli_len(patch_deg);
    total_size = recovery_size*n_stress_comp;
    
    % Create raw system of equations
    pTerms = polininterp2d(local_x, local_y, patch_deg);
    pTerms_W = bsxfun(@times, pTerms, W);
    M_single = pTerms' * pTerms_W;
    M = block_matrix(M_single, recovery_size, n_stress_comp);
    G_aux = pTerms_W'*stress_data(target_points,:);
    G = G_aux(:);

    % Internal equilibrium. Partial derivative terms
    if c_inter 
        % derivatives of weighting function
        dWds = (-12 + 24*S - 12*S.^2);
        dWdx = dWds.* local_x / r_patch^2;
        dWdy = dWds.* local_y / r_patch^2;
        % constraint matrix
        pTerms_dWdx = bsxfun(@times, pTerms, dWdx);
        dMdx_single = pTerms' * pTerms_dWdx;
        M = block_matrix(M_single, recovery_size, n_stress_comp);
        dMdx = block_matrix(dMdx_single, recovery_size, n_stress_comp);
        dGdx_aux = pTerms_dWdx'*stress_data(target_points,:);
        dGdx = dGdx_aux(:);
        % constraint vector
        pTerms_dWdy = bsxfun(@times, pTerms, dWdy);
        dMdy_single = pTerms' * pTerms_dWdy;
        dMdy = blkdiag(dMdy_single, dMdy_single, dMdy_single);
        dGdy_aux = pTerms_dWdy'*stress_data(target_points,:);
        dGdy = dGdy_aux(:);
    end

    % Inilialize variables in case Lagrange for boundary
    Lagrange = false;
    Madd = [];
    Cadd = [];
    if c_bound
        % Check if the patch support contains any boundary
        [boundary_xy, normal_vector, boundary_id] = ...
            find_boundary(xy, r_patch, geometry);
        n_boundary = numel(boundary_id);
        if n_boundary > 0
            % Local coordinates of the boundary points.
            XBP = (boundary_xy(:,1) - xy(1)) / r_patch;
            YBP = (boundary_xy(:,2) - xy(2)) / r_patch;
            % weighting function evaluated at boundary points (W')
            SBP = (XBP.^2 + YBP.^2).^0.5;
            % If output point is located at boundary, Lagrange enforcement
            Lagrange = SBP <= bound_tol;
            % If outpout point is located at a corner, only enforce one
            % boundary constraint as Lagrange
            if sum(Lagrange) == 2
                corner_points = find(Lagrange);
                [delete_id] = select_constrained_boundary(boundary_xy(corner_points,:), boundary_id(corner_points));
                delete_point = corner_points(delete_id);
                Lagrange(delete_point) = [];
                SBP(delete_point) = [];
                XBP(delete_point) = [];
                YBP(delete_point) = [];
                boundary_xy(delete_point, :) = [];
                normal_vector(delete_point, :) = [];
                boundary_id(delete_point) = [];
                n_boundary = n_boundary - 1;
            end
            WpFactor = 1;
            WBP = abs(1./SBP - 6*SBP + 8*SBP.^2 - 3*SBP.^3) * WpFactor;
            if c_inter
                dWPds = (-1./SBP.^3 - 6./SBP + 16 - 9*SBP) * WpFactor;
                dWBPdx = dWPds .* XBP / r_patch^2;
                dWBPdy = dWPds .* YBP / r_patch^2;
            end
            for i = 1:n_boundary
                %Prepare matrices
                CA = normal_vector(i, 1);
                SA = normal_vector(i, 2);
                RStress = [...
                    CA^2  , SA^2 ,  2*SA*CA   ;...
                    SA^2  , CA^2 , -2*SA*CA   ;...
                    -SA*CA, SA*CA,  CA^2-SA^2];
                pBP = polininterp2d(XBP(i), YBP(i), patch_deg, 1);
                PTotalNode = zeros(n_stress_comp, total_size);
                for i_comp = 1:n_stress_comp
                    Index = recovery_size * (i_comp - 1) + 1;
                    PTotalNode(i_comp, Index:Index+recovery_size-1) = pBP;
                end
                RPxx = RStress(1,:) * PTotalNode; 
                RPxy = RStress(3,:) * PTotalNode;
                % Read boundary constraint
                [traction, c_type] = neumann_conditions(boundary_xy(i,:), boundary_id(i));
                switch c_type
                    case 'traction'
                        PtRtRP   = RPxx'*RPxx   + RPxy'*RPxy;
                        PtRtSigt = RPxx'*traction(1) + RPxy'*traction(2);
                    case 'symmetry' % ==> no tangent traction = 0
                        PtRtRP = RPxy' * RPxy;
                    otherwise
                        error('Boundary node type %s is not defined.', c_type)
                end
                % Add equations to M and G
                if Lagrange(i)
                    if  c_type ~= 'symmetry'
                        Madd = [Madd ; RPxx ; RPxy];
                        Cadd = [Cadd ; traction(1) ; traction(2)];
                    else
                        Madd = [Madd ; RPxy];
                        Cadd = [Cadd ; 0];
                    end
                else
                    M = M + WBP(i)*PtRtRP;
                    if c_type ~= 'symmetry'
                        G = G + WBP(i)*PtRtSigt;
                    end
                    if c_inter
                        dMdx = dMdx + dWBPdx(i)*PtRtRP;
                        dMdy = dMdy + dWBPdy(i)*PtRtRP;
                        if c_type ~= 'symmetry'
                            dGdx = dGdx + dWBPdx(i)*PtRtSigt;
                            dGdy = dGdy + dWBPdy(i)*PtRtSigt;
                        end
                    end
                end
            end
        end
    end
    % You can draw patch support data with the following function
    % plot_patch_support(xy_patch, xy_out, r_patch, geometry, bound_id)
    
    % Internal equilibrium enforced at patch center with coords xy and
    % local coords [0, 0]
    if c_inter
        % Evaluate body forces for r.h.s
        b = volume_load(xy, [], []);
        % pTerms and its partial derivatives at origin
        [p0, dp0dx, dp0dy] = polininterp2d(0, 0, patch_deg, 1);
        P0 = blkdiag(p0, p0, p0);
        dP0dx = blkdiag(dp0dx, dp0dx, dp0dx);
        dP0dy = blkdiag(dp0dy, dp0dy, dp0dy);
        % constraint matrix and vector
        Minv = inv(M);
        T1x = dP0dx - c_inter_full*P0*Minv*dMdx;
        T1y = dP0dy - c_inter_full*P0*Minv*dMdy;
        MConstr=[T1x(1, :) + T1y(3, :);...
                 T1y(2, :) + T1x(3, :)];
        T2x = c_inter_full*P0*Minv*dGdx;
        T2y = c_inter_full*P0*Minv*dGdy;
        GConstr = [-T2x(1, :) - T2y(3,:) - b(1);...
                   -T2y(2, :) - T2x(3,:) - b(2)];
        % Assemble constraints in global system
        M = [M       , MConstr';...
             MConstr , zeros(size(MConstr,1), size(MConstr,1))];
        G = [G; GConstr];
    end

    % Solve the system
    if any(Lagrange)
        if size(Madd,2) < size(M,2)
            Madd(size(Madd,1),size(M,2)) = 0;
        end
        M = [M Madd';Madd sparse(size(Madd,1),size(Madd,1))];
        G = [G; Cadd];
    end
    solution = M \ G;
    % Gather only constant coefficient solutions for each component
    % (solution at local coord [0, 0])
    stress_out(i_out,:) = solution(1 : recovery_size : ...
        ((n_stress_comp - 1) * recovery_size + 1))';
end
end

function full_matrix = block_matrix(M, sizeM, repetitions)
total_size = sizeM * repetitions;
full_matrix = zeros(total_size);
for i_rep = 1:repetitions
    id = sizeM * (i_rep - 1) + 1;
    full_matrix(id : id + sizeM - 1, id : id + sizeM - 1) = M;
end
end
function [delete_id] = select_constrained_boundary(boundary_xy, boundary_id)
    % WARNING: not all cases are currently implemented. Below is a
    % tentative case list
    [~, c_type] = neumann_conditions(boundary_xy, boundary_id);
    n_traction_bounds = sum(c_type == 'traction');
    switch n_traction_bounds
        case 0
            % 2 symmetry surfaces with same normal = 1 symmetry surface (same)
            % 2 symmetry surfaces with different normal = enforced 0 traction
            delete_id = 2;
        case 1
        % 1 symmetry + 1 traction = 1 traction surface (symmetry is either redundant or impossible to accomplish)
            delete_id = c_type ~= 'traction';
        case 2
            % 1 traction + 1 zero traction = 1 traction (either redundant or impossible to accomplish)
            % 2 non-zero traction = first traction (Probably impossible to accomplish. Combination could be added)
            delete_id = 2;
        otherwise
            error('An output point belongs to more than 2 boundaries. This is not implemented.')
    end
end