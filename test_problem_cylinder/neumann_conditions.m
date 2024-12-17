function [traction, c_type] = neumann_conditions(xy, bound_id)
% Inputs:
% - xy [n_point x 2]: query points for imposed tractions
% - bound_id: boundary side ID
%
% Outputs:
% - traction [n_point x 2]: normal and tangent enforced tractions
% - c_type: cell with constraint type for each point (options: "traction" or "symmetry")


% thick-walled cylinder problem: 
% - Pressure enforce at internal arc segment (id - 2)
% - Free boundary at external arc segment (id - 4) 
% - Symmetry conditions at straight segments (id - 1, 3)
n_point = size(xy, 1);
traction = zeros(n_point, 2);
c_type = categorical(nan(n_point, 1));
c_type = addcats(c_type,{'symmetry', 'traction'});
for i = 1:n_point
    switch bound_id(i)
        case {1, 3}
            traction(i, :) = [0, 0];
            c_type(i) = 'symmetry';
        case 4
            traction(i, :) = [0, 0];
            c_type(i) = 'symmetry';
        case 2
            traction(i, :) = [-1, 0];
            c_type(i) = 'traction';
        otherwise
            error('Boundary id %d not defined for current geometry', bound_id(i))
    end
end