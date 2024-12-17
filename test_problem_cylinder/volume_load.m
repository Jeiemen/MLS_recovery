function b = volume_load(xy, E, nu)
% Inputs:
% - xy [n_point x 2]: query points for imposed tractions
% - E: Young's modulus
% - nu: Poisson's coefficient
%
% Outputs:
% - b [n_point x 2]: volume forces at given point


% thick-walled cylinder problem: 
% - No volume loads

if size(xy,1)~=2
    xy = xy';
end

NumOfPoints=size(xy,2);
b=[];
x=xy(1,:); y=xy(2,:);
b=[0;0]*ones(1,NumOfPoints);
