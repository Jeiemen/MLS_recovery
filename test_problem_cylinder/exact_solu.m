function sigma_ex = exact_solu(xy)
% thick-walled cylinder under internal pressure. Analytical solution for
% 2D plane strain hypothesis

% Problem parameters
a = 5; % internal radius
b = 20; % external radius
P = 1; % Internal pressure

% input coordinates
x = xy(1,:);
y = xy(2,:);
r = (x.^2 + y.^2).^0.5;
theta=atan2(y,x);

% Stress solution in Polar coordinates
c=b/a;
Sr = (P/(c^2-1)) * (1-(b./r).^2);
St = (P/(c^2-1)) * (1+(b./r).^2);

% Stress rotation to Cartesian coordinates
s_xx = Sr .* cos(theta).^2 + St .* sin(theta).^2;
s_yy = Sr .* sin(theta).^2 + St .* cos(theta).^2;
s_xy = (Sr - St) .* sin(theta) .* cos(theta);
sigma_ex = [s_xx ; s_yy ; s_xy];
