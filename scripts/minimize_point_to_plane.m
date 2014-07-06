
% Minimize the point to plane metric according to
% Kok Lim Low : Linear Least Squares Optimization for Point-to-Plane
% ICP Surface Registration
% Also check 
% Ef?cient Variants of the ICP Algorithm by Szymon Rusinkiewicz
% 
% Author: Tolga Birdal
%
function [x]=minimize_point_to_plane(Src, Dst, Normals)

b = dot(Dst-Src, Normals, 2);
A1 = cross(Src, Normals);
A2 = Normals;
A=[A1 A2];
x = (A\b)';

% A variant of Gelfand et. al. 2003 would use
% but we stick to the original one for now
% Checkout sample_pc_stable for more details
% C = (A'*A);
% b = A'*b;
% x = (C\b)';


end
