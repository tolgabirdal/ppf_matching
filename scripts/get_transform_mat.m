
% Given the parameters as [\theta_x, \theta_y, \theta_z, t_x, t_y, t_z]
% Returns the pose matrix [R | t]
% Author: Tolga Birdal
%
function M=get_transform_mat(par)

r=par(1:3);
t=par(4:6);
Rx=[1 0 0 ;
    0 cos(r(1)) -sin(r(1)) ;
    0 sin(r(1)) cos(r(1)) ];

Ry=[cos(r(2)) 0 sin(r(2)) ;
    0 1 0;
    -sin(r(2)) 0 cos(r(2))];

Rz=[cos(r(3)) -sin(r(3)) 0;
    sin(r(3)) cos(r(3)) 0;
    0 0 1];

M=[Rx*Ry*Rz t'];

end
