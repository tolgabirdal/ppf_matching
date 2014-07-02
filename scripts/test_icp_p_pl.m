
function []=test_icp_p_pl()

close all;

% Generate a random transformation
% q = q_getRandom();
% t = rand(1,3);
% R = quat2dcm(q);
% [R1 R2 R3] = quat2angle( q );
t=[0,0,0];
R1 = deg2rad(15);
R2 = deg2rad(0);
R3 = deg2rad(0);
x=[R1 R2 R3 t];
PoseGT = get_transform_mat(x);

PoseGT

% Prepare point cloud
[vertices, faces] = read_ply('C:\Users\tolga\Documents\GitHub\ppf_matching\data\parasaurolophus_6700_2.ply');

SrcPC = vertices;

% apply transformation
vertices_trans=movepoints(PoseGT, vertices);
normals = compute_normal(vertices_trans, faces);

DstPC = vertices_trans;
DstN = normals';

figure,plot3(SrcPC(:,1), SrcPC(:,2), SrcPC(:,3),'r.');
hold on,plot3(DstPC(:,1), DstPC(:,2), DstPC(:,3),'b.');
axis equal; pause;

FinalPose = icp_p_pl(SrcPC, DstPC, DstN)

% Display the final pose

%figure,plot3(SrcPC(:,1), SrcPC(:,2), SrcPC(:,3),'r.');
%hold on,plot3(DstPC(:,1), DstPC(:,2), DstPC(:,3),'b.');
%axis equal; %pause;

end

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
