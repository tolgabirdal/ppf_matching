
% Test function for efficient ICP registration
% Author: Tolga Birdal
function []=test_icp_p_pl()

close all;

% Generate a tough random transformation
t=[0.1,1.2,1.1];
R1 = deg2rad(14);
R2 = deg2rad(37);
R3 = deg2rad(-31);
x=[R1 R2 R3 t];
PoseGT = get_transform_mat(x);

disp('Ground truth pose: ');
PoseGT

% Prepare point cloud
[vertices, faces] = read_ply('parasaurolophus_6700_2.ply');

SrcPC = vertices;
normals = compute_normal(vertices, faces);
SrcN = normals';

% apply transformation
vertices_trans=movepoints(PoseGT, vertices);
normals = compute_normal(vertices_trans, faces);

DstPC = vertices_trans;
DstN = normals';

figure,plot3(SrcPC(:,1), SrcPC(:,2), SrcPC(:,3),'r.');
hold on,plot3(DstPC(:,1), DstPC(:,2), DstPC(:,3),'b.');
axis equal; 
% pause; % Enable to pause before registration

% Corrput the Point Cloud:

% addd some noise to the scene
DstN = DstN + 0.025*rand(size(DstN));
DstPC = DstPC + 0.025*rand(size(DstPC));

% simulate a partial observation
DstN = DstN(1:length(DstN)/2,:);
DstPC = DstPC(1:length(DstPC)/2,:);

% FinalPose = icp_mod_point_plane(SrcPC, SrcN, DstPC, DstN, 0.001, 10000, -1, 1, 250, 1);
tic ();
%[Pose]=icp_mod_point_plane_pyr(SrcPC, SrcN, DstPC, DstN, Tolerence, MaxIterations, RejectionScale, NumNeighborsCorr, NumLevels, SampleType, DrawRegistration)
FinalPose = icp_mod_point_plane_pyr(SrcPC, SrcN, DstPC, DstN, 0.05, 100, 3, 1, 8, 0, 1);
toc();

disp('Estimated Pose: ');
FinalPose

% Display the final pose
Src_Moved=movepoints(FinalPose, SrcPC);

plot3(Src_Moved(:,1), Src_Moved(:,2), Src_Moved(:,3),'r.');
hold on,plot3(DstPC(:,1), DstPC(:,2), DstPC(:,3),'b.');
axis equal;

end
