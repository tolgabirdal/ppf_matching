
% by Tolga Birdal
% A sample test file for the demo of sample_pc_stable

[vertices, faces] = read_ply('parasaurolophus_6700_2.ply');

SrcPC = vertices;
normals = compute_normal(vertices, faces);
SrcN = normals';

% Perform stable sampling
[SrcSample, SrcSampleNormals]=sample_pc_stable(SrcPC, SrcN, 1200);

figure, 
subplot(1,2,1);
plot3(SrcPC(:,1), SrcPC(:,2), SrcPC(:,3), 'b.');
legend('Original Point Cloud','Location','NorthWest');
axis equal;
subplot(1,2,2);
hold on, plot3(SrcSample(:,1), SrcSample(:,2), SrcSample(:,3), 'r.');
legend('Sampled Points');
axis equal;