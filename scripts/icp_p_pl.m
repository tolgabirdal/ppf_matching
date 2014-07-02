% by Tolga Birdal
% Implementation of Kok Kim Low's ICP registration
% based on the linearization of the rotations.

function [Pose]=icp_p_pl(SrcPC, DstPC, DstN)

TolP=0.025;
Iterations=100;
    
n = length(SrcPC);

meanSrc = mean(SrcPC);
SrcPC = bsxfun(@minus, SrcPC, meanSrc);

meanDst = mean(DstPC);
DstPC = bsxfun(@minus, DstPC, meanDst);

trel = meanDst-meanSrc;

% Set initial rigid parameters
% 3 angles, 3 translations
% x=[0 0 0 0 0 0];

% Distance error in last itteration
fval_old=inf;

% Change in distance error between two itterations
fval_perc=0;
fval_min=0;

% Array which contains the transformed points
Src_Moved=SrcPC;
%SrcN_Moved=SrcN;

kdtreeobj = KDTreeSearcher(DstPC,'distance','euclidean');

minPose=[];
fval_min=inf;

figure;
i=0;
% perform trimmed icp
while( (~(fval_perc<(1+TolP) && fval_perc>(1-TolP))) && i<Iterations)
    
    % Calculate closest point for all points
    [j,d]=knnsearch(kdtreeobj,Src_Moved,'k',1);
    Dst_Match=DstPC(j,:);
    DstN_Match=DstN(j,:);
    
    x = minimize_point_to_plane(SrcPC, Dst_Match, DstN_Match);
    
    % Make the transformation matrix
    Pose=get_transform_mat(x);
    
    % Transform the Points
    Src_Moved=movepoints(Pose, SrcPC);
    
    fval = sum( (Src_Moved - Dst_Match).^2, 2);
    fval = sum(fval)./length(fval);
    
    % Calculate change in error between itterations
    fval_perc=fval/fval_old;
    
    % Store error value
    fval_old=fval;
    
    if (fval < fval_min)
        fval_min = fval;
        minPose = Pose;
    end
    
    plot3(Src_Moved(:,1), Src_Moved(:,2), Src_Moved(:,3),'r.');
    hold on, plot3(DstPC(:,1), DstPC(:,2), DstPC(:,3),'b.');
    hold off;
    axis equal;
    %pause(1);
    drawnow;
    
    i=i+1;
end

Pose = minPose;

Pose(1:3, 4) = Pose(1:3, 4)+trel;

end

% Minimize the point to plane metric according to :
% Kok Lim Low : Linear Least Squares Optimization for Point-to-Plane
% ICP Surface Registration
function [x]=minimize_point_to_plane(Src, Dst, Normals)

b1 = dot(Dst, Normals, 2);
b2 = dot(Src, Normals, 2);
b=b1-b2;

A1 = cross(Src, Normals);
A2 = Normals;
A=[A1 A2];

x = (A\b)';

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
