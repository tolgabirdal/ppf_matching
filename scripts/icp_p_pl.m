
% by Tolga Birdal
% A combined implementation of many state of the art ICP ideas
% First, as described by .... the points are normalized. 
% Next, robust weighting as in .... is employed. 
% Bi-unique correspondences are obtained : 
% Finally, a linear approximation of point-to-plane ICP is used to compute
% the pose.
% This function can be used in coarse to fine matching in order to optimize
% the speed

function [Pose]=icp_p_pl(SrcPC, SrcN, DstPC, DstN)

TolP=0.025;
Iterations=100;
useRobustRejection = 0;
outlierScale = 3;
numNeighbors = 1; % Allow for 5 multiple assignments
useCov = 0;

n = length(SrcPC);

% shift to origin
meanSrc = mean(SrcPC);
SrcPC = bsxfun(@minus, SrcPC, meanSrc);
meanDst = mean(DstPC);
DstPC = bsxfun(@minus, DstPC, meanDst);

% compute average dist from origin
avgDist = (sqrt(sum(SrcPC.^2, 2)) + sqrt(sum(DstPC.^2, 2)))*0.5;
avgDist = sum(avgDist(:));

% scale to unit sphere
scale = n / avgDist;
SrcPC = SrcPC.*scale;
DstPC = DstPC.*scale;

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

while( (~(fval_perc<(1+TolP) && fval_perc>(1-TolP))) && i<Iterations)
    
    % Calculate closest point for all points
    [j,d]=knnsearch(kdtreeobj,Src_Moved,'k',numNeighbors);
    
    % Implement Picky ICP
    newI = (1:length(j))';
    newJ = j;
    rejectInd=[];
    
    % Step 1 of Picky ICP : Robustly reject outliers
    if (useRobustRejection>0)
        % I don't like this as it is very dependent on the dataset
        dMin = d(:,1);
        med = median(dMin);
        sigma = 1.4826 * mad(dMin,1); % Robust estimate of stddev
        threshold = (outlierScale*sigma+med);
        rejectInd = (dMin>threshold);
        rejectInd = find(rejectInd(:));
        j(rejectInd, :)=NaN;
        newI(rejectInd)=NaN;
        newJ=j;
        d(rejectInd)=NaN;
        
        newJ = newJ(~isnan(newJ(:,1)), :);
        newI = newI(~isnan(newI));
        d = d(~isnan(newJ(:,1)), :);
        
        mapSrc = newI;
        mapDst = newJ;
    end
        
    
    if (numNeighbors>1)
        % Resolve bi-unique correspondences (BC-ICP)
        n = length(newJ);
        mapSrc = NaN*ones(length(newJ), 1);
        mapDst = NaN*ones(length(newJ), 1);
        map = NaN*ones(length(newJ), 1);
        dists = Inf*ones(length(newJ), 1);
        
        for it=1:n
            jc = newJ(it,:);
            %curDist = dists(it);
            dc = d(it,:);  
                        
            % Sort the distance between pi and every point in q in ascending order;
            % Note that due to knnsearch, they are sorted already
            % dc = d(it,:);            
            % [sd, ind] = sort(dc);
            % jc = jc(ind);
            
            %mapSrc(it)=jc(1);
            %mapDst(it)=it;
            
          %  assigned =0 ;
            for k=1:numNeighbors
                jck = jc(k);
                dck = dc(k);
                curDist = dists(jck);
                
                aj = map(jck);
                if (isnan(aj))
                    if (dck < curDist)
                        mapDst(it) = it;
                        mapSrc(it) = jck;
                        map(jck)=1;
                        dists(jck) = dck;
                        break;
                    end
              %      assigned = 1;
                end
            end
            
            % no bi-unique correspondences are found
%             if (~assigned)
%                 mapDst(i)=NaN;
%             end
        end
        
        %numNeighbors=numNeighbors-1;
        
    else
        % Step 2 of Picky ICP:
        % Among the resulting corresponding pairs, if more than one scene point p_i
        % is assigned to the same model point m_j, then select p_i that corresponds
        % to the minimum distance.
        
        [duplicates, index]=get_duplicates(newJ);
        
        for di=1:length(duplicates)
            dup = duplicates(di);
            
            % Such search could be done much faster i.e. when implemented in C
            % Using, say, a hashtable
            indJ = (find(newJ==dup));
            
            dists = d(indJ);
            [minD, indD] = min(dists);
            
            newJ(indJ) = NaN;
            newI(indJ) = NaN;
            newJ( indJ(1) ) = indJ(indD);
            newI( indJ(1) ) = dup;
        end
            
        mapSrc = newI;
        mapDst = newJ;
    end
    
    mapSrc = mapSrc(~isnan(mapSrc));
    mapDst = mapDst(~isnan(mapDst));
        
    Src_Match=SrcPC(mapSrc,:);
    Dst_Match=DstPC(mapDst,:);
    DstN_Match=DstN(mapDst,:);
    
%     fval = sum( (Src_Match - Dst_Match).^2, 2);
%     outliers = (fval>1.4826 * mad(fval,1));
%     Src_Match=Src_Match(outliers,:);
%     Dst_Match=Dst_Match(outliers,:);
%     DstN_Match=DstN_Match(outliers,:);
%   
    if (useCov)
        x = minimize_point_to_plane_cov(Src_Match, Dst_Match, DstN_Match);
    else
        x = minimize_point_to_plane(Src_Match, Dst_Match, DstN_Match);
    end
    
    % Make the transformation matrix
    Pose=get_transform_mat(x);
    
    % Transform the Points
    Src_Moved=movepoints(Pose, SrcPC);
    
    fval = sum( (Src_Match - Dst_Match).^2, 2);
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
%     for t=1:4:length(Src_Match)
%         hold on, plot3([Src_Match(t,1) Dst_Match(t,1)], [Src_Match(t,2) Dst_Match(t,2)], [Src_Match(t,3) Dst_Match(t,3)],'g-');
%     end
    hold off;
    axis equal;
    %pause(1);
    drawnow;
   % pause;
    
    i=i+1;
end

%Pose = minPose;
Pose(1:3, 4) = Pose(1:3, 4)./scale + meanDst' - Pose(1:3, 1:3)*meanSrc';

end

function [duplicates, index]=get_duplicates(X)
    uniqueX = unique(X);
    countOfX = hist(X,uniqueX);
    index = (countOfX~=1);
    duplicates = uniqueX(index);
end

% Minimize the point to plane metric according to :
% Kok Lim Low : Linear Least Squares Optimization for Point-to-Plane
% ICP Surface Registration
function [x]=minimize_point_to_plane(Src, Dst, Normals)

%b1 = dot(Dst, Normals, 2);
%b2 = dot(Src, Normals, 2);
%b=b1-b2;
b = dot(Dst-Src, Normals, 2);
A1 = cross(Src, Normals);
A2 = Normals;
A=[A1 A2];

x = (A\b)';

end

% Gelfand et. al. 2003
function [x]=minimize_point_to_plane_cov(Src, Dst, Normals)

b = dot(Dst-Src, Normals, 2);
A1 = cross(Src, Normals);
A2 = Normals;
A=[A1 A2];

C = (A'*A);
b = A'*b;

x = (C\b)';

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


