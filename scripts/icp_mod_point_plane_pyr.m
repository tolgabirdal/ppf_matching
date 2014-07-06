
%% icp_mod_point_plane_pyr
% by Tolga Birdal
% This file is an attempt towards fast and robust ICP. Here I combine many
% of the proposed variants into a single implementation to make ICP a 
% practical tool for matching, pose refinement, SLAM etc.
% The the details are explained in more detailed in ICP.pdf, but let me 
% briefly explain the variants I use in here.
% 1. A smart sampling option is provided (Gelfand et. al.). To use it set
% SampleType=1. However note that this function is provided seperately at:
% http://www.mathworks.com/matlabcentral/fileexchange/47138-stable-sampling-of-point-clouds-for-icp-registration
% 2. Robust outlier filtering is incorporated
% 3. Linear point-to-plane metric is minimized
% 4. A correspondence filtering as in picky icp was utilized
% 5. Multi-Resolution scheme is utilized (Coarse-to-Fine ICP)
% Morover, a Hartley-Zissermann type of normalization is done.
% Here is a sample usage:
% FinalPose = icp_mod_point_plane_pyr(SrcPC, SrcN, DstPC, DstN, 0.05, 100, 3, 1, 8, 0, 1);
% Paraameters:
%   SrcPC :         Nx3 array of input point cloud model
%   SrcN :          Nx3 array of the normals of the model
%   DstPC :         Nx3 array of the target point cloud (scene)
%   DstN :          Nx3 array of the normals of the target
%   Tolerence:      If there is insufficient change between iterations, ICP is
%                   terminated at this level. This insufficiency is determined 
%                   by Tolerence.
%   MaxIterations:  Maximum number of iterations. Note that this parameter is
%                   updated through the pyramid. Less iterations are
%                   carried out in the coarser levels
%   RejectionScale: Threshold for outlier rejection. RejectionScale*stddev
%                   is set as the threshold of rejection. A typically used
%                   value is 3.
%   NumNeighborsCorr: Not used for now
%   NumLevels:      Maximum levels of the pyramid (depth of the search)
%   SampleType:     0 for uniform sampling, 1 for stable sampling as in 
%                   Gelfand et. al.
%   DrawRegistration: 1 for visualization, 0 for no visualizations
%   Author: Tolga Birdal
%%
function [Pose]=icp_mod_point_plane_pyr(SrcPC, SrcN, DstPC, DstN, Tolerence, MaxIterations, RejectionScale, NumNeighborsCorr, NumLevels, SampleType, DrawRegistration)

% Assigne variables. Leave as default if not provided
TolP = 0.0001;
Iterations = 100;
useRobustRejection = 0;
outlierScale = 3;
%numNeighbors = 1; % Allow for N multiple assignments
visualize = 1;
nK = 1;
numPyr = 8;
samplePCType = 0;
%

% Assign parameters
if (nargin > 4), TolP = Tolerence;                  end
if (nargin > 5), Iterations = MaxIterations;        end
if (nargin > 6)
    useRobustRejection = (RejectionScale>0);
    outlierScale = RejectionScale;
end
if (nargin > 7)
    numNeighbors = NumNeighborsCorr;
    
    if (numNeighbors>0)
        nK = numNeighbors;
    end
end
if (nargin > 8), numPyr = NumLevels;      end
if (nargin > 9), samplePCType = SampleType;      end
if (nargin > 10), visualize = DrawRegistration;      end

n = length(SrcPC);

% Hartley-Zissermann Scaling:
meanSrc = mean(SrcPC);
SrcPC = bsxfun(@minus, SrcPC, meanSrc);
meanDst = mean(DstPC);
DstPC = bsxfun(@minus, DstPC, meanDst);

% compute average dist from origin
avgDist = (sum(sqrt(sum(SrcPC.^2, 2))) + sum(sqrt(sum(DstPC.^2, 2))))*0.5;
%avgDist = sum(avgDist(:));

% scale to unit sphere
scale = n / avgDist;
SrcPC = SrcPC.*scale;
DstPC = DstPC.*scale;

SrcPCOrig = SrcPC;
DstPCOrig = DstPC;
SrcNOrig = SrcN;
DstNOrig = DstN;

Pose = [eye(3) zeros(3,1) ] ;
PoseInit = Pose;

if (visualize)
    figure;
    %set_plot_nice('Registration Result on Synthetic Data', 'X', 'Y', 1)
end
IterationsInit = Iterations;
TolPInit = TolP;

% Construct the KD-tree out of destination points (scene)
DstPC = DstPCOrig;
DstN = DstNOrig;
kdtreeobj = KDTreeSearcher(DstPC,'distance','euclidean');

% Start the registration
for level=numPyr-1:-1:0
    
    % Obtain the parameters in this pyramid level
    numSamples = uint32(n./(2^level));
    TolP = TolPInit.*(2^(2*level));
    Iterations = uint32(IterationsInit./(2^(level)));
    
    % Obtain the sampled point clouds for this level
    Pose=PoseInit;
    SrcPC = movepoints(Pose, SrcPCOrig);
    rotPose = Pose;
    rotPose(:,4)=[0,0,0]';
    SrcN = movepoints(rotPose, SrcNOrig);   
    
    if (samplePCType==1)
        [SrcPC, SrcN] = sample_pc_stable(SrcPC, SrcN, numSamples);
    else
        [SrcPC, SrcN] = sample_pc_uniform(SrcPC, SrcN, numSamples);
    end    
    
    % Distance error in last itteration
    fval_old=inf;
    
    % Change in distance error between two itterations
    fval_perc=0;
    
    % Array which contains the transformed points
    Src_Moved=SrcPC;
        
    %minPose=[]; % we will keep track of the pose with minimum residual
    fval_min=inf;
   
    i=0;
    
    % Start main loop of ICP
    while( (~(fval_perc<(1+TolP) && fval_perc>(1-TolP))) && i<Iterations)
        
        % Calculate closest point for all points
        [j,d]=knnsearch(kdtreeobj,Src_Moved,'k',nK);
        
        % Implement Picky ICP or BC-ICP
        newI = (1:length(j))';
        newJ = j;
        
        % Step 1 of ICP : Robustly reject outliers
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
        end
        
        % Step 2 of Picky ICP:
        % Among the resulting corresponding pairs, if more than one scene point p_i
        % is assigned to the same model point m_j, then select p_i that corresponds
        % to the minimum distance.
        
        duplicates=get_duplicates(newJ);
        
        for di=1:length(duplicates)
            dup = duplicates(di);
            
            % Such search could be done much faster i.e. when implemented in C
            % Using, say, a hashtable
            indJ = (find(newJ==dup));
            
            dists = d(indJ);
            [~, indD] = min(dists);
            
            tempI = newI( indJ(indD) );
            newJ(indJ) = NaN;
            newI(indJ) = NaN;
            newJ( indJ(1) ) = dup;
            newI( indJ(1) ) = tempI;
        end
        
        mapSrc = newI;
        mapDst = newJ;
        
        mapSrc = mapSrc(~isnan(mapSrc));
        mapDst = mapDst(~isnan(mapDst));
        
        Src_Match=SrcPC(mapSrc,:);
        Dst_Match=DstPC(mapDst,:);
        DstN_Match=DstN(mapDst,:);
        
        x = minimize_point_to_plane(Src_Match, Dst_Match, DstN_Match);
        
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
            %minPose = Pose;
        end
        
        if (visualize)
            plot3(Src_Moved(:,1), Src_Moved(:,2), Src_Moved(:,3),'r.');
            hold on, plot3(DstPC(:,1), DstPC(:,2), DstPC(:,3),'b.');
            srcMap = Src_Moved(mapSrc, :);
            for t=1:2:length(srcMap)
                hold on, plot3([srcMap(t,1) Dst_Match(t,1)], [srcMap(t,2) Dst_Match(t,2)], [srcMap(t,3) Dst_Match(t,3)],'g-');
            end
            hold off;
            axis equal;
            %pause(1);
            drawnow;
            % pause;
        end
        
        % release some memory, just in case
        clear Src_Match;
        clear Dst_Match;
        clear DstN_Match;
        
        i=i+1;
    end
    
    PoseInit = [Pose; 0 0 0 1]*[PoseInit; 0 0 0 1];
    PoseInit = PoseInit(1:3, 1:4);
end

Pose = PoseInit;
%Pose = minPose;
Pose(1:3, 4) = Pose(1:3, 4)./scale + meanDst' - Pose(1:3, 1:3)*meanSrc';

end


