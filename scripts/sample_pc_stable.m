
% by Tolga Birdal
% This function implements the point sampling strategy from 
% Gelfand et. al. 2003
% The algorithm targets a point cloud sampling of the model for
% registration using the ICP algorithm. 
% The method is described in the paper
% Geometrically Stable Sampling for the ICP Algorithm, 3DIM 2003
% The idea is to sample the points by constraining the eigenvectors
% of the covariance matrix of the torque and force.
% The details and relation to the paper are given in the comments

function [SrcSample, SrcSampleNormals]=sample_pc_stable(Src, NormalsSrc, numPoints)

n = length(Src);

% Pre-processing as in Hartley Zissermann (map to unit sphere):

% shift to origin
meanSrc = mean(Src);
SrcPC = bsxfun(@minus, Src, meanSrc);

% compute average dist from origin
avgDist = sqrt(sum(SrcPC.^2, 2)) ;
avgDist = sum(avgDist(:));

% scale to unit sphere
scale = n / avgDist;
SrcPC = SrcPC.*scale;

% Figure out how bad the situation would be if we just used 
% random sampling. To do this, we estimate what C would be 
% using just a small uniform sampling of points.

nSample = max(numPoints / 10, 10);
indices = 1+floor((rand(n,1))*n);
indices = indices(1:nSample);

SrcRand = SrcPC(indices, :);
NormalsSrcRand = NormalsSrc(indices, :);

C = estimate_cov(SrcRand, NormalsSrcRand);

% sorted eigenvectors of cov matrix
[V, D]=eig(C);
[~, ind]=sort(diag(D), 'descend');
Vs = V(:, ind);
%kappa = d(1)/d(6);  % If we were to compute this, kappa would not be ~1.

L=cell(6);
Lind=cell(6);
totals = zeros(6,1);

% Step B1: form v
v= [ cross(SrcPC, NormalsSrc) NormalsSrc];

% Step B2: form Lk
for k=1:6
    xk = Vs(:, k)';
    Xk = repmat(xk, n, 1);
    mag = (dot(v, Xk,2));
    mag = mag.*mag;
    [~, ind] = sort(mag, 'descend');
    L{k} = v(ind, :);
    Lind{k} = ind;
end

% Step B3 and B4 : Constrain eigenvectors and choose the points

% keep track of point indices
lCounts=ones(6,1);
sampledIndices = zeros(numPoints*2, 1);

for i=1:n
    
    % Choose the next point from sorted list that has the smallest total
    [~, indT] = min (totals);
    cL = L{indT};
    vi = cL(lCounts(indT), :)';
    sampledIndices(i)=Lind{indT}(lCounts(indT));
    lCounts(indT) = lCounts(indT) + 1;
    
    % Update totals
    Vi = repmat(vi, 1, 6);
    err = dot (Vi, Vs);
    errMag = err.*err;
    totals = totals + errMag';
    
    kappa = totals(1)./totals(6);
    
    % if we sampled enough points and the condition number is
    % satisfactory (around 1)
    % Note that kappa constraint is usually satisfied and we end up
    % with desired number of points
    if (i>=numPoints && kappa>0.95 && kappa<1.15)
        numPoints=i;
        break;
    end
end

% Really sample them.
sampledIndices=sampledIndices(1:numPoints);
SrcSample = Src(sampledIndices, :);
SrcSampleNormals = NormalsSrc(sampledIndices, :);

end

% Gelfand et. al. 2003
function [C]=estimate_cov(Src, Normals)

A1 = cross(Src, Normals);
A2 = Normals;
A=[A1 A2];
C = (A'*A);

end
