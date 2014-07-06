
% by Tolga Birdal
% Uniform sampling of point clouds

function [SrcSample, SrcSampleNormals]=sample_pc_uniform(Src, NormalsSrc, numPoints)

n = length(Src);

% Really sample them.
sampledIndices=1:(n/numPoints):n;

SrcSample = Src(sampledIndices, :);
SrcSampleNormals = NormalsSrc(sampledIndices, :);

end
