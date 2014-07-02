
%SrcPC = [SrcPC ones(length(SrcPC),1)];
%DstPC = [DstPC ones(length(DstPC),1)];
%DstN = [DstN zeros(length(DstN),1)];

n = length(SrcPC);

meanSrc = mean(SrcPC);
SrcPC = bsxfun(@minus, SrcPC, meanSrc);

meanDst = mean(DstPC);
DstPC = bsxfun(@minus, DstPC, meanDst);

accum = sqrt(sum(SrcPC.^2, 2)) + sqrt(sum(DstPC.^2, 2));
accum = sum(accum(:));

factor = 2 * (n) / accum;

%SrcPC = SrcPC./factor;
%DstPC = DstPC./factor;