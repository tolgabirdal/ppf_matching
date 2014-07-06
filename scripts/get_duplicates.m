
% Returns the duplicates of a vector
% 
% Author: Tolga Birdal

function duplicates = get_duplicates(X)
uniqueX = unique(X);
countOfX = hist(X,uniqueX);
index = (countOfX~=1);
duplicates = uniqueX(index);
end