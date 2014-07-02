function x = getRandom_uniform(min, max)
%GETRANDOM_UNIFORM Random numbers in uniform distribution.
%   X = GETRANDOM_UNIFORM(MIN, MAX) returns one random number in the
%   [min,max] interval.

%   Author: Damien Teney

x = min + rand(1) * (max - min);