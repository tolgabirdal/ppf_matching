function q = q_getRandom()
%Q_GETRANDOM Uniformly distributed unit quaternion.
%   Q = Q_GETRANDOM() returns a random unit quaternion (1x4), uniformely
%   distributed on the 3-sphere.

%   Author: Damien Teney

% Reference: see (1)-(4) in
% http://mathworld.wolfram.com/HyperspherePointPicking.html.

x(1:4) = [+inf +inf +inf +inf];

while x(1) * x(1) + x(2) * x(2) >= 1 || x(3) * x(3) + x(4) * x(4) >= 1
  for i = 1:4
    x(i) = getRandom_uniform(-1, 1);
  end
end

tmp = (1 - x(1) * x(1) - x(2) * x(2)) / (x(3) * x(3) + x(4) * x(4));

q(1) = x(4) * sqrt(tmp);
q(2) = x(1);
q(3) = x(2);
q(4) = x(3) * sqrt(tmp);

% Old method below
%{
x1 = rand();
x2 = rand() * 2 * pi;
x3 = rand() * 2 * pi;

q = [sin(x2) * sqrt(1 - x1) ...
     cos(x2) * sqrt(1 - x1) ...
     sin(x3) * sqrt(x1) ...
     cos(x3) * sqrt(x1)];
%}