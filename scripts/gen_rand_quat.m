
siz=1;

d   = randn( [4, prod( siz )] );
n   = sqrt( sum( d.^2, 1 ));
dn  = bsxfun( @rdivide, d, n );
neg = dn(1,:) < 0;
dn(:,neg) = -dn(:,neg);
q   = [dn(1,:), dn(2,:), dn(3,:), dn(4,:) ]
%q   = reshape( q, siz );