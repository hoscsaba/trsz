function yy = simi_spline(x,y,p,pf,xx)
%
% yy = simi_spline(x,y,p,pf,xx)
% közelítõ spline
% az (x,y) pontokra illeszt a p közelítõ paraméterrel
% pf - peremfeltétel
%      ha véges, akkor pf(1) = g'(x_1) és pf(2) = g'(x_N)
%      ha inf, akkor g''(x_1) = 0 illetve g''(x_N) = 0
% visszatérési értéke yy = g(xx)
%
% Vaik István, 2005.
% 
N = length(x);
if length(y) ~= N || length(p) ~= N
    exit;
end;

% x szerint növekvõ sorba rendezés
[x, s] = sort(x);
y = y(s);
p = p(s);


dx = x(2:N) - x(1:N-1);

A = diag( dx/6, -1) + diag( [ dx(1)/3 ; ( dx(1:N-2) + dx(2:N-1) )/3; dx(N-1)/3 ] ) + diag( dx/6,1);
B = -diag( 1./dx, -1 ) + diag([ 1/dx(1); 1./dx(1:N-2) + 1./dx(2:N-1); 1/dx(N-1) ] ) - diag( 1./dx, 1);

b = [-pf(1); zeros( N-2,1); -pf(2)];

if pf(1) == inf
    A(1,1) = 1;
    A(1,2) = 0;
    B(1,1) = 0;
    B(1,2) = 0;
    b(1) = 0;
end;

if pf(2) == inf
    A(N,N) = 1;
    A(N,N-1) = 0;
    B(N,N) = 0;
    B(N,N-1) = 0;
    b(N) = 0;
end;

D = diag(p);

gdd = (A+B/D*B)\(b-B*y);
g = y + D\B*gdd;



M = length(xx);
yy = zeros(M,1);
for k=1:M
    i = max( find( x < xx(k) ) );
    if isempty(i)
        if xx(k) == x(1)
            yy(k) = g(1);
        else
            continue;
        end;
    else
        if xx(k) > x(N)
            continue;
        else
            yy(k) = g(i)*( dx(i) - ( xx(k)-x(i) ) ) / dx(i);
            yy(k) = yy(k) + g(i+1)*( xx(k)-x(i) )/ dx(i);
            yy(k) = yy(k) + gdd(i)*( -(xx(k)-x(i))^3 + 3*(xx(k)-x(i))^2*dx(i) - 2*(xx(k)-x(i))*dx(i)^2 )/ 6 / dx(i);
            yy(k) = yy(k) + gdd(i+1)*( (xx(k)-x(i))^3 - (xx(k)-x(i))*dx(i)^2 )/6/dx(i);
        end;
    end;
end;


[xx, s] = sort(xx);
yy = yy(s);

% ialso = min( find( xx >= x(1) ) );
% ifelso = max( find( xx <= x(N) ) );
% 
% plot(xx(ialso:ifelso), yy(ialso:ifelso), '-b', x, y, 'xr', x, g, 'og');