function cso = solve(cso,pf)

bct1=pf{1}{1};
bcv1=pf{1}{2};
bct2=pf{2}{1};
bcv2=pf{2}{2};

%fprintf('\n\n Cs� l�p�s:  PFE tipus:  %s  �rt�k: %g    PFV tipus:  %s  �rt�k: %g',bct1,bcv1,bct2,bcv2);

%nu = cso.tranziens_agelem_2csp.nu; 
ro = cso.tranziens_agelem_2csp.ro;

D  = cso.D; a  = cso.a; 
N  = cso.N; g  = 9.81;  lambda = cso.lambda;
p  = cso.p; v  = cso.v;  dzdx=cso.dzdx;
dt = cso.dt; h=cso.h; dx=cso.dx;

for i=1:N+1
    if isnan(v(i))
        fprintf('\n----------------------------------');
        disp(v);
        error('PANIK!!!!');
    end
end

% frissites a belso pontokban
for i=2:N            
    AA = -ro*a*g*dzdx(i-1) - lambda/2/D*ro*a*v(i-1)*abs(v(i-1));
    BB =  ro*a*g*dzdx(i) + lambda/2/D*ro*a*v(i+1)*abs(v(i+1));
    ar = p(i-1)+ro*a*v(i-1); br = p(i+1)-ro*a*v(i+1);
    vuj(i) = (dt*(AA-BB)+ar-br)/2/ro/a;
    puj(i) = (dt*(AA+BB)+ar+br)/2;
%    fprintf('\n    i=%d  dt*(AA-BB)=%+6.4e  ar-br=%+6.4e   vr=%+6.4e   vu=%+6.4e   pr=%+6.4e   pu=%+6.4e',i,dt*(AA-BB),ar-br,v(i),vuj(i),p(i),puj(i));
end

% Peremfeltetelek beallitasa.
% A dzdx(1) �gy j� ahogy van!!!!!!
perem_1= p(2)-ro*a*v(2) + dt*( ro*a*g*dzdx(1) + ro*a*lambda/2/D*v(2)*abs(v(2)));
perem_N= p(N)+ro*a*v(N) - dt*( ro*a*g*dzdx(N) + ro*a*lambda/2/D*v(N)*abs(v(N)));

switch bct1
    case 'p'
        puj(1) = bcv1-cso.h(1)*ro*9.81*0;
        vuj(1) = (puj(1)-perem_1)/ro/a;
    case 'v'
        vuj(1) = bcv1;
        puj(1) = perem_1 + ro*a*vuj(1);
    otherwise
        error('Ismeretlen peremfeltetel');        
end

switch bct2
    case 'p'  
        puj(N+1) = bcv2-cso.h(end)*ro*9.81*0;
        vuj(N+1) = -(puj(N+1)-perem_N)/ro/a;                                        
%        fprintf('\n v(N)=%+6.4e  p(N)=%+6.4e  perem_N=%+6.4e   v(N)=%+6.4e  p(N)=%+6.4e',v(N),p(N),perem_N,vuj(N+1),puj(N+1));
    case 'v'        
        vuj(N+1) = bcv2;        
        puj(N+1) = perem_N-ro*a*vuj(N+1);
    otherwise
        error('Ismeretlen peremfeltetel');
end

% Kavitacio ellenorzese
% pgoz=-427.66*(ro/1000)^3 + 1309.4*(ro/1000)^2 - 1336.4*ro/1000 + 454.61;
% for i=1:N+1
%     if puj(i)<pgoz
%         puj(i)=pgoz;
%             warning(['!! Figyelmeztetes: gazkivalas a ',...
%                 cso.tranziens_agelem_2csp.nev,' csoben az ',num2str(i),...
%                 '-edik oszt�spontban!']);
%     end
% end

cso.p=puj; cso.v=vuj;
cso.v1=cso.v(1);       cso.p1=cso.p(1);
cso.vN=cso.v(cso.N+1); cso.pN=cso.p(cso.N+1);
cso.t=cso.t+dt;
