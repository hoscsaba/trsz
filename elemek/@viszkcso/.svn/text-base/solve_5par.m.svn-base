function viszkcso = solve(viszkcso,pf)

%fprintf('\nA viszkcso solve függvény meghívása...\n')
%pause

bct1=pf{1}{1};
bcv1=pf{1}{2};
bct2=pf{2}{1};
bcv2=pf{2}{2};

viszkcso.v1 = viszkcso.v(1);
viszkcso.p1 = viszkcso.p(1);
viszkcso.a1 = viszkcso.a(1);
viszkcso.d1 = viszkcso.d(1);
viszkcso.vN = viszkcso.v(viszkcso.N+1);
viszkcso.pN = viszkcso.p(viszkcso.N+1);
viszkcso.aN = viszkcso.a(viszkcso.N+1);
viszkcso.dN = viszkcso.d(viszkcso.N+1);

%fprintf('\n\n Csõ lépés:  PFE tipus:  %s  érték: %g    PFV tipus:  %s  érték: %g',bct1,bcv1,bct2,bcv2);

%nu = cso.tranziens_agelem_2csp.nu;
ro = viszkcso.tranziens_agelem_2csp.ro;

p0 = viszkcso.p0;

s0  = viszkcso.s;
a  = viszkcso.a;
N  = viszkcso.N;
g  = 9.81;
nu = viszkcso.nu;
p  = viszkcso.p;
v  = viszkcso.v;
dzdx = viszkcso.dzdx;
dt = viszkcso.dt;
h = viszkcso.h;
dx = viszkcso.dx;
t = viszkcso.t;
x = viszkcso.x;
L = viszkcso.L;
d = viszkcso.d;
d0 = viszkcso.d0;

E1 = viszkcso.E1;
E2 = viszkcso.E2;
eta2 = viszkcso.eta2;

epsz = viszkcso.epsz;
epsz1 = viszkcso.epsz2;
epsz1v = viszkcso.epsz1v;

for i=1:N+1
    if isnan(v(i))
        fprintf('\n----------------------------------');
        disp(v);
        error('PÁNIK!!!!');
    end
end

% frissites a belso pontokban

% E3 = E1;
% 
% eta1 = 1e5;
% E1 = 1e5;
%     eta1 = -8.1e5;
%     E1 = -6.44e4;
    
% Az ot parameter beolvasasa filebol


par = textread('params.dat','%f',5);

E1 = par(1);
E2 = par(2);
E3 = par(3);
eta1 = par(4);
eta2 = par(5);

for i = 2:N

    %%%%% Atalakitott valtozat %%%%%

    xL = x(i-1); xR = x(i+1);
    hL = h(i-1); hR = h(i+1); hP = h(i);
    
    vL = v(i-1); vR = v(i+1);    
    pL = p(i-1); pR = p(i+1);    
    aL = a(i-1); aR = a(i+1);
    
    dL = d(i-1); dR = d(i+1);
    epszL = epsz(i-1); epszR = epsz(i+1);
    epsz1L = epsz1(i-1); epsz1R = epsz1(i+1);
    
    epsz1vL = epsz1v(i-1); epsz1vR = epsz1v(i+1);
    
    tP = (xR - xL - t*(vR - aR)...
        + t*(vL + aL))/((vL + aL) - (vR - aR));
    dt = tP - t;
    xP = xL + (vL + aL)*dt;    
    
    alpL = d0/s0*(2*epszL+1);
    alpR = d0/s0*(2*epszR+1);
     
    gammaL = 2*epszL+1;
    gammaR = 2*epszR+1;
    
    %GL = -2*E2/eta2*epsz2L*exp(-E2/eta2*dt) + alpL/eta2*(pL - p0)*exp(-E2/eta2*dt);
    %GR = -2*E2/eta2*epsz2R*exp(-E2/eta2*dt) + alpR/eta2*(pR - p0)*exp(-E2/eta2*dt);
    
    K1 = eta1*eta2/E2; K2 = (eta2*E1/E2+eta1); K3 = E1;
    m1 = (K2 + sqrt(K2^2-4*K1*K3)) / (2*K1);
    m2 = sqrt(K2^2-4*K1*K3) / K1;
    
    epsz1sL = epsz1L; %epsz1vL = 0;
    szig3L = 0.5*alpL*(pL-p0);
    
    GL = (2*E1*K1^2*m2)^-1*exp(-m1*dt)*( -eta1*(epsz1vL*(K2^2 - 2*K1*K3 + K2*K1*m2 + ...
         exp(m2*dt)*K2*K1 - exp(m2*dt)*(K2^2 - 2*K1*K3)) + ...
         (K2 - exp(m2*dt)*K2 + (1 + exp(m2*dt))*m2*K1)*(epsz1sL*K3 - szig3L)) + ...
         2*E1*K1*( epsz1vL*K2 + 2*epsz1sL*K3 + epsz1vL*m2*K1 + exp(m2*dt)*epsz1vL*m2*K1 - ...
         exp(m2*dt)*(epsz1vL*K2 + 2*epsz1sL*K3 - 2*szig3L) - 2*szig3L));
 
    epsz1sR = epsz1R; %epsz1vR = 0;
    szig3R = 0.5*alpR*(pR-p0);
    
    GR = (2*E1*K1^2*m2)^-1*exp(-m1*dt)*( -eta1*(epsz1vR*(K2^2 - 2*K1*K3 + K2*K1*m2 + ...
         exp(m2*dt)*K2*K1 - exp(m2*dt)*(K2^2 - 2*K1*K3)) + ...
         (K2 - exp(m2*dt)*K2 + (1 + exp(m2*dt))*m2*K1)*(epsz1sR*K3 - szig3R)) + ...
         2*E1*K1*( epsz1vR*K2 + 2*epsz1sR*K3 + epsz1vR*m2*K1 + exp(m2*dt)*epsz1vR*m2*K1 - ...
         exp(m2*dt)*(epsz1vR*K2 + 2*epsz1sR*K3 - 2*szig3R) - 2*szig3R));
    
    JL = -dt*( g*(hP-hL)/(xP-xL) +32*nu/(dL^2)*vL + aL/gammaL*GL );
    JR = -dt*( g*(hP-hR)/(xP-xR) +32*nu/(dR^2)*vR - aR/gammaR*GR );
    
    vP = 1/(ro*aR + ro*aL)*( ro*aR*vR + ro*aL*vL + pL - pR + ro*aL*JL + ro*aR*JR );
    
    pP = pL + ro*aL*vL - ro*aL*vP + ro*aL*JL;
    
%         vuj(i) = 1/(ro*(a(i-1)+a(i+1)))*(p(i-1)-p(i+1) + ro*(a(i-1)*v(i-1)+a(i+1)*v(i+1)) -ro*delta_t...
%         *(32*nu/d(i-1)^2*(v(i-1)*a(i-1)+v(i+1)*a(i+1)) + (jobb_1 -2*jobb_2)*exp(-E2/eta2*delta_t)));
% 
%     puj(i) = p(i-1) - ro*a(i-1)*(vuj(i)-v(i-1)) - ro*a(i-1)*delta_t*(32*nu/d(i-1)^2*v(i-1) + bet(i-1)...
%         *alp(i-1)/eta2*(p(i-1)-1e5) - 2*bet(i-1)*E2/eta2*epsz2(i-1))*exp(-E2/eta2*delta_t);
    
    alp = d0/s0*(2*epsz(i)+1);
    %epsz2P = epsz2(i)*exp(-E2/eta2*dt) + 0.5*alp/E2*(p(i) - p0)*(1 - exp(-E2/eta2*dt));
    
    epsz1s = epsz1(i); %epsz1v = epsz1v(i);
    szig3 = 0.5*alp*(pP-p0);
    
    epsz1P = (2*K3*K1*m2)^-1*exp(-m1*dt)*( epsz1s*K3*K1*m2 + exp(m2*dt)*epsz1s*K3*K1*m2 - ...
         K1*m2*szig3 - exp(m2*dt)*K1*m2*szig3 + 2*exp(m1*dt)*K1*m2*szig3 + ...
         (-1 + exp(m2*dt))*(2*epsz1v(i)*K1*K3 + epsz1s*K2*K3 - K2*szig3) );
    
    epsz1vP = (epsz1P-epsz1(i))/dt; 
     
%     (2*K3*K1*m2)^-1
%     - K1*m2*szig3
%     - exp(m2*dt)*K1*m2*szig3
%     + 2*exp(m1*dt)*K1*m2*szig3
%     + (-1 + exp(m2*dt))*(2*epsz1v*K1*K3 + epsz1s*K2*K3 - K2*szig3)
%     exp(-m1*dt)
%     epsz1P
%     pause
    
    epsz12P = (2*E1*K3*K1*m2)^-1*exp(-m1*dt)*( 2*E1*epsz1s*K3*K1*m2 + 2*exp(m2*dt)*E1*epsz1s*K3*K1*m2 + ...
         epsz1vP*eta1*K3*K1*m2 + exp(m2*dt)*epsz1vP*eta1*K3*K1*m2 - 2*E1*K1*m2*szig3 - ...
         2*exp(m2*dt)*E1*K1*m2*szig3 + 4*exp(m1*dt)*E1*K1*m2*szig3 + ...
         (-1 + exp(m2*dt))*( -eta1*K3*(epsz1vP*K2 + 2*epsz1s*K3 - 2*szig3) + 2*E1*(2*epsz1vP*K1*K3 + ...
         epsz1s*K2*K3 - K2*szig3) ));
    
    epszP = epsz12P + 0.5*alp*(pP - p0)/E3;
    alpP = d0/s0*(2*epszP+1);
    aP = sqrt(E3*(2*epszP+1)/(ro*alpP));
    %aP = sqrt(E1*s0/(ro*d0));
    %aP = a(i); %%%
    
    dP = d0*(1 + epszP);
    
    vuj(i) = vP; puj(i) = pP; auj(i) = aP;
    epsz1uj(i) = epsz1P; epszuj(i) = epszP; epsz1vuj(i) = epsz1vP;
    duj(i) = dP;
    tuj(i) = tP; xuj(i) = xP;
    
end

% Peremfeltetelek beallitasa.

xuj(1) = 0;
tuj(1) = t - x(2)/(v(2) - a(2));
dt = tuj(1) - t;

hP = 0; hR = 0; % Ezt meg be kell epiteni

alpR = d0/s0*(2*epsz(2)+1);
epsz1sR = epsz1(2); epsz1vR = epsz1v(2);
szig3R = 0.5*alpR*(p(2)-p0);
%GR = -2*E2/eta2*epsz2(2)*exp(-E2/eta2*dt) + alpR/eta2*(p(2) - p0)*exp(-E2/eta2*dt);
GR = (2*E1*K1^2*m2)^-1*exp(-m1*dt)*( -eta1*(epsz1vR*(K2^2 - 2*K1*K3 + K2*K1*m2 + ...
        exp(m2*dt)*K2*K1 - exp(m2*dt)*(K2^2 - 2*K1*K3)) + ...
        (K2 - exp(m2*dt)*K2 + (1 + exp(m2*dt))*m2*K1)*(epsz1sR*K3 - szig3R)) + ...
        2*E1*K1*( epsz1vR*K2 + 2*epsz1sR*K3 + epsz1vR*m2*K1 + exp(m2*dt)*epsz1vR*m2*K1 - ...
        exp(m2*dt)*(epsz1vR*K2 + 2*epsz1sR*K3 - 2*szig3R) - 2*szig3R));
JR = -dt*( g*(hP-hR)/(xP-xR) +32*nu/(d(2)^2)*v(2) + a(2)/(2*epsz(2)+1)*GR );

switch bct1  % Peremfeltétel a csõ elején
    case 'p'
        puj(1) = bcv1;
        vuj(1) = v(2) + JR + 1/(ro*a(2))*puj(1) - 1/(ro*a(2))*p(2);
    case 'v'
        vuj(1) = bcv1;
        puj(1) = ro*a(2)*( vuj(1) - v(2) - JR ) + p(2);
    otherwise
        error('Ismeretlen peremfeltétel');
end

alp = d0/s0*(2*epsz(1)+1);
%epsz2uj(1) = epsz2(1)*exp(-E2/eta2*dt) + 0.5*alp/E2*(p(1) - p0)*(1 - exp(-E2/eta2*dt));
epsz1s = epsz1(1); %epsz1v = 0;
szig3 = 0.5*alp*(puj(1)-p0);
epsz1uj(1) = (2*K3*K1*m2)^-1*exp(-m1*dt)*( epsz1s*K3*K1*m2 + exp(m2*dt)*epsz1s*K3*K1*m2 - ...
        K1*m2*szig3 - exp(m2*dt)*K1*m2*szig3 + 2*exp(m1*dt)*K1*m2*szig3 + ...
        (-1 + exp(m2*dt))*(2*epsz1v(1)*K1*K3 + epsz1s*K2*K3 - K2*szig3) );
    
epsz12uj(1) = (2*E1*K3*K1*m2)^-1*exp(-m1*dt)*( 2*E1*epsz1s*K3*K1*m2 + 2*exp(m2*dt)*E1*epsz1s*K3*K1*m2 + ...
        epsz1v(1)*eta1*K3*K1*m2 + exp(m2*dt)*epsz1v(1)*eta1*K3*K1*m2 - 2*E1*K1*m2*szig3 - ...
        2*exp(m2*dt)*E1*K1*m2*szig3 + 4*exp(m1*dt)*E1*K1*m2*szig3 + ...
        (-1 + exp(m2*dt))*( -eta1*K3*(epsz1v(1)*K2 + 2*epsz1s*K3 - 2*szig3) + 2*E1*(2*epsz1v(1)*K1*K3 + ...
        epsz1s*K2*K3 - K2*szig3) ));
    
epszuj(1) = epsz12uj(1) + 0.5*alp*(puj(1) - p0)/E3;
alpP = d0/s0*(2*epszuj(1)+1);
auj(1) = sqrt(E3*(2*epszuj(1)+1)/(ro*alpP));
%auj(1) = sqrt(E1*s0/(ro*d0));
duj(1) = d0*(1 + epszuj(1));
epsz1vuj(1) = (epsz1uj(1)-epsz1(1))/dt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xuj(N+1) = L;
tuj(N+1) = t + (L - x(N))/(v(N) + a(N));
dt = tuj(N+1) - t;

hP = 0; hR = 0; % Ezt meg be kell epiteni

alpL = d0/s0*(2*epsz(N)+1);
%GL = -2*E2/eta2*epsz2(N)*exp(-E2/eta2*dt) + alpL/eta2*(p(N) - p0)*exp(-E2/eta2*dt);
epsz1sL = epsz1(N); epsz1vL = epsz1v(N);
szig3L = 0.5*alpL*(p(N)-p0);
GL = (2*E1*K1^2*m2)^-1*exp(-m1*dt)*( -eta1*(epsz1vL*(K2^2 - 2*K1*K3 + K2*K1*m2 + ...
        exp(m2*dt)*K2*K1 - exp(m2*dt)*(K2^2 - 2*K1*K3)) + ...
        (K2 - exp(m2*dt)*K2 + (1 + exp(m2*dt))*m2*K1)*(epsz1sL*K3 - szig3L)) + ...
        2*E1*K1*( epsz1vL*K2 + 2*epsz1sL*K3 + epsz1vL*m2*K1 + exp(m2*dt)*epsz1vL*m2*K1 - ...
        exp(m2*dt)*(epsz1vL*K2 + 2*epsz1sL*K3 - 2*szig3L) - 2*szig3L));
JL = -dt*( g*(hP-hR)/(xP-xR) +32*nu/(d(N)^2)*v(N) + a(N)/(2*epsz(N)+1)*GL );

switch bct2  % Peremfeltétel a csõ végén
    case 'p'
        puj(N+1) = bcv2;
        %vuj(N+1) = (perem_N - puj(N+1))/(ro*a(N));
        vuj(N+1) = JL + v(N) - 1/(ro*a(N))*p(N+1) + 1/(ro*a(N))*p(N);
    case 'v'
        vuj(N+1) = bcv2;
        %puj(N+1) = perem_N-ro*a(N)*vuj(N+1);
        puj(N+1) = ro*a(N)*( v(N) - vuj(N+1) + JL ) + p(N);
    otherwise
        error('Ismeretlen peremfeltétel');
end

alp = d0/s0*(2*epsz(N+1)+1);
%epsz2uj(N+1) = epsz2(N+1)*exp(-E2/eta2*dt) + 0.5*alp/E2*(p(N+1) - p0)*(1 - exp(-E2/eta2*dt));

epsz1s = epsz1(N+1); %epsz1v = 0;
szig3 = 0.5*alp*(p(N+1)-p0);
epsz1uj(N+1) = (2*K3*K1*m2)^-1*exp(-m1*dt)*( epsz1s*K3*K1*m2 + exp(m2*dt)*epsz1s*K3*K1*m2 - ...
        K1*m2*szig3 - exp(m2*dt)*K1*m2*szig3 + 2*exp(m1*dt)*K1*m2*szig3 + ...
        (-1 + exp(m2*dt))*(2*epsz1v(N+1)*K1*K3 + epsz1s*K2*K3 - K2*szig3) );
epsz12uj(N+1) = (2*E1*K3*K1*m2)^-1*exp(-m1*dt)*( 2*E1*epsz1s*K3*K1*m2 + 2*exp(m2*dt)*E1*epsz1s*K3*K1*m2 + ...
    epsz1v(N+1)*eta1*K3*K1*m2 + exp(m2*dt)*epsz1v(N+1)*eta1*K3*K1*m2 - 2*E1*K1*m2*szig3 - ...
    2*exp(m2*dt)*E1*K1*m2*szig3 + 4*exp(m1*dt)*E1*K1*m2*szig3 + ...
    (-1 + exp(m2*dt))*( -eta1*K3*(epsz1v(N+1)*K2 + 2*epsz1s*K3 - 2*szig3) + 2*E1*(2*epsz1v(N+1)*K1*K3 + ...
    epsz1s*K2*K3 - K2*szig3) ));

epszuj(N+1) = epsz12uj(N+1) + 0.5*alp*(puj(N+1) - p0)/E3;
alpP = d0/s0*(2*epszuj(N+1)+1);
auj(N+1) = sqrt(E3*(2*epszuj(N+1)+1)/(ro*alpP));
%auj(N+1) = sqrt(E1*s0/(ro*d0));
duj(N+1) = d0*(1 + epszuj(N+1));
epsz1vuj(N+1) = (epsz1uj(N+1)-epsz1(N+1))/dt;

% Interpoláció végrehajtása

t_min = min(tuj(:));

data.t = t; data.tuj = tuj;
data.x = x; data.xuj = xuj;
data.p = p; data.puj = puj;
data.v = v; data.vuj = vuj;
data.epsz1 = epsz1; data.epsz1uj = epsz1uj;
data.epsz = epsz; data.epszuj = epszuj;
data.d = d; data.duj = duj;
data.a = a; data.auj = auj;

out = intpol(t_min,data);

tuj =[];
tuj = t_min;

% Uj interpolacio

puj = out.p;
vuj = out.v;
epsz1uj = out.epsz2;
epszuj = out.epsz;
duj = out.d;
auj = out.a;

dt = tuj - t;

delta_d = duj - d;
delta_epsz = epszuj-epsz;

d_seb = delta_d/dt;
epsz_seb = delta_epsz/dt;

% Kavitáció ellenõrzése

pgoz = - 427.66*(ro/1000)^3 + 1309.4*(ro/1000)^2 - 1336.4*ro/1000 + 454.61;

for i=1:N+1
    if puj(i)<pgoz
        puj(i)=pgoz;
        if (viszkcso.warnings==1),
            fprintf('\n!! Figyelmeztetés: gõzkiválás a %s csõben az %d-edik osztáspontban!\n',viszkcso.tranziens_agelem_2csp.nev,i);
            fprintf('a gõz nyomása: %d',puj(i));
        end
    end
end

viszkcso.p = puj;
viszkcso.v = vuj;
viszkcso.t = viszkcso.t+dt;

% A viszkoelasztikus paraméterek frissítése:

viszkcso.a = auj;
viszkcso.epsz = epszuj;
viszkcso.epsz2 = epsz1uj;
viszkcso.epsz1v = epsz1vuj;
viszkcso.d = duj;
viszkcso.A = duj.^2*pi/4;

viszkcso.dseb = d_seb;

viszkcso.dt = dt;

% FONTOS!!! UJ IDOLEPES KISZAMITASA!!!

dtuj(1) = x(2)/(-vuj(2) + auj(2));
dtuj(N+1) = (L - x(N))/(vuj(N) + auj(N));
for i = 2:N
     xL = x(i-1); xR = x(i+1);
     
     vL = vuj(i-1); vR = vuj(i+1);    
     aL = auj(i-1); aR = auj(i+1);

     tPuj = (xR - xL - viszkcso.t*(vR - aR)...
         + viszkcso.t*(vL + aL))/((vL + aL) - (vR - aR));
     dtuj(i) = tPuj - viszkcso.t;
end

viszkcso.dtuj = min(dtuj);

function [tomb_int] = visco_interpol(t_sz,x_sz,tomb,max_ii)

%%% Fontos: az algoritmusban hiba van, ellenorizni kell !!!

% t_sz(:,j:j+1) - az új és az elõzõ idõlépéshez tartozó "t"
%                 értékek
% x_sz(:,j:j+1) - az új és az elõzõ idõlépéshez tartozó "x"
%                 értékek

t_min = min(t_sz(:,2));

% a t_min idõsíkkal való metszéspontok megállapítása

for i = 1:max_ii+1
        ke = 1;
        kv = 2;

    % Egyenes egyenlet: r = r0 + v*tau

    for k = ke:kv

        % K+ irányban (páros index)

        if k == 1
            if i == max_ii+1
                vpt = i;
            else
                vpt = i+1;
            end
            mp = 2*i;
        elseif k == 2
            if i == 1
                vpt = i;
            else
                vpt = i-1;
            end
            mp = 2*i-1;
        end

        % K- irányban (páratlan index)

        vr(1) = x_sz(vpt,2)- x_sz(i,1);
        vr(2) = t_sz(vpt,2)- t_sz(i,1);

        tau = (t_min - t_sz(i,1))/vr(2);
        x_m(mp)= x_sz(i,1) + vr(1)*tau;

        % Interpolációs együtthatók

        s(1) = x_m(mp) - x_sz(i,1);
        s(2) = t_min - t_sz(i,1);

        v_absz = sqrt(vr(1)^2 + vr(2)^2);
        s_absz = sqrt(s(1)^2 + s(2)^2);

        L(1) = s_absz/v_absz;
        L(2) = 1 - L(1);

        for l = 1:6
            tomb_mp(mp,l) = L(1)*tomb(vpt,2,l) + L(2)*tomb(i,1,l);
        end
    end
end

% Interpoláció t_min mentén

for i = 1:max_ii+1
    dist = [];
    tt = 2*(max_ii+1);
    dist(1:tt) = x_sz(i,1)-x_m(1:tt);
    [a hely_a] = min(abs(dist));
    a = dist(hely_a);
    if abs(a) <= 1e-10
        tomb_int(i,1:6) = tomb_mp(hely_a,1:6);
    else

        dist = [dist(1:(hely_a-1)) dist((hely_a+1):tt)];
        b = a;
        while a*b > 0
            [b hely_b] = min(abs(dist));
            b = dist(hely_b);
            if (a*b > 0)
                dist = [dist(1:(hely_b-1)) dist((hely_b+1):end)];
            end
        end
        d = abs(a)+abs(b);
        L_a = abs(a)/d;
        L_b = abs(b)/d;
        tomb_int(i,1:6) = tomb_mp(hely_a,1:6)*L_b + tomb_mp(hely_b,1:6)*L_a;
    end
end

function out = intpol(tmin,inp)

t = inp.t; tuj = inp.tuj;
x = inp.x; xuj = inp.xuj;

p = inp.p; puj = inp.puj;
v = inp.v; vuj = inp.vuj;
epsz2 = inp.epsz1; epsz2uj = inp.epsz1uj;
epsz = inp.epsz; epszuj = inp.epszuj;
d = inp.d; duj = inp.duj;
a = inp.a; auj = inp.auj;

Nmax = length(tuj);

k = 1;
for i = 2:Nmax
    tP = tuj(i);
    xP = xuj(i); xL = x(i-1); pP = puj(i); pL = p(i-1); vP = vuj(i); vL = v(i-1);
    epsz2P = epsz2uj(i); epsz2L = epsz2(i-1); epszP = epszuj(i); epszL = epsz(i-1);
    dP = duj(i); dL = d(i-1); aP = auj(i); aL = a(i-1);
    
    xIL(k) = xP*(tmin - t)/(tP - t) + xL*(tP - tmin)/(tP - t);
    pIL(k) = pP*(tmin - t)/(tP - t) + pL*(tP - tmin)/(tP - t);
    vIL(k) = vP*(tmin - t)/(tP - t) + vL*(tP - tmin)/(tP - t);
    epsz2IL(k) = epsz2P*(tmin - t)/(tP - t) + epsz2L*(tP - tmin)/(tP - t);
    epszIL(k) = epszP*(tmin - t)/(tP - t) + epszL*(tP - tmin)/(tP - t);
    dIL(k) = dP*(tmin - t)/(tP - t) + dL*(tP - tmin)/(tP - t);
    aIL(k) = aP*(tmin - t)/(tP - t) + aL*(tP - tmin)/(tP - t);
    
    k = k+1;
end

k = 1;
for i = 1:Nmax-1
    
    tP = tuj(i);
    xP = xuj(i); xR = x(i+1); pP = puj(i); pR = p(i+1);  vP = vuj(i); vR = v(i+1);
    epsz2P = epsz2uj(i); epsz2R = epsz2(i+1); epszP = epszuj(i); epszR = epsz(i+1);
    dP = duj(i); dR = d(i+1); aP = auj(i); aR = a(i+1);
    
    temp = xP*(tmin - t)/(tP - t) + xR*(tP - tmin)/(tP - t);
    
    if ~ismember(temp,xIL)
        xIR(k) = xP*(tmin - t)/(tP - t) + xR*(tP - tmin)/(tP - t);   
        pIR(k) = pP*(tmin - t)/(tP - t) + pR*(tP - tmin)/(tP - t);
        vIR(k) = vP*(tmin - t)/(tP - t) + vR*(tP - tmin)/(tP - t);
        epsz2IR(k) = epsz2P*(tmin - t)/(tP - t) + epsz2R*(tP - tmin)/(tP - t);
        epszIR(k) = epszP*(tmin - t)/(tP - t) + epszR*(tP - tmin)/(tP - t);
        dIR(k) = dP*(tmin - t)/(tP - t) + dR*(tP - tmin)/(tP - t);
        aIR(k) = aP*(tmin - t)/(tP - t) + aR*(tP - tmin)/(tP - t);
        
        k = k+1;
    end
end

xILR = [xIL xIR];
pILR = [pIL pIR];
vILR = [vIL vIR];
epsz2ILR = [epsz2IL epsz2IR];
epszILR = [epszIL epszIR];
dILR = [dIL dIR];
aILR = [aIL aIR];

pI = interp1(xILR,pILR,x,'linear',NaN);
vI = interp1(xILR,vILR,x,'linear',NaN);
epsz2I = interp1(xILR,epsz2ILR,x,'linear',NaN);
epszI = interp1(xILR,epszILR,x,'linear',NaN);
dI = interp1(xILR,dILR,x,'linear',NaN);
aI = interp1(xILR,aILR,x,'linear',NaN);

% Perem a cso elejen
pI(1) = puj(1)*(tmin - t)/(tuj(1) - t) + p(1)*(tuj(1) - tmin)/(tuj(1) - t);
vI(1) = vuj(1)*(tmin - t)/(tuj(1) - t) + v(1)*(tuj(1) - tmin)/(tuj(1) - t);
epsz2I(1) = epsz2uj(1)*(tmin - t)/(tuj(1) - t) + epsz2(1)*(tuj(1) - tmin)/(tuj(1) - t);
epszI(1) = epszuj(1)*(tmin - t)/(tuj(1) - t) + epsz(1)*(tuj(1) - tmin)/(tuj(1) - t);
dI(1) = duj(1)*(tmin - t)/(tuj(1) - t) + d(1)*(tuj(1) - tmin)/(tuj(1) - t);
aI(1) = auj(1)*(tmin - t)/(tuj(1) - t) + a(1)*(tuj(1) - tmin)/(tuj(1) - t);

% Perem a cso vegen
pI(Nmax) = puj(Nmax)*(tmin - t)/(tuj(Nmax) - t) + p(Nmax)*(tuj(Nmax) - tmin)/(tuj(Nmax) - t);
vI(Nmax) = vuj(Nmax)*(tmin - t)/(tuj(Nmax) - t) + v(Nmax)*(tuj(Nmax) - tmin)/(tuj(Nmax) - t);
epsz2I(Nmax) = epsz2uj(Nmax)*(tmin - t)/(tuj(Nmax) - t) + epsz2(Nmax)*(tuj(Nmax) - tmin)/(tuj(Nmax) - t);
epszI(Nmax) = epszuj(Nmax)*(tmin - t)/(tuj(Nmax) - t) + epsz(Nmax)*(tuj(Nmax) - tmin)/(tuj(Nmax) - t);
dI(Nmax) = duj(Nmax)*(tmin - t)/(tuj(Nmax) - t) + d(Nmax)*(tuj(Nmax) - tmin)/(tuj(Nmax) - t);
aI(Nmax) = auj(Nmax)*(tmin - t)/(tuj(Nmax) - t) + a(Nmax)*(tuj(Nmax) - tmin)/(tuj(Nmax) - t);

out.p = pI; out.v = vI; out.epsz2 = epsz2I;
out.epsz = epszI; out.d = dI; out.a = aI;

