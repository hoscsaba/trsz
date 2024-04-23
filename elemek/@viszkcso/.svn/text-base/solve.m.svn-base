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
epsz2 = viszkcso.epsz2;

for i=1:N+1
    if isnan(v(i))
        fprintf('\n----------------------------------');
        disp(v);
        error('PÁNIK!!!!');
    end
end

% frissites a belso pontokban

for i = 2:N

    %%%%% Atalakitott valtozat %%%%%

    xL = x(i-1); xR = x(i+1);
    hL = h(i-1); hR = h(i+1); hP = h(i);
    
    vL = v(i-1); vR = v(i+1);    
    pL = p(i-1); pR = p(i+1);    
    aL = a(i-1); aR = a(i+1);
    
    dL = d(i-1); dR = d(i+1);
    epszL = epsz(i-1); epszR = epsz(i+1);
    epsz2L = epsz2(i-1); epsz2R = epsz2(i+1);
    
    tP = (xR - xL - t*(vR - aR)...
        + t*(vL + aL))/((vL + aL) - (vR - aR));
    dt = tP - t;
    xP = xL + (vL + aL)*dt;    
    
    alpL = d0/s0*(2*epszL+1);
    alpR = d0/s0*(2*epszR+1);
     
    gammaL = 2*epszL+1;
    gammaR = 2*epszR+1;
    
    GL = -2*E2/eta2*epsz2L*exp(-E2/eta2*dt) + alpL/eta2*(pL - p0)*exp(-E2/eta2*dt);
    GR = -2*E2/eta2*epsz2R*exp(-E2/eta2*dt) + alpR/eta2*(pR - p0)*exp(-E2/eta2*dt);
    
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
    epsz2P = epsz2(i)*exp(-E2/eta2*dt) + 0.5*alp/E2*(p(i) - p0)*(1 - exp(-E2/eta2*dt));
    epszP = epsz2P + 0.5*alp*(pP - p0)/E1;
    alpP = d0/s0*(2*epszP+1);
    %aP2 = sqrt(E1*(2*epszP+1)/(ro*alpP));
    %Kt = pP/((pP-p(i))/dt);
    
    aP = sqrt( (2*epszP+1)/(ro*alpP*(1/E1+1/E2*(1-exp(-E2/eta2*dt)))) );    
 
    dP = d0*(1 + epszP);
    
    vuj(i) = vP; puj(i) = pP; auj(i) = aP;
    epsz2uj(i) = epsz2P; epszuj(i) = epszP; duj(i) = dP;
    tuj(i) = tP; xuj(i) = xP;
    
end

% A peremfeltetelek elott ki kell szamolni a eleje/vege idoket!

tuj(1) = t - x(2)/(v(2) - a(2));
tuj(N+1) = t + (L - x(N))/(v(N) + a(N));
t_min = min(tuj);

% Peremfeltetelek beallitasa.

xuj(1) = 0;
tuj(1) = t - x(2)/(v(2) - a(2));
dt1 = tuj(1) - t;
dt = t_min - t;
x_interp = dt/dt1*dx;

v_2 = interp1(x,v,x_interp,'linear');
p_2 = interp1(x,p,x_interp,'linear');
a_2 = interp1(x,a,x_interp,'linear');
epsz2_2 = interp1(x,epsz2,x_interp,'linear');
epsz_2 = interp1(x,epsz,x_interp,'linear');
d_2 = interp1(x,d,x_interp,'linear');

hP = 0; hR = 0; % Ezt meg be kell epiteni

%alpR = d0/s0*(2*epsz(2)+1);
%GR = -2*E2/eta2*epsz2(2)*exp(-E2/eta2*dt) + alpR/eta2*(p(2) - p0)*exp(-E2/eta2*dt);
%JR = -dt*( g*(hP-hR)/(xP-xR) +32*nu/(d(2)^2)*v(2) + a(2)/(2*epsz(2)+1)*GR );
alpR = d0/s0*(2*epsz_2+1);
GR = -2*E2/eta2*epsz2_2*exp(-E2/eta2*dt) + alpR/eta2*(p_2 - p0)*exp(-E2/eta2*dt);
JR = -dt*( g*(hP-hR)/(xP-xR) +32*nu/(d_2^2)*v_2 - a_2/(2*epsz_2+1)*GR );

switch bct1  % Peremfeltétel a csõ elején
    case 'p'
        puj(1) = bcv1;
        %vuj(1) = v(2) + JR + 1/(ro*a(2))*puj(1) - 1/(ro*a(2))*p(2);
        vuj(1) = v_2 + JR + 1/(ro*a_2)*puj(1) - 1/(ro*a_2)*p_2;       
    case 'v'
        vuj(1) = bcv1;
        %puj(1) = ro*a(2)*( vuj(1) - v(2) - JR ) + p(2);
        puj(1) = ro*a_2*( vuj(1) - v_2 - JR ) + p_2;
    otherwise
        error('Ismeretlen peremfeltétel');
end

alp = d0/s0*(2*epsz(1)+1);
epsz2uj(1) = epsz2(1)*exp(-E2/eta2*dt) + 0.5*alp/E2*(p(1) - p0)*(1 - exp(-E2/eta2*dt));
epszuj(1) = epsz2uj(1) + 0.5*alp*(puj(1) - p0)/E1;
alpP = d0/s0*(2*epszuj(1)+1);
%auj(1) = sqrt(E1*(2*epszuj(1)+1)/(ro*alpP));
auj(1) = sqrt( (2*epszuj(1)+1)/(ro*alpP*(1/E1+1/E2*(1-exp(-E2/eta2*dt)))) ); 
duj(1) = d0*(1 + epszuj(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xuj(N+1) = L;
tuj(N+1) = t + (L - x(N))/(v(N) + a(N));
dtN1 = tuj(N+1) - t;
dt = t_min - t;
x_interp = L - dt/dtN1*dx;

v_N = interp1(x,v,x_interp,'linear');
p_N = interp1(x,p,x_interp,'linear');
a_N = interp1(x,a,x_interp,'linear');
epsz2_N = interp1(x,epsz2,x_interp,'linear');
epsz_N = interp1(x,epsz,x_interp,'linear');
d_N = interp1(x,d,x_interp,'linear');

hP = 0; hR = 0; % Ezt meg be kell epiteni

alpL = d0/s0*(2*epsz_N+1);
GL = -2*E2/eta2*epsz2_N*exp(-E2/eta2*dt) + alpL/eta2*(p_N - p0)*exp(-E2/eta2*dt);
JL = -dt*( g*(hP-hR)/(xP-xR) +32*nu/(d_N^2)*v_N + a_N/(2*epsz_N+1)*GL );

switch bct2  % Peremfeltétel a csõ végén
    case 'p'
        puj(N+1) = bcv2;
        %vuj(N+1) = JL + v(N) - 1/(ro*a(N))*puj(N+1) + 1/(ro*a(N))*p(N);
        vuj(N+1) = JL + v_N - 1/(ro*a_N)*puj(N+1) + 1/(ro*a_N)*p_N;
    case 'v'
        vuj(N+1) = bcv2;
        %puj(N+1) = ro*a(N)*( v(N) - vuj(N+1) + JL ) + p(N);
        puj(N+1) = ro*a_N*( v_N - vuj(N+1) + JL ) + p_N;
    otherwise
        error('Ismeretlen peremfeltétel');
end

alp = d0/s0*(2*epsz(N+1)+1);
epsz2uj(N+1) = epsz2(N+1)*exp(-E2/eta2*dt) + 0.5*alp/E2*(p(N+1) - p0)*(1 - exp(-E2/eta2*dt));
epszuj(N+1) = epsz2uj(N+1) + 0.5*alp*(puj(N+1) - p0)/E1;
alpP = d0/s0*(2*epszuj(N+1)+1);
%auj(N+1) = sqrt(E1*(2*epszuj(N+1)+1)/(ro*alpP));
auj(N+1) = sqrt( (2*epszuj(N+1)+1)/(ro*alpP*(1/E1+1/E2*(1-exp(-E2/eta2*dt)))) ); 

duj(N+1) = d0*(1 + epszuj(N+1));

% Interpoláció végrehajtása

t_min = min(tuj(:));

data.t = t; data.tuj = tuj;
data.x = x; data.xuj = xuj;
data.p = p; data.puj = puj;
data.v = v; data.vuj = vuj;
data.epsz2 = epsz2; data.epsz2uj = epsz2uj;
data.epsz = epsz; data.epszuj = epszuj;
data.d = d; data.duj = duj;
data.a = a; data.auj = auj;

out = intpol(t_min,data);

if max(isnan(puj)) == 1
    viszkcso.tranziens_agelem_2csp.nev
    bcv1
    puj
    error('Leallok, mert elegem van')
end

tuj =[];
tuj = t_min;

puj(2:N) = out.p(2:N);
vuj(2:N) = out.v(2:N);
epsz2uj(2:N) = out.epsz2(2:N);
epszuj(2:N) = out.epsz(2:N);
duj(2:N) = out.d(2:N);
auj(2:N) = out.a(2:N);

dt = tuj - t;

delta_d = duj - d;
delta_epsz = epszuj-epsz;

d_seb = delta_d/dt;
epsz_seb = delta_epsz/dt;


% Peremfeltetelek

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
%viszkcso.t = viszkcso.t+viszkcso.dt;

% A viszkoelasztikus paraméterek frissítése:

viszkcso.a = auj;
viszkcso.epsz = epszuj;
viszkcso.epsz2 = epsz2uj;
viszkcso.d = duj;
viszkcso.A = duj.^2*pi/4;
viszkcso.dseb = d_seb;

viszkcso.dt = dt;

viszkcso.all = {viszkcso.x.'; viszkcso.v.'; viszkcso.a.'; ...
                viszkcso.t; viszkcso.N; viszkcso.L};
           
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
            
rajz = 0;

if rajz == 1
    figure(1)
    
%     plot(x,duj/2,'+')
%     grid on
%     hold on
%     plot(x,-duj/2,'+')
%     plot([x(1) x(end)],[d0/2 d0/2],'r')
%     plot([x(1) x(end)],[-d0/2 -d0/2],'r')
%     hold off

    subplot(3,1,1)
    plot(x,puj,'r','LineWidth',2)
    title('Nyomas')
    grid on

    subplot(3,1,2)
    plot(x,vuj,'g','LineWidth',2)
    title('Sebesseg')
    YLim([0 0.4]);
    XLim([0 1.0]);
    grid on
    
    subplot(3,1,3)
    plot(x,duj,'b','LineWidth',2)
    title('Atmero')
    grid on
    pause
end

function out = intpol(tmin,inp)

t = inp.t; tuj = inp.tuj;
x = inp.x; xuj = inp.xuj;

p = inp.p; puj = inp.puj;
v = inp.v; vuj = inp.vuj;
epsz2 = inp.epsz2; epsz2uj = inp.epsz2uj;
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

