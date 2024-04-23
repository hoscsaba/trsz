function vena = solve(vena,pf)

%fprintf('\nA VENA solve függvény meghívása...\n')
%pause

bct1=pf{1}{1};
bcv1=pf{1}{2};
bct2=pf{2}{1};
bcv2=pf{2}{2};

%bct1 = 'v';
%bcv1 = 0;

vena.v1 = vena.v(1);
vena.p1 = vena.p(1);
vena.a1 = vena.a(1);
vena.d1 = vena.d(1);
vena.vN = vena.v(vena.N+1);
vena.pN = vena.p(vena.N+1);
vena.aN = vena.a(vena.N+1);
vena.dN = vena.d(vena.N+1);

%fprintf('\n\n Csõ lépés:  PFE tipus:  %s  érték: %g    PFV tipus:  %s  érték: %g',bct1,bcv1,bct2,bcv2);

%nu = cso.tranziens_agelem_2csp.nu;
ro = vena.tranziens_agelem_2csp.ro;

s  = vena.s;
a  = vena.a;
N  = vena.N;
g  = 9.81;
nu = vena.nu;
p  = vena.p;
v  = vena.v;
dzdx = vena.dzdx;
dt = vena.dt;
h = vena.h;
dx = vena.dx;
t = vena.t;
x = vena.x;
L = vena.L;
d = vena.d;
d0 = vena.d0;

E1 = vena.E1;
E2 = vena.E2;
eta2 = vena.eta2;

alp = vena.alp;
bet = vena.bet;
gam = vena.gam;
epsz = vena.epsz;
epsz2 = vena.epsz2;

for i=1:N+1
    if isnan(v(i))
        fprintf('\n----------------------------------');
        disp(v);
        error('PÁNIK!!!!');
    end
end

% frissites a belso pontokban

for i = 2:N

    % A pontos "t" és "x" kiszámítása:

    t_p(i) = (x(i+1) - x(i-1) - t*(v(i+1) - a(i+1))...
        + t*(v(i-1) + a(i-1)))/((v(i-1) + a(i-1)) - (v(i+1) - a(i+1)));

    delta_t = t_p(i) - t;

    x_p(i) = x(i-1) + (v(i-1) + a(i-1))*(t_p(i)-t);

    p_e = p(i)-1e5;
    
%     alphauj(i) = D/s*(2*epsz(i)+1);
%     gammauj(i) = 2*epsz2(i)*exp(-E2/eta2*delta_t) + alphauj(i)*1/E2*p_e...
%         *(1-exp(-E2/eta2*delta_t)) + 1/E1*alphauj(i)*p_e + 1;
%     auj(i) = sqrt((gammauj(i)*E1)/(ro*alphauj(i)));
%     betauj(i) = sqrt(E1/(ro*gammauj(i)*alphauj(i)));

    jobb_1 = a(i-1)*bet(i-1)*alp(i-1)/eta2*(p(i-1)-1e5) - a(i+1)*bet(i+1)*alp(i+1)/eta2*(p(i+1)-1e5);
    jobb_2 = a(i-1)*bet(i-1)*E2/eta2*epsz2(i-1) - a(i+1)*bet(i+1)*E2/eta2*epsz2(i+1);

    vuj(i) = 1/(ro*(a(i-1)+a(i+1)))*(p(i-1)-p(i+1) + ro*(a(i-1)*v(i-1)+a(i+1)*v(i+1)) -ro*delta_t...
        *(32*nu/d(i-1)^2*v(i-1)*a(i-1)+32*nu/d(i+1)^2*v(i+1)*a(i+1) + (jobb_1 -2*jobb_2)*exp(-E2/eta2*delta_t)));

    puj(i) = p(i-1) - ro*a(i-1)*(vuj(i)-v(i-1)) - ro*a(i-1)*delta_t*(32*nu/d(i-1)^2*v(i-1) + bet(i-1)...
        *alp(i-1)/eta2*(p(i-1)-1e5) - 2*bet(i-1)*E2/eta2*epsz2(i-1))*exp(-E2/eta2*delta_t);
    
    epsz2_uj(i) = epsz2(i)*exp(-E2/eta2*delta_t) + 0.5*alp(i)*1/E2*p_e*(1-exp(-E2/eta2*delta_t));
    epszuj(i) = epsz2_uj(i) + 0.5*1/E1*alp(i)*(puj(i)-1e5);

    alphauj(i) = d0/s*(2*epszuj(i)+1);
    gammauj(i) = 2*epsz2_uj(i)*exp(-E2/eta2*delta_t) + alphauj(i)*1/E2*p_e...
        *(1-exp(-E2/eta2*delta_t)) + 1/E1*alphauj(i)*(puj(i)-1e5) + 1;
    auj(i) = sqrt((gammauj(i)*E1)/(ro*alphauj(i)));
    betauj(i) = sqrt(E1/(ro*gammauj(i)*alphauj(i)));
    
    duj(i) = epszuj(i)*d0 + d0;    % Fontos: ide mindenkeppen d0-t kell irni...
    

    %gammauj(i) = 2*epsz2_e(i)*exp(-E2/eta2*delta_t) + alphauj(i)*1/E2*p_e...
    %    *(1-exp(-E2/eta2*delta_t)) + 1/E1*alphauj(i)*(puj(i)-1e5) + 1;
    %auj(i) = sqrt((gammauj(i)*E1)/(ro*alphauj(i)));
    
    %             epsz_2(pt,tt) = epsz_2_e*exp(-E_2/eta_2*delta_t) + 0.5*alpha*1/E_2*p_e*(1-exp(-E_2/eta_2*delta_t));
    %             epsz(pt,tt) = epsz_2(pt,tt) + 0.5*1/E_1*alpha*p(pt,tt);
    %
    %             d(pt,tt) = epsz(pt,tt)*D_0 + D_0;
    %
    %             % Új hangsebesség számítása
    %
    %             A = 2*epsz_2_e*exp(-E_2/eta_2*delta_t) + alpha*1/E_2*p_e...
    %                 *(1-exp(-E_2/eta_2*delta_t)) + 1/E_1*alpha*p(pt,tt) + 1;
    %             a(pt,tt) = sqrt((A*E_1)/(rho*alpha));
    %
    %    fprintf('\n    i=%d  dt*(AA-BB)=%+6.4e  ar-br=%+6.4e   vr=%+6.4e   vu=%+6.4e   pr=%+6.4e   pu=%+6.4e',i,dt*(AA-BB),ar-br,v(i),vuj(i),p(i),puj(i));
end

% Peremfeltetelek beallitasa.

%val1 = ro*vena.a*v - p;
%val2 = - ro*vena.a*dt*(32*vena.lambda/vena.D^2*v - vena.beta*vena.alpha/vena.eta2*...
%       p*exp(-vena.E2/vena.eta2*dt) + 2*vena.beta*vena.E2/vena.eta2*epsz2_e*exp(-vena.E2/vena.eta2*dt));
%val = val1 + val2;

%perem_1= p(2)-ro*a*v(2) + dt*( ro*a*g*dzdx(1) + ro*a*lambda/2/D*v(2)*abs(v(2)));
%perem_N= p(N)+ro*a*v(N) - dt*( ro*a*g*dzdx(N) + ro*a*lambda/2/D*v(N)*abs(v(N)));

% perem_1 = ro*a(2)*v(2)-p(2) - ro*auj(2)*delta_t*(32*lambda/D^2*v(2) - betauj(2)*alphauj(2)/eta2*...
%     (p(2)-1e5)*exp(-E2/eta2*delta_t) + 2*betauj(2)*E2/eta2*epsz2_e(2)*exp(-E2/eta2*delta_t));

% perem_N = ro*a(N)*v(N)+p(N) - ro*auj(N)*delta_t*(32*lambda/D^2*v(N) + betauj(N)*alphauj(N)/eta2*...
%     (p(N)-1e5)*exp(-E2/eta2*delta_t) - 2*betauj(N)*E2/eta2*epsz2_e(N)*exp(-E2/eta2*delta_t));

x_p(1) = 0;
t_p(1) = t - x(2)/(v(2) - a(2));

delta_t = t_p(1) - t;

p_e = p(1)-1e5;
% alphauj(1) = D/s*(2*epsz(1)+1);
% 
% gammauj(1) = 2*epsz2(1)*exp(-E2/eta2*delta_t) + alphauj(1)*1/E2*p_e...
%     *(1-exp(-E2/eta2*delta_t)) + 1/E1*alphauj(1)*(p(1)-1e5) + 1;
% auj(1) = sqrt((gammauj(1)*E1)/(ro*alphauj(1)));
% betauj(1) = sqrt(E1/(ro*gammauj(1)*alphauj(1)));

perem_1 = ro*a(2)*v(2)-p(2) - ro*a(2)*delta_t*(32*nu/d(2)^2*v(2) - bet(2)*alp(2)/eta2*...
    (p(2)-1e5)*exp(-E2/eta2*delta_t) + 2*bet(2)*E2/eta2*epsz2(2)*exp(-E2/eta2*delta_t));

switch bct1  % Peremfeltétel a csõ elején
    case 'p'
        puj(1) = bcv1;
        vuj(1) = (puj(1) + perem_1)/(ro*a(2));
    case 'v'
        vuj(1) = bcv1;
        puj(1) = ro*a(2)*vuj(1) - perem_1;
    otherwise
        error('Ismeretlen peremfeltétel');
end

epsz2_uj(1) = epsz2(1)*exp(-E2/eta2*delta_t) + 0.5*alp(1)*1/E2*p_e*(1-exp(-E2/eta2*delta_t));
epszuj(1) = epsz2_uj(1) + 0.5*1/E1*alp(1)*(puj(1)-1e5);

alphauj(1) = d0/s*(2*epszuj(1)+1);

gammauj(1) = 2*epsz2(1)*exp(-E2/eta2*delta_t) + alphauj(1)*1/E2*p_e...
    *(1-exp(-E2/eta2*delta_t)) + 1/E1*alphauj(1)*(puj(1)-1e5) + 1;
auj(1) = sqrt((gammauj(1)*E1)/(ro*alphauj(1)));
betauj(1) = sqrt(E1/(ro*gammauj(1)*alphauj(1)));

duj(1) = epszuj(1)*d0 + d0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_p(N+1) = L;
t_p(N+1) = t + (L - x(N))/(v(N) + a(N));

p_e = p(N+1)-1e5;

perem_N = ro*a(N)*v(N)+p(N) - ro*a(N)*delta_t*(32*nu/d(N)^2*v(N) + bet(N)*alp(N)/eta2*...
    (p(N)-1e5)*exp(-E2/eta2*delta_t) - 2*bet(N)*E2/eta2*epsz2(N)*exp(-E2/eta2*delta_t));

switch bct2  % Peremfeltétel a csõ végén
    case 'p'
        puj(N+1) = bcv2;
        vuj(N+1) = (perem_N - puj(N+1))/(ro*a(N));
    case 'v'
        vuj(N+1) = bcv2;
        puj(N+1) = perem_N-ro*a(N)*vuj(N+1);
    otherwise
        error('Ismeretlen peremfeltétel');
end

epsz2_uj(N+1) = epsz2(N+1)*exp(-E2/eta2*delta_t) + 0.5*alp(N+1)*1/E2*p_e*(1-exp(-E2/eta2*delta_t));
epszuj(N+1) = epsz2_uj(N+1) + 0.5*1/E1*alp(N+1)*(puj(N+1)-1e5);

alphauj(N+1) = d0/s*(2*epszuj(N+1)+1);

gammauj(N+1) = 2*epsz2(N+1)*exp(-E2/eta2*delta_t) + alphauj(N+1)*1/E2*p_e...
    *(1-exp(-E2/eta2*delta_t)) + 1/E1*alphauj(N+1)*(puj(N+1)-1e5) + 1;
auj(N+1) = sqrt((gammauj(N+1)*E1)/(ro*alphauj(N+1)));
betauj(N+1) = sqrt(E1/(ro*gammauj(N+1)*alphauj(N+1)));

duj(N+1) = epszuj(N+1)*d0 + d0;

% Interpoláció végrehajtása

t_min = min(t_p(:));

t_sz(:,1) = t_p(:);
t_sz(:,2) = t(:);

x_sz(:,1) = x_p(:);
x_sz(:,2) = x(:);

tomb(:,1,1) = p(:);
tomb(:,2,1) = puj(:);
tomb(:,1,2) = v(:);
tomb(:,2,2) = vuj(:);
tomb(:,1,3) = epsz2(:);
tomb(:,2,3) = epsz2_uj(:);
tomb(:,1,4) = epsz(:);
tomb(:,2,4) = epszuj(:);
tomb(:,1,5) = d(:);
tomb(:,2,5) = duj(:);
tomb(:,1,6) = a(:);
tomb(:,2,6) = auj(:);

[tomb_int] = visco_interpol(t_sz,x_sz,tomb,N);

t_uj = t_min;
puj(:) = tomb_int(:,1);
vuj(:) = tomb_int(:,2);
epsz2_uj(:) = tomb_int(:,3);
epszuj(:) = tomb_int(:,4);
duj(:) = tomb_int(:,5);
auj(:) = tomb_int(:,6);

dt = t_uj - t;

delta_d = duj - d;
delta_epsz = epszuj-epsz;

d_seb = delta_d/dt;
epsz_seb = delta_epsz/dt;

% Kavitáció ellenõrzése

pgoz = - 427.66*(ro/1000)^3 + 1309.4*(ro/1000)^2 - 1336.4*ro/1000 + 454.61;

for i=1:N+1
    if puj(i)<pgoz
        puj(i)=pgoz;
        if (vena.warnings==1),
            fprintf('\n!! Figyelmeztetés: gõzkiválás a %s csõben az %d-edik osztáspontban!\n',vena.tranziens_agelem_2csp.nev,i);
            fprintf('a gõz nyomása: %d',puj(i));
        end
    end
end

vena.p = puj;
vena.v = vuj;
vena.t = vena.t+dt;

% A viszkoelasztikus paraméterek frissítése:

vena.a = auj;
%vena.alp = alphauj; %%%%%%%%% !!!
vena.bet = betauj;
vena.gam = gammauj;
vena.epsz = epszuj;
vena.epsz2 = epsz2_uj;
vena.d = duj;
vena.A = duj.^2*pi/4;

vena.dseb = d_seb;


% vena.alp1 = alp(1);
% vena.bet1 = bet(1);
% vena.gam1 = gam(1);
% vena.epsz2_1 = epsz2_uj(1);
% 
% vena.alpN = alp(N+1);
% vena.betN = bet(N+1);
% vena.gamN = gam(N+1);
% vena.epsz2_N = epsz2_uj(N+1);

vena.dt = dt;

function [tomb_int] = visco_interpol(t_sz,x_sz,tomb,max_ii)

%%% Fontos: az algoritmusban hiba van, ellenorizni kell !!!

% t_sz(:,j:j+1) - az új és az elõzõ idõlépéshez tartozó "t"
%                 értékek
% x_sz(:,j:j+1) - az új és az elõzõ idõlépéshez tartozó "x"
%                 értékek

t_min = min(t_sz(:,2));

% a t_min idõsíkkal való metszéspontok megállapítása

for i = 1:max_ii+1
    %if i == 1
        %ke = 1;
        %kv = 1;
    %elseif i == max_ii+1
        %ke = 2;
        %kv = 2;
    %else
        ke = 1;
        kv = 2;
    %end
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
            %if (a*b > 0 && hely_b ~= 1)
            if (a*b > 0)
                dist = [dist(1:(hely_b-1)) dist((hely_b+1):end)];
            end
            %if hely_b == 1, break;
        end
        d = abs(a)+abs(b);
        L_a = abs(a)/d;
        L_b = abs(b)/d;
        tomb_int(i,1:6) = tomb_mp(hely_a,1:6)*L_b + tomb_mp(hely_b,1:6)*L_a;
    end
end


p_int = 1;
v_int = 1;


