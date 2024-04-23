function out = karp(viszkcso,dt,dtmin)

%fprintf('A karpv függvény meghívása. "dt" értéke: %d \n',dt)
%pause


% FIGYELEM!!! Mivel a vénás szakasz végén jelentkeznek a problémák,
% könnyen lehet, hogy ebben a szubrutinban van a hiba!!
%disp(dt);

ro = viszkcso.tranziens_agelem_2csp.ro;

if dt < 0
    if abs(dt) < dtmin;
        dt = dtmin;
    end
end

if dt <= 0
    %if abs(dt)<dtmin;
    %    dt=abs(dt);
    %else
    %disp(dt);
    %error('Nagy a baj...');

    dt = abs(dt);
    delta_t = viszkcso.dt;

    N  = viszkcso.N;
    p_uj = viszkcso.p(N+1);
    v_uj = viszkcso.v(N+1);
    a_uj = viszkcso.a(N+1);
    dd_uj = viszkcso.d(N+1);
    p_e = viszkcso.pN;
    v_e = viszkcso.vN;
    a_e = viszkcso.aN;
    dd_e = viszkcso.dN;

    % Lineáris interpoláció (j és j-1 között)

    L1 = dt/delta_t;
    L2 = (delta_t-dt)/delta_t;
    p = L1*p_e + L2*p_uj;
    v = L1*v_e + L2*v_uj;
    a = L1*a_e + L2*a_uj;
    dd = L1*dd_e + L2*dd_uj;

    val1 = 0;
    val2 = 0;
    val = ro*a*v + p;

    A = dd^2*pi/4;

    %error('Csak semmi pánik...');

    %end
elseif dt > 0

    x = viszkcso.x;
    L = viszkcso.L;
    aa = viszkcso.a;
    vv = viszkcso.v;
    N = viszkcso.N;

    delta_t = (L - x(N))/(vv(N) + aa(N));

    % A peremen interpolálunk. (A csõ végén kapcsolódó csõ lép...)

    %x_interp = viszkcso.L-viszkcso.dx*dt/viszkcso.dt;
    %x_interp = viszkcso.L-viszkcso.dx*dt/viszkcso.ddtN;
    x_interp = viszkcso.L-viszkcso.dx*dt/delta_t;
    p = interp1(viszkcso.x,viszkcso.p,x_interp);
    v = interp1(viszkcso.x,viszkcso.v,x_interp);
    %dzdx=cso.dzdx(cso.N);
    %ro = viszkcso.tranziens_agelem_2csp.ro;

    %fprintf('x_interp: %5d\n',x_interp);
    %fprintf('viszkcso.L: %5d\n',viszkcso.L);

    if isnan(v)
        fprintf('\n\n')
        disp(viszkcso.x);
        disp(viszkcso.v);
        disp(cso.a*dt);
        error('PÁNIK!!!')
    end

    %val1=p+ro*cso.a*v;
    %val2= - dt*( ro*cso.a*9.81*(cso.h(cso.N+1)-cso.h(cso.N))/cso.dx + ro*cso.a*cso.lambda/2/cso.D*v*abs(v));
    %val=val1+val2;

    % alpha = viszkcso.alphaN;
    % beta = viszkcso.betaN;
    % gamma = viszkcso.gammaN;
    % a = viszkcso.aN;
    % epsz2_e = viszkcso.epsz2N;


    alp = interp1(viszkcso.x,viszkcso.alp,x_interp);
    bet = interp1(viszkcso.x,viszkcso.bet,x_interp);
    gam = interp1(viszkcso.x,viszkcso.gam,x_interp);
    a = interp1(viszkcso.x,viszkcso.a,x_interp);
    epsz2_e = interp1(viszkcso.x,viszkcso.epsz2,x_interp);
    dd = interp1(viszkcso.x,viszkcso.d,x_interp);

    %gamma = 1 + (alpha*p)/viszkcso.E1 + (alpha*p)/viszkcso.E2*(1 - exp(-viszkcso.E2/viszkcso.eta2*dt));
    %beta = sqrt(viszkcso.E1/(alpha*ro*gamma));
    %a = sqrt(viszkcso.E1*gamma/(alpha*ro));
    %aa = sqrt(viszkcso.gamma*viszkcso.E1/(ro*viszkcso.alpha));

    val1 = ro*a*v + p;
    val2 = - ro*a*dt*(32*viszkcso.nu/dd^2*v + bet*alp/viszkcso.eta2*...
        (p-1e5)*exp(-viszkcso.E2/viszkcso.eta2*dt) - 2*bet*viszkcso.E2/viszkcso.eta2*epsz2_e*exp(-viszkcso.E2/viszkcso.eta2*dt));
    val = val1 + val2;

    % tesztval1 = 32*viszkcso.lambda/viszkcso.D^2*v;
    % tesztval2 = beta*alpha/viszkcso.eta2*p*(exp(viszkcso.E2/viszkcso.eta2*dt));
    % tesztval3 = 2*beta*viszkcso.E2/viszkcso.eta2*epsz2_e*exp(-viszkcso.E2/viszkcso.eta2*dt);
    %
    % fprintf('A viszkcso.dt változó értéke: %5d\n', viszkcso.dt);
    % fprintf('A p változó értéke: %5d\n', p);
    % fprintf('A v változó értéke: %5d\n', v);
    % fprintf('A val1 változó értéke: %5d\n', val1);
    % fprintf('A val2 változó értéke: %5d\n', val2);
    % fprintf('A val változó értéke: %5d\n', val);
    % fprintf('A ro*viszkcso.a*dt változó értéke: %5d\n', ro*a*dt);
    % fprintf('A exp(-viszkcso.E2/viszkcso.eta2*dt) változó értéke: %5d\n', exp(-viszkcso.E2/viszkcso.eta2*dt));
    % fprintf('A tesztval1 változó értéke: %5d\n', tesztval1);
    % fprintf('A tesztval2 változó értéke: %5d\n', tesztval2);
    % fprintf('A tesztval3 változó értéke: %5d\n', tesztval3);
    % fprintf('A ro változó értéke: %5d\n', ro);

    if isnan([val1,val2])
        fprintf('\n\n-------------------------------------');
        fprintf('\n  val1=%g  val1=%g');
        fprintf('\n-------------------------------------\n\n');
        error('Pánik indul....')
    end

    %dd = viszkcso.dd(N+1);
    A = dd^2*pi/4;

end

% Az elsõ 4 szám a fontos, a többi csak infó

%disp(a);

%out=[1e5,  viszkcso.a/viszkcso.A, val,  1, val, p, v, val1, val2];
%out=[1e5, a/viszkcso.A, val,  1, val, p, v, val1, val2];
out=[1e5, a/A, val,  1, val, p, v, val1, val2];

%fprintf('\n %s csõ állapota:\n',cso.tranziens_agelem_2csp.nev);
%disp([viszkcso.x; viszkcso.p/1e5; viszkcso.v]);
%fprintf('\n Ide kell C+ mentén interpolálni: dt=%g  cso.dt=%g  cso.dx=%g  cso.a=%g  x=%g,  p=%g,  v=%g\n',dt,cso.dt,cso.dx,cso.a,x_interp,p/1e5,v);
%pause