function out = karm(viszkcso,dt,dtmin)

%fprintf('A karmv függvény meghívása...\n')
%pause

% Az idõlépés helyességének ellenõrzése

ro = viszkcso.tranziens_agelem_2csp.ro;

if dt < 0
    if abs(dt) < dtmin;
        dt = dtmin;
    end
end

if dt <= 0
    %if abs(dt)<dtmin;
    %    dt=dtmin;
    %fprintf('\nWe experienced some problems...\n');
    %else
    %disp(dt);
    %error('Nagy a baj...');

    dt = abs(dt);
    delta_t = viszkcso.dt;

    N  = viszkcso.N;
    p_uj = viszkcso.p(1);
    v_uj = viszkcso.v(1);
    a_uj = viszkcso.a(1);
    d_uj = viszkcso.d(1);
    p_e = viszkcso.p1;
    v_e = viszkcso.v1;
    a_e = viszkcso.a1;
    d_e = viszkcso.d1;

    %fprintf('p_e: %8.4f \t p_uj: %8.4f\n', p_e, p_uj);
    %fprintf('v_e: %8.4f \t v_uj: %8.4f\n', v_e, v_uj);
    %fprintf('a_e: %8.4f \t a_uj: %8.4f\n', a_e, a_uj);

    % Lineáris interpoláció (j és j-1 között)

    L1 = dt/delta_t;
    L2 = (delta_t-dt)/delta_t;
    p = L1*p_e + L2*p_uj;
    v = L1*v_e + L2*v_uj;
    a = L1*a_e + L2*a_uj;
    dd = L1*d_e + L2*d_uj;

    val1 = 0;
    val2 = 0;
    val = ro*a*v - p;

    A = dd^2*pi/4;

    %end
elseif dt > 0

    % A peremen (i = 1) interpolálunk (a csõ elejéhez kapcsolódó csõ lép...)

    x = viszkcso.x;
    aa = viszkcso.a;
    vv = viszkcso.v;

    delta_t = x(2)/(aa(2) - vv(2));

    %x_interp = viszkcso.dx*dt/viszkcso.dt;
    %x_interp = viszkcso.dx*dt/viszkcso.ddt1;

    x_interp = viszkcso.dx*dt(1)/delta_t;
    p = interp1(viszkcso.x,viszkcso.p,x_interp);
    v = interp1(viszkcso.x,viszkcso.v,x_interp);
    %dzdx=cso.dzdx(1);
    ro = viszkcso.tranziens_agelem_2csp.ro;

%     aa
%     vv
%     delta_t
%     x_interp
%     v
%     pause

    if isnan(v)
        fprintf('\n\n --> %s\n',viszkcso.tranziens_agelem_2csp.nev);
        disp(viszkcso.x);
        disp(viszkcso.v);
        disp(viszkcso.a*dt);
        error('PÁNIK!!!')
    end

    % alpha = viszkcso.alpha1;
    % beta = viszkcso.beta1;
    % gamma = viszkcso.gamma1;
    % a = viszkcso.a1;
    % epsz2_e = viszkcso.epsz21;

    alp = interp1(viszkcso.x,viszkcso.alp,x_interp);
    bet = interp1(viszkcso.x,viszkcso.bet,x_interp);
    gam = interp1(viszkcso.x,viszkcso.gam,x_interp);
    a = interp1(viszkcso.x,viszkcso.a,x_interp);
    epsz2_e = interp1(viszkcso.x,viszkcso.epsz2,x_interp);
    dd = interp1(viszkcso.x,viszkcso.d,x_interp);

    %val1 = p-ro*cso.a*v;
    %val2 = dt*( ro*cso.a*9.81*(cso.h(2)-cso.h(1))/cso.dx + ro*cso.a*cso.lambda/2/cso.D*v*abs(v));

    %val1 = viszkcso.E1/(viszkcso.beta*viszkcso.alpha)*v - p;

    val1 = ro*a*v - p;
    val2 = - ro*a*dt*(32*viszkcso.nu/dd^2*v - bet*alp/viszkcso.eta2*...
        (p-1e5)*exp(-viszkcso.E2/viszkcso.eta2*dt) + 2*bet*viszkcso.E2/viszkcso.eta2*epsz2_e*exp(-viszkcso.E2/viszkcso.eta2*dt));
    %val = val1 + val2;
    val = -(val1 + val2);

    % fprintf('A viszkcso.dt változó értéke: %5d\n', viszkcso.dt);
    % fprintf('A p változó értéke: %5d\n', p);
    % fprintf('A v változó értéke: %5d\n', v);
    % fprintf('A val1 változó értéke: %5d\n', val1);
    % fprintf('A val2 változó értéke: %5d\n', val2);
    % fprintf('A val változó értéke: %5d\n', val);
    % fprintf('A ro*viszkcso.a*dt változó értéke: %5d\n', ro*a*dt);
    % fprintf('A exp(-viszkcso.E2/viszkcso.eta2*dt) változó értéke: %5d\n',
    % exp(-viszkcso.E2/viszkcso.eta2*dt));

    if isnan([val1,val2])
        fprintf('\n\n---------------------------------------------------------------------------------');
        fprintf('\n  val1=%g  | val2=%g -> dt=%g,  dzdx-es tag=%g,  lambda-s tag= %g ',val1,val2,dt,ro*cso.a*9.81*dzdx,ro*cso.a*cso.lambda/2/cso.D*v*abs(v));
        if isnan(ro*cso.a*cso.lambda/2/cso.D*v*abs(v))
            fprintf('\n\t  ro=%g  a=%g  lambda=%g  D=%g  v=%g ',ro,cso.a,cso.lambda,cso.D,v);
        end
        fprintf('\n-----------------------------------------------------------------------------------\n\n');
        error('Pánik indul....')
    end

    %dd = viszkcso.dd(1);
    A = dd^2*pi/4;
end
% Az elsõ 4 szám a fontos, a többi csak infó
% disp(a);
% out=[1e5, -viszkcso.a/viszkcso.A, val, -1, val, p, v, val1, val2];
% out=[1e5, -a/viszkcso.A, val, -1, val, p, v, val1, val2];
out=[1e5, -a/A, val, -1, val, p, v, val1, val2];

%fprintf('\n %s csõ állapota:\n',cso.tranziens_agelem_2csp.nev)
%disp([viszkcso.x; viszkcso.p/1e5; viszkcso.v]);
%fprintf('\n Ide kell C- mentén interpolálni: dt=%g  cso.dt=%g  cso.dx=%g  cso.a=%g  x=%g,  p=%g,  v=%g\n',dt,cso.dt,cso.dx,cso.a,x_interp,p/1e5,v);
%pause