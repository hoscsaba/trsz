function out = karp(viszkcso,dt,dtmin)

%fprintf('A karpv függvény meghívása. "dt" értéke: %d \n',dt)
%pause

ro = viszkcso.tranziens_agelem_2csp.ro;

if dt < 0
    if abs(dt) < dtmin;
        dt = dtmin;
    end
end

if dt <= 0
    %
    %    dt=abs(dt);
    %else
    %disp(dt)
    %error('Nagy a baj...');
    if abs(dt)<dtmin;
        dt = abs(dt);
    else
        fprintf('\nNegatív idõlépés (karp): %8.3e hely: %s\n',dt,viszkcso.tranziens_agelem_2csp.nev);
        error('Hiba!')
    end
    
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

elseif dt > 0

    g = 9.81;
    
    x = viszkcso.x;
    h = viszkcso.h;
    L = viszkcso.L;
    aa = viszkcso.a;
    vv = viszkcso.v;
    N = viszkcso.N;
    nu = viszkcso.nu;
    p0 = viszkcso.p0;
    d0 = viszkcso.d0;
    s0 = viszkcso.s;
    E2 = viszkcso.E2;
    eta2 = viszkcso.eta2;

    %delta_t = (L - x(N))/(vv(N) + aa(N));
    delta_t = viszkcso.dtuj;
%     if tcso > takt
%         fprintf('A(z) %s cso eloresietett!',viszkcso.tranziens_agelem_2csp.nev)
%         error('Vegzetes hiba!')
%     end

    if dt > delta_t*1.01
        disp(dt)
        disp(delta_t)
        error('Extrapolacio!!!')
    end
        
    % A peremen interpolálunk. (A csõ végén kapcsolódó csõ lép...)

    x_interp = viszkcso.L-viszkcso.dx*dt/delta_t;
    
    p = interp1(viszkcso.x,viszkcso.p,x_interp);
    v = interp1(viszkcso.x,viszkcso.v,x_interp);

    %dzdx=cso.dzdx(cso.N);
    %ro = viszkcso.tranziens_agelem_2csp.ro;

    if isnan(v)
        fprintf('\n\n')
        v
        x_interp
        delta_t
        dt
        disp(viszkcso.tranziens_agelem_2csp.nev);
        disp(viszkcso.x);
        disp(viszkcso.v);
        disp(viszkcso.a*dt);
        error('PÁNIK!!!')
    end

    hL = interp1(viszkcso.x,viszkcso.h,x_interp);    
    a = interp1(viszkcso.x,viszkcso.a,x_interp);
    epsz2 = interp1(viszkcso.x,viszkcso.epsz2,x_interp);
    epsz = interp1(viszkcso.x,viszkcso.epsz,x_interp);
    d = interp1(viszkcso.x,viszkcso.d,x_interp);

    alp = d0/s0*(2*epsz+1);
    gam = 2*epsz+1;
    G = -2*E2/eta2*epsz2*exp(-E2/eta2*dt) + alp/eta2*(p - p0)*exp(-E2/eta2*dt);
  
    hP = h(N+1);
    xP = L; xL = x_interp;
    
    val1 = ro*a*v + p;
    val2 = - ro*a*dt*( g*(hP-hL)/(xP-xL) + 32*nu/d^2*v + a/gam*G );

    val = val1 + val2;

    if max(isnan([val1,val2])) == 1
        fprintf('\n\n-------------------------------------');
        dt
        gam
        (xP-xL)
        fprintf('\n  val1=%g  val2=%g',val1,val2);
        fprintf('\n-------------------------------------\n\n');
        error('Pánik indul....')
    end

    A = d^2*pi/4;
        
end

% Az elsõ 4 szam a fontos, a többi csak info

out=[1e5, a/A, val,  1, val, p, v, val1, val2];

