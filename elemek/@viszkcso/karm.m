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
    if abs(dt)<dtmin;
        dt = abs(dt);
    else
        fprintf('\nNegatív idõlépés (karm): %8.3e hely: %s\n',dt,viszkcso.tranziens_agelem_2csp.nev);
        error('Hiba!')
    end
    %else
    %disp(dt)
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
    d = L1*d_e + L2*d_uj;

    val1 = 0;
    val2 = 0;
    val = -ro*a*v + p;

    A = d^2*pi/4;

    %end
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

    %delta_t = x(2)/(aa(2) - vv(2));
    delta_t = viszkcso.dtuj;
    
    if dt > delta_t*1.01
        disp(dt)
        disp(delta_t)
        error('Extrapolacio!!!')
    end
    
    % A peremen interpolálunk. (A csöelejehez kapcsolodo csö lep...)

    x_interp = viszkcso.dx*dt/delta_t;
    
    p = interp1(viszkcso.x,viszkcso.p,x_interp);
    v = interp1(viszkcso.x,viszkcso.v,x_interp);
    
    %dzdx=cso.dzdx(cso.N);
    %ro = viszkcso.tranziens_agelem_2csp.ro;

    if isnan(v)
        fprintf('\n\n')
        disp(viszkcso.x);
        disp(viszkcso.v);
        disp(cso.a*dt);
        error('PÁNIK!!!')
    end

    hR = interp1(viszkcso.x,viszkcso.h,x_interp);
    a = interp1(viszkcso.x,viszkcso.a,x_interp);
    epsz2 = interp1(viszkcso.x,viszkcso.epsz2,x_interp);
    epsz = interp1(viszkcso.x,viszkcso.epsz,x_interp);
    d = interp1(viszkcso.x,viszkcso.d,x_interp);

    alp = d0/s0*(2*epsz+1);
    gam = 2*epsz+1;
    
    G = -2*E2/eta2*epsz2*exp(-E2/eta2*dt) + alp/eta2*(p - p0)*exp(-E2/eta2*dt);
       
    hP = h(1); %hR = 0;
    xP = 0; xR = x_interp;

    val1 = p - ro*a*v;
    val2 = ro*a*dt*( g*(hP-hR)/(xP-xR) + 32*nu/d^2*v - a/gam*G );  

    val = val1 + val2;
   
    if isnan([val1,val2])
        fprintf('\n\n---------------------------------------------------------------------------------');
        fprintf('\n  val1=%g  | val2=%g -> dt=%g,  dzdx-es tag=%g,  lambda-s tag= %g ',val1,val2,dt,ro*cso.a*9.81*dzdx,ro*cso.a*cso.lambda/2/cso.D*v*abs(v));
        if isnan(ro*cso.a*cso.lambda/2/cso.D*v*abs(v))
            fprintf('\n\t  ro=%g  a=%g  lambda=%g  D=%g  v=%g ',ro,cso.a,cso.lambda,cso.D,v);
        end
        fprintf('\n-----------------------------------------------------------------------------------\n\n');
        error('Pánik indul....')
    end

    A = d^2*pi/4;
    
end

out=[1e5, -a/A, val, -1, val, p, v, val1, val2];
