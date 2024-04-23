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
    E1 = viszkcso.E1;
    E2 = viszkcso.E2;
    eta2 = viszkcso.eta2;

    delta_t = viszkcso.dtuj;

    if dt > delta_t*1.01
        disp(dt)
        disp(delta_t)
        error('Extrapolacio!!!')
    end
    
    % Az ot parameter beolvasasa filebol
   
    par = textread('params.dat','%f',5);
    
    E1 = par(1);
    E2 = par(2);
    E3 = par(3);
    eta1 = par(4);
    eta2 = par(5);
   
    K1 = eta1*eta2/E2; K2 = (eta2*E1/E2+eta1); K3 = E1;
    m1 = (K2 + sqrt(K2^2-4*K1*K3)) / (2*K1);
    m2 = sqrt(K2^2-4*K1*K3) / K1;    
    
    %delta_t = (L - x(N))/(vv(N) + aa(N));

    % A peremen interpolálunk. (A csõ végén kapcsolódó csõ lép...)

    x_interp = viszkcso.L-viszkcso.dx*dt/delta_t;
    p = interp1(viszkcso.x,viszkcso.p,x_interp);
    v = interp1(viszkcso.x,viszkcso.v,x_interp);
    %dzdx=cso.dzdx(cso.N);
    %ro = viszkcso.tranziens_agelem_2csp.ro;

    if isnan(v)
        fprintf('\n\n')
        disp(viszkcso.x);
        disp(viszkcso.v);
        disp(viszkcso.a*dt);
        error('PÁNIK!!!')
    end

    hL = interp1(viszkcso.x,viszkcso.h,x_interp);    
    a = interp1(viszkcso.x,viszkcso.a,x_interp);
    epsz1 = interp1(viszkcso.x,viszkcso.epsz2,x_interp);
    epsz1v = interp1(viszkcso.x,viszkcso.epsz1v,x_interp);
    epsz = interp1(viszkcso.x,viszkcso.epsz,x_interp);
    d = interp1(viszkcso.x,viszkcso.d,x_interp);

    %gamma = 1 + (alpha*p)/viszkcso.E1 + (alpha*p)/viszkcso.E2*(1 - exp(-viszkcso.E2/viszkcso.eta2*dt));
    %beta = sqrt(viszkcso.E1/(alpha*ro*gamma));
    %a = sqrt(viszkcso.E1*gamma/(alpha*ro));
    %aa = sqrt(viszkcso.gamma*viszkcso.E1/(ro*viszkcso.alpha));

    alp = d0/s0*(2*epsz+1);
    gam = 2*epsz+1;
    %G = -2*E2/eta2*epsz2*exp(-E2/eta2*dt) + alp/eta2*(p - p0)*exp(-E2/eta2*dt);
    
    epsz1s = epsz1; %epsz1v = 0;
    szig3 = 0.5*alp*(p-p0);
    G = (2*E1*K1^2*m2)^-1*exp(-m1*dt)*( -eta1*(epsz1v*(K2^2 - 2*K1*K3 + K2*K1*m2 + ...
        exp(m2*dt)*K2*K1 - exp(m2*dt)*(K2^2 - 2*K1*K3)) + ...
        (K2 - exp(m2*dt)*K2 + (1 + exp(m2*dt))*m2*K1)*(epsz1s*K3 - szig3)) + ...
        2*E1*K1*( epsz1v*K2 + 2*epsz1s*K3 + epsz1v*m2*K1 + exp(m2*dt)*epsz1v*m2*K1 - ...
        exp(m2*dt)*(epsz1v*K2 + 2*epsz1s*K3 - 2*szig3) - 2*szig3));
    
    hP = h(N+1);
    xP = L; xL = x_interp;
    
    val1 = ro*a*v + p;
    val2 = - ro*a*dt*( g*(hP-hL)/(xP-xL) + 32*nu/d^2*v + a/gam*G );
    
%     val2 = - ro*a*dt*(32*viszkcso.nu/dd^2*v + bet*alp/viszkcso.eta2*...
%         (p-1e5)*exp(-viszkcso.E2/viszkcso.eta2*dt) - 2*bet*viszkcso.E2/viszkcso.eta2*epsz2_e*exp(-viszkcso.E2/viszkcso.eta2*dt));
    val = val1 + val2;

    if isnan([val1,val2])
        fprintf('\n\n-------------------------------------');
        fprintf('\n  val1=%g  val1=%g');
        fprintf('\n-------------------------------------\n\n');
        error('Pánik indul....')
    end

    A = d^2*pi/4;

end

% Az elsõ 4 szam a fontos, a többi csak info

out=[1e5, a/A, val,  1, val, p, v, val1, val2];
