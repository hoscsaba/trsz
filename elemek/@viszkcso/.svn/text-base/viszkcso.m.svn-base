function viszkcso = viszkcso(varargin)

viszkcso.d0      = varargin{4};   % Átmérõ [m]
viszkcso.nu     = varargin{5};   % Lambda
viszkcso.s      = varargin{6};   % Falvastagság [m]
viszkcso.L      = varargin{7};   % Csõhossz [m]
viszkcso.N      = varargin{8}-1; % Csõdarabok száma (op-1)
ro              = varargin{9};
viszkcso.E1     = varargin{13};
viszkcso.E2     = varargin{14};
viszkcso.eta2   = varargin{15};
viszkcso.p0 = 1e5;
p0 = viszkcso.p0;

viszkcso.dx = viszkcso.L/viszkcso.N;

A  = viszkcso.d0^2*pi/4;
Q = varargin{10}*A;

trag2 = tranziens_agelem_2csp(varargin{1},varargin{2},varargin{3},Q);

viszkcso.p = varargin{11};



for i = 1:viszkcso.N+1
    viszkcso.x(i)=(i-1)*viszkcso.dx;
    viszkcso.v(i)=varargin{10};
    viszkcso.alp(i) = viszkcso.d0/viszkcso.s;
    viszkcso.gam(i) = (1 + viszkcso.alp(i)*(viszkcso.p(i) - p0)/viszkcso.E1);
    viszkcso.bet(i) = sqrt(viszkcso.E1/(viszkcso.alp(i)*ro*viszkcso.gam(i)));
    viszkcso.a(i)  = viszkcso.E1/(viszkcso.bet(i)*viszkcso.alp(i)*ro);
    viszkcso.ddt(i) = viszkcso.dx/viszkcso.a(i);
    viszkcso.epsz(i) = 0;
    viszkcso.epsz2(i) = 0;
    viszkcso.epsz1v(i) = 0;
    viszkcso.dseb(i) = 0;
    viszkcso.d(i) = viszkcso.d0;
    dt(i) = viszkcso.dx/viszkcso.a(i);    
    viszkcso.A(i) = A;
end

agnev = varargin{1};
fnev = [agnev '.rst'];

if exist(fnev)
    data = load(fnev);
    p1_old = data(end,2);
    v1_old = data(end,4);
    d1_old = data(end,10);
    p2_old = data(end,3);
    v2_old = data(end,5);
    d2_old = data(end,11);
    viszkcso.p(1:end) = linspace(p1_old,p2_old,length(viszkcso.p));
    viszkcso.v(1:end) = linspace(v1_old,v2_old,length(viszkcso.v));
    viszkcso.d(1:end) = linspace(d1_old,d2_old,length(viszkcso.d));
end

viszkcso.dtki = 0;
viszkcso.dt = min(dt);

% Uj idolepes meghatarozasa

x = viszkcso.x; v = viszkcso.v; a = viszkcso.a; L = viszkcso.L;
N = viszkcso.N;

dtuj(1) = x(2)/(-v(2) + a(2));
dtuj(N+1) = (L - x(N))/(v(N) + a(N));
for i = 2:N
     xL = x(i-1); xR = x(i+1);
     
     vL = v(i-1); vR = v(i+1);    
     aL = a(i-1); aR = a(i+1);

     tPuj = (xR - xL - 0*(vR - aR)...
         + 0*(vL + aL))/((vL + aL) - (vR - aR));
     dtuj(i) = tPuj - 0;
end

viszkcso.dtuj = min(dtuj);

%aa = max(viszkcso.a); % !!! min/max?
%viszkcso.dt = viszkcso.dx/aa;

viszkcso.h = varargin{12};

%viszkcso.tipus = varargin{21};

for i=1:viszkcso.N
    viszkcso.dzdx(i)=(viszkcso.h(i+1)-viszkcso.h(i))/viszkcso.dx;
end
viszkcso.dzdx(viszkcso.N+1)=0;

viszkcso.t = 0;
viszkcso.v1 = viszkcso.v(1);
viszkcso.p1 = viszkcso.p(1);
viszkcso.alp1  = viszkcso.alp(1);
viszkcso.bet1   = viszkcso.bet(1);
viszkcso.gam1  = viszkcso.gam(1);
viszkcso.a1  = viszkcso.a(1);
viszkcso.d1 = viszkcso.d(1);
%viszkcso.dt1 = viszkcso.dt(1);
viszkcso.epsz_1 = viszkcso.epsz(1);
viszkcso.epsz2_1 = viszkcso.epsz2(1);

viszkcso.vN = viszkcso.v(viszkcso.N+1);
viszkcso.pN = viszkcso.p(viszkcso.N+1);
viszkcso.alpN  = viszkcso.alp(viszkcso.N+1);
viszkcso.betN   = viszkcso.bet(viszkcso.N+1);
viszkcso.gamN  = viszkcso.gam(viszkcso.N+1);
viszkcso.aN  = viszkcso.a(viszkcso.N+1);
viszkcso.dN = viszkcso.d(viszkcso.N+1);
%viszkcso.dtN = viszkcso.dt(viszkcso.N+1);
viszkcso.epsz_N = viszkcso.epsz(viszkcso.N+1);
viszkcso.epsz2_N = viszkcso.epsz2(viszkcso.N+1);

viszkcso.all = {viszkcso.x.'; viszkcso.v.'; viszkcso.a.'; ...
                viszkcso.t; viszkcso.N; viszkcso.L};

viszkcso.warnings = 1;

viszkcso.pf_eleje_tipus='karakterisztika';
viszkcso.pf_vege_tipus='karakterisztika';

viszkcso = class(viszkcso,'viszkcso',trag2);

viszkcso.tranziens_agelem_2csp.ro = ro;

viszkcso.tranziens_agelem_2csp.p =[viszkcso.p1 viszkcso.pN];