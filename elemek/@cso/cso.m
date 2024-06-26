function cso = cso(varargin)

% GAZKIVALAST LEHESSEN OPCIOKENT MEGADNI A BEMENETI FAJLBAN

cso.D      = varargin{4}; % Atmero [m]
cso.lambda = varargin{5}; % Lambda
cso.s      = varargin{6}; % Falvastagsag [m]   
cso.L      = varargin{7}; % Csohossz [m]
cso.Ecso   = varargin{8}; % Cso rug. mod. [Pa]    
cso.N      = varargin{9}-1; % Cs�darabok sz�ma (op-1)
ro         = varargin{10};
cso.Ef     = varargin{11};

cso.dx = cso.L/cso.N;

cso.A  = cso.D^2*pi/4;
Q = varargin{12}*cso.A;

% Adatfajlba kiiras
cso.dtki = 0.0;

trag2 = tranziens_agelem_2csp(varargin{1},varargin{2},varargin{3},Q);  

Er=1/(1/cso.Ef+cso.D/cso.s/cso.Ecso);
cso.a  = sqrt(Er/ro);
cso.dt = cso.dx/cso.a;

cso.p=varargin{13};
for i=1:cso.N+1
  cso.x(i)=(i-1)*cso.dx; 
  cso.v(i)=varargin{12};
end

cso.h=varargin{14};

for i=1:cso.N
    cso.dzdx(i)=(cso.h(i+1)-cso.h(i))/cso.dx;
end
cso.dzdx(cso.N+1)=0;

cso.t =0;
cso.v1=cso.v(1); 
cso.p1=cso.p(1);
cso.vN=cso.v(cso.N+1); 
cso.pN=cso.p(cso.N+1);

cso.warnings=0;

cso.pf_eleje_tipus='karakterisztika';
cso.pf_vege_tipus='karakterisztika';

cso = class(cso,'cso',trag2);

cso.tranziens_agelem_2csp.ro = ro;
cso.tranziens_agelem_2csp.p =[cso.p1 cso.pN];
