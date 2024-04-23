function sziv = szivattyu(varargin)

%sziv.cspe = varargin{2};
%sziv.cspv = varargin{3};
sziv.jgQ = varargin{4};
sziv.jgH = varargin{5};
sziv.Ds  = varargin{6};	
sziv.Dn  = varargin{7};
ro       = varargin{8};
m0       = varargin{9};
sziv.tranziens = varargin{10};

% Itt inicializ.A��ljuk az ��sszes k��s.B��bbi adatot.
sziv.n = 1;
sziv.jgP = sziv.jgQ.*sziv.jgH*9.81*ro;
sziv.teta = 1;
sziv.tki  = 0;
sziv.res  = [];
sziv.Mb   = 1; 
sziv.Mi   = 1; 
sziv.nb   = 1;
sziv.nsz  = 1;
sziv.tbe  = 1;
sziv.tind = 1;
sziv.M    = 0;
sziv.Mm   = 0;

switch sziv.tranziens
    
    case 0  % stacioner
        sziv.n    = 1;
    
    case 1  % kifutas
        sziv.jgP  = varargin{11}*1000;
        sziv.n    = varargin{12}/60;
        sziv.teta = varargin{13};
        sziv.tki  = varargin{14};
    
    case 2  % inditas, nincs Mm(n) jg.
        sziv.jgP  = varargin{11}*1000;
        sziv.n    = varargin{12}/60;
        sziv.teta = varargin{13};
        sziv.Mb   = varargin{14};
        sziv.nb   = varargin{15}/60;
        sziv.nsz  = varargin{16}/60;
    
    case 4  % frekivaltos inditas, kozelito motor jg.
        sziv.jgP  = varargin{11}*1000;
        sziv.n    = varargin{12}/60;
        sziv.teta = varargin{13};
        sziv.Mb   = varargin{14};
        sziv.nb   = varargin{15}/60;
        sziv.Mi   = varargin{16};
        sziv.nsz  = varargin{17}/60;
        sziv.tbe  = varargin{18};
        sziv.tind = varargin{19};
    
    case 5 % szintkapcsolos szivattyu
        sziv.jgP   = varargin{11}*1000;
        sziv.n     = varargin{12}/60;
        sziv.hbe   = varargin{13};
        sziv.hki   = varargin{14};
        sziv.uzem  = varargin{15};
        sziv.init  = 0;
        sziv.aknae = 0; 
        
    case 6
        sziv.jgP  = varargin{11}*1000; % teljesitmeny
        sziv.nmax = varargin{12}/60;   % maximalis fordulatszam
        sziv.pcsp = varargin{13};      % az eloirt nyomasu pont NEVE
        sziv.pny  = varargin{14};      % az eloirt nyomas
        sziv.ParP = varargin{15};      % aranyos tag
        sziv.ParI = varargin{16};      % integralo tag
        sziv.nn   = varargin{17}/60;   % kezdeti fordulatszam
        sziv.n    = sziv.nmax;
        sziv.intdelay=2;
        sziv.intszumma = 0;
        sziv.kezd = 1;
        sziv.akt = 1;
        sziv.adat(1)=0;
        sziv.ido(1) = 0;
        sziv.told = 0;
        sziv.Huj = (sziv.pny - 1e5)/ro/9.81;
        sziv.Hregi = sziv.Huj;
        sziv.pp = sziv.pny;
        sziv.deltap = 0;
end

sziv.n0 = sziv.n;


% Belso valtozok iteraciokhoz
sziv.nn=sziv.n;
sziv.MM=sziv.M;
sziv.MMm=sziv.Mm;
% Affin jelleggorbe szamitasa
sziv.qq=[sziv.jgQ/sziv.n];
sziv.hh=[sziv.jgH/sziv.n^2];
sziv.jgP=[sziv.jgP]; % csak a biztonsag kerdveert...
sziv.mm=sziv.jgP/(2*pi*sziv.n)/sziv.n^2;

% Inditas eseten nullazni kell a forulatszamot
if sziv.tranziens==4, sziv.n=1/60; sziv.nn=1/60; end

for i=1:length(sziv.qq)-1
    sziv.a(i)=(sziv.hh(i+1)-sziv.hh(i))/(sziv.qq(i+1)-sziv.qq(i));
    sziv.b(i)=sziv.hh(i)-sziv.a(i)*sziv.qq(i);
%    sziv.c(i)=(sziv.mm(i+1)-sziv.mm(i))/(sziv.qq(i+1)-sziv.qq(i));
%    sziv.d(i)=sziv.mm(i)-sziv.c(i)*sziv.qq(i);	      
end
mer=1;
sziv.a=[sziv.a -mer];
sziv.b=[sziv.b mer*sziv.qq(length(sziv.qq))];

if sziv.tranziens==2 || sziv.tranziens==3, sziv.n=0.1; end

% Elem l��trehoz��sa
trag2 = tranziens_agelem_2csp(varargin{1},varargin{2},varargin{3},m0/ro);
sziv = class(sziv,'szivattyu',trag2);

%sziv.tranziens_agelem_2csp.p(2) = sziv.tranziens_agelem_2csp.p(1)+interp1(sziv.qq,sziv.hh,m0/ro)*1000*9.81;

% Ha csak a s��r��s��get akarjuk ��t��ll��tani:
sziv.tranziens_agelem_2csp.ro = ro;

% Ha minden folyad��k param��tert (n��v,ro,nu,mu,B):
%csov=setfluid(csov,folynev);

