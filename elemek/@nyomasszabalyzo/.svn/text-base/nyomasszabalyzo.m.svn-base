function nysz = nyomasszabalyzo(varargin)

m0=varargin{5};
ro=varargin{4};
nysz.jge   = varargin{6};
nysz.jgk   = varargin{7};
nysz.szab  = varargin{8};

switch nysz.szab
    case 'p'
        nysz.szcsp1 = varargin{9};
        nysz.szcsp2 =-1;
        nysz.ajel  = varargin{10}/1e5;
        nysz.pmax  = varargin{11}/1e5;
        nysz.szabP = varargin{12};
        nysz.szabI = varargin{13};
        nysz.szabD = varargin{14};
        nysz.szcspnev1=varargin{15};
        nysz.szcspnev2='';
        nysz.A=1;
        nysz.e=varargin{16};
        nysz.vmax=varargin{17};
        nysz.tbe=varargin{18};
        
    case 'dp'        
        nysz.szcsp1 = varargin{9};
        nysz.szcsp2 =varargin{10};
        nysz.ajel  = varargin{11}/1e5;
        nysz.pmax  = varargin{12}/1e5;
        nysz.szabP = varargin{13};
        nysz.szabI = varargin{14};
        nysz.szabD = varargin{15};
        nysz.szcspnev1=varargin{16};
        nysz.szcspnev2=varargin{17};
        nysz.A=1;
        nysz.e=varargin{18};
        nysz.vmax=varargin{19};
        nysz.tbe=varargin{20};        
end

nysz.tt=[-4 -3 -2 -1]*1e-3;
nysz.xx=[0 0 0 0];
nysz.euj=nysz.e;
nysz.er =nysz.e;
nysz.err=nysz.e;
nysz.ii=0;

% Elem létrehozása
trag2 = tranziens_agelem_2csp(varargin{1},varargin{2},varargin{3},m0/ro);
nysz = class(nysz,'nyomasszabalyzo',trag2);

% Ha csak a sûrûséget akarjuk átállítani:
nysz.tranziens_agelem_2csp.ro = ro;
% Ha minden folyadék paramétert (név,ro,nu,mu,B):
%csov=setfluid(csov,folynev);
