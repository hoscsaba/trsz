function fojt = fojtas(varargin)

fojt.K0  = varargin{4};
fojt.K1  = varargin{5};
ro       = varargin{6};
m0       = varargin{7};

% Ez a sor így marad
trag2 = tranziens_agelem_2csp(varargin{1},varargin{2},varargin{3},m0/ro);

fojt = class(fojt,'fojtas',trag2);

% Ha csak a sûrûséget akarjuk átállítani:

fojt.tranziens_agelem_2csp.ro = ro;
% Ha minden folyadék paramétert (név,ro,nu,mu,B):
%csov=setfluid(csov,folynev);
