function vcssz = visszacsapo_szelep(varargin)

m0=varargin{5};
ro=varargin{4};
vcssz.dummy=1;

% Elem létrehozása
trag2 = tranziens_agelem_2csp(varargin{1},varargin{2},varargin{3},m0/ro);
vcssz = class(vcssz,'visszacsapo_szelep',trag2);

% Ha csak a sûrûséget akarjuk átállítani:
vcssz.tranziens_agelem_2csp.ro = ro;
% Ha minden folyadék paramétert (név,ro,nu,mu,B):
%csov=setfluid(csov,folynev);
