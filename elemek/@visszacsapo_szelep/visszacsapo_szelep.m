function vcssz = visszacsapo_szelep(varargin)

m0=varargin{5};
ro=varargin{4};
vcssz.dummy=1;

% Elem l�trehoz�sa
trag2 = tranziens_agelem_2csp(varargin{1},varargin{2},varargin{3},m0/ro);
vcssz = class(vcssz,'visszacsapo_szelep',trag2);

% Ha csak a s�r�s�get akarjuk �t�ll�tani:
vcssz.tranziens_agelem_2csp.ro = ro;
% Ha minden folyad�k param�tert (n�v,ro,nu,mu,B):
%csov=setfluid(csov,folynev);
