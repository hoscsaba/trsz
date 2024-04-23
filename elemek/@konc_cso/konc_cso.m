function csov = konc_cso(varargin)
  
csov.D      = varargin{4}; % Átmérõ [m]
csov.L      = varargin{5}; % Hossz [m]        
csov.lambda = varargin{6}; % Lambda
ro          = varargin{7};
m0          = varargin{8};
csov.A   = csov.D^2*pi/4;

trag2 = tranziens_agelem_2csp(varargin{1},varargin{2},varargin{3},m0/ro);  
csov = class(csov,'konc_cso',trag2);

% Ha csak a sûrûséget akarjuk átállítani:
cso.tranziens_agelem_2csp.ro = ro;
% Ha minden folyadék paramétert (név,ro,nu,mu,B):
%csov=setfluid(csov,folynev);