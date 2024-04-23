function vfojt = vez_fojtas(varargin)
  
vfojt.A      = varargin{5};
vfojt.jgpsz  = varargin{6};
vfojt.jgpszt = varargin{7};
ro           = varargin{4};
m0           = varargin{8};
vfojt.epsK   = varargin{9};
vfojt.K      = varargin{10};
vfojt.t_vf   = varargin{11};
vfojt.epst   = varargin{12};

trag2 = tranziens_agelem_2csp(varargin{1},varargin{2},varargin{3},m0/ro);
vfojt = class(vfojt,'vez_fojtas',trag2);

% Ha csak a sûrûséget akarjuk átállítani:
vfojt.tranziens_agelem_2csp.ro = ro;
% Ha minden folyadék paramétert (név,ro,nu,mu,B):
%csov=setfluid(csov,folynev);
