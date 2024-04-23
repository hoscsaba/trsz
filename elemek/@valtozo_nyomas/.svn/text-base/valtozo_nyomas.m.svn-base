function nyom = valtozo_nyomas(varargin)

ro   = varargin{4};
nyom.tt   = varargin{6};
nyom.pp   = varargin{7};
nyom.jgpsz= length(nyom.tt);
nyom.p = interp1(nyom.tt,nyom.pp,0);

trag1 = tranziens_agelem_1csp(varargin{1},varargin{2},0);
nyom = class(nyom,'valtozo_nyomas',trag1);
nyom.tranziens_agelem_1csp.ro = ro;
