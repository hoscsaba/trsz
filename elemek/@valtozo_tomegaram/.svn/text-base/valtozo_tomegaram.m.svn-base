function ta = valtozo_tomegaram(varargin)

ro = varargin{4};
m0 = varargin{5}; 
ta.tt   = varargin{6};
ta.mm   = varargin{7};
ta.jgpsz= length(ta.tt);
ta.m = interp1(ta.tt,ta.mm,0);

trag1 = tranziens_agelem_1csp(varargin{1},varargin{2},0);
ta = class(ta,'valtozo_tomegaram',trag1);
ta.tranziens_agelem_1csp.ro = ro;
