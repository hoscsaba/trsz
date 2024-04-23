function nyom = nyomas(varargin)

nyom.p=varargin{3};

trag1 = tranziens_agelem_1csp(varargin{1},varargin{2},0);
nyom = class(nyom,'nyomas',trag1);
