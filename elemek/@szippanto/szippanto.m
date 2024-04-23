function szip = szippanto(varargin)

szip.pnyit  = varargin{3}; % Politropikus kitevo
szip.K=1;
ro=varargin{4};
szip.V=0;
szip.nyitva='nem';

trag1 = tranziens_agelem_1csp(varargin{1},varargin{2},0);
szip = class(szip,'szippanto',trag1);

szip.tranziens_agelem_1csp.ro = ro;
