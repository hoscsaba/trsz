function lu = legust(varargin)

%lu.p  = varargin{3};
ro    = varargin{4}; % Nyomas
m0    = varargin{5}; % Tomegaram
lu.n  = varargin{6}; % Politropikus kitevo
lu.V0 = varargin{7}; % Kezdeti gáztérfogat
lu.p0 = varargin{8}; % Kezdeti gáznyomás
lu.A  = varargin{9}; % Felület
lu.l  = varargin{10}; % Alappont és befolyás magasságkülönbsége
lu.H  = varargin{11}; % Teljes belsõ magasság		     

lu.p0 = lu.p0 + ro*9.81*(lu.A*lu.H-lu.V0)/lu.A;

lu.extrap = 0;

lu.ize = lu.p0*(lu.V0)^(lu.n);
lu.V   = lu.V0;  lu.p   = lu.p0;
lu.Vr  = lu.V0;  lu.pr  = lu.p0; 
lu.Vrr = lu.V0;  lu.prr = lu.p0; 
lu.dte=1e-5; lu.dtee=1e-5;

trag1 = tranziens_agelem_1csp(varargin{1},varargin{2},m0/ro);
lu = class(lu,'legust',trag1);

lu.tranziens_agelem_1csp.ro = ro;
%lu.tranziens_agelem_1csp.p = varargin{9};