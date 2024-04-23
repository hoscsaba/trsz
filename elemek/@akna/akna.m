function ak = akna(varargin)

%ak.nev = varargin{1};
%ak.csp = varargin{2};
ak.A  = varargin{3}; % Alapterulet
ak.hmin=varargin{4}; % Fenekszint
ak.hmax=varargin{5}; % Fedlapszint
ak.y0 = varargin{6}+ak.hmin; % Kezdeti szint
ro    = varargin{7}; % Nyomas
m0    = varargin{8}; % Kezdeti tomegaram
ak.rajz  = varargin{9}; % Rajz

ak.y = ak.y0;
ak.yr = ak.y0;
ak.dte = 1e-5;
ak.Q = 0;
ak.Qr = 0;
ak.Qrr = 0;

ak.ures = 0;
ak.tele = 0;

trag1 = tranziens_agelem_1csp(varargin{1},varargin{2},m0/ro);
ak = class(ak,'akna',trag1);
ak.tranziens_agelem_1csp.ro = ro;
