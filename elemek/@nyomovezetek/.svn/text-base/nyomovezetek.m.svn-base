function csov = nyomovezetek(varargin)

%csov.nev = varargin{1};
csov.cspe   = varargin{2}; % Csomopont eleje
csov.cspv   = varargin{3}; % Csomopont vege
csov.D      = varargin{4}; % Atmero [m]
csov.L      = varargin{5}; % Hossz [m]        
csov.lambda = varargin{6}; % Lambda
csov.ze     = varargin{7};
csov.zv     = varargin{8};
ro          = varargin{9};
m0          = varargin{10};
csov.A   = csov.D^2*pi/4;

csov.init = 0;
csov.aknae = 0;
csov.aknav = 0;
csov.sziv{1} = 0;

trag2 = tranziens_agelem_2csp(varargin{1},varargin{2},varargin{3},m0/ro);  
csov = class(csov,'nyomovezetek',trag2);
