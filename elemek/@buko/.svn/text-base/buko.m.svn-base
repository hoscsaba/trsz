function buko = buko(varargin)
  
buko.h0     = varargin{4}; % bukoszint
buko.Cd     = varargin{5}; % atfolyasi tenyezo
buko.B      = varargin{6}; % jellemzo meret
buko.kitevo = varargin{7}; % vizszint kitevo
ro          = varargin{8};
m0          = varargin{9};

buko.is_working=0;

trag2 = tranziens_agelem_2csp(varargin{1},varargin{2},varargin{3},m0/ro);  
buko = class(buko,'buko',trag2);

buko.tranziens_agelem_2csp.ro = ro;
