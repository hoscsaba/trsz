function [out,szip]=setport(flag,szip,varargin)
  
ro  = szip.tranziens_agelem_1csp.ro;
csp = szip.tranziens_agelem_1csp.csp;
p   = szip.tranziens_agelem_1csp.p;

switch flag
  case 3
   xr = varargin{1};
   t  = varargin{2}(1);
   dt = varargin{2}(2);
   mar = varargin{3};
   p=get(mar,'p',szip.tranziens_agelem_1csp.csp);
   szip.V=szip.V-dt*xr;
%   fprintf('\n t=%6.3e  %4s  p=%6.3f  szip.Q=%+6.4e  szip.V=%+6.4e',t,szip.nyitva,p/1e5,szip.tranziens_agelem_1csp.Q,szip.V);
   out{1}=1;
   
 case 1
  xr = varargin{1};
  t  = varargin{2}{1}{1}(1);
  dt = varargin{2}{1}{1}(2);
  mar = varargin{3};
  
  p=get(mar,'p',szip.tranziens_agelem_1csp.csp);

  if p<szip.pnyit
    szip.nyitva='igen';
    out{1}={0,0,szip.K/ro,szip.pnyit,[csp,-1e5]};
  elseif (szip.V<=0)
    szip.nyitva='nem';
    out{1}={1,0,0,0,[csp,0]};
    szip.V=0;
  else
    out{1}={0,0,szip.K/ro,szip.pnyit,[csp,-1e5]};
  end
  
end
