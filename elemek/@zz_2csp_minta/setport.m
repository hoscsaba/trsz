function [out,vcssz]=setport(flag,vcssz,varargin)
g   = 9.81;
ro  = vcssz.tranziens_agelem_2csp.ro;
csp = vcssz.tranziens_agelem_2csp.csp;
%xr  = sziv.tranziens_agelem_2csp.Qr;

if flag==-1 % inicializálás
  out{1}=1;
  
elseif flag==3

  % Ide jon a belso valtozok atallitasa, pl. szivattyú kifutás
  % esetén a núj
  xr=varargin{1};
  t  = varargin{2}(1);
  dt = varargin{2}(2);
  
%  sziv.n=sziv.nn;
%  sziv.res=[sziv.res; t sziv.n];
  out{1}=1;
%  fprintf('\nt=%5.3e   dt=%5.3e   n=%5.3e  Q=%+6.3e',t,dt,sziv.n,xr);
  
else
  
  xr = varargin{1};
  t  = varargin{2}{1}{1}(1);
  dt = varargin{2}{1}{1}(2);
  
  if xr>0, alfa=1e-5;
  else alfa=1e10; end

  out{1} = {alfa, 0, 0, 0, [csp(1),-1], [csp(2),1] };

end
