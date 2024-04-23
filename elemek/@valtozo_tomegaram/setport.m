function [out,tomegaram]=setport(flag,tomegaram,varargin)

% flag beállításai:
% -1: inicializálás, csak az out vektr hossza kell
%  0: stacioner futtatas
%  1: instacioner futtatas, dt=varargin{1}

ro  = tomegaram.tranziens_agelem_1csp.ro;
csp = tomegaram.tranziens_agelem_1csp.csp;

switch flag
 case -1
  out{1}=1;
 case 0
  out{1}={1,0,0,tomegaram.m,[csp,0.0]};
 case 1
%  xr = varargin{1};
  t  = varargin{2}{1}{1}(1);
  dt = varargin{2}{1}{1}(2);
  if t > tomegaram.tt(end)
      tmod = mod(t,tomegaram.tt(end));
      tomegaram.m = interp1(tomegaram.tt,tomegaram.mm,tmod);
  else
      tomegaram.m = interp1(tomegaram.tt,tomegaram.mm,t);
  end
  %out{1} = {0,0,0,tomegaram.Q,[csp,-1e5]};
  out{1} = {1,0,0,tomegaram.m,[csp,0.0]};
 case 3
  out{1}=1;
 otherwise
  error('Ismeretlen opcio');
end