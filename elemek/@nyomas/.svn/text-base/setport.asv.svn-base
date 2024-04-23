function [out,nyomas]=setport(flag,nyomas,varargin)

% flag beállításai:
% -1: inicializálás, csak az out vektr hossza kell
%  0: stacioner futtatas
%  1: instacioner futtatas, dt=varargin{1}

ro  = nyomas.tranziens_agelem_1csp.ro;
csp = nyomas.tranziens_agelem_1csp.csp;

switch flag
    case -1
        out{1}=1;
    case 0
        out{1}={0,0,0,nyomas.p,[csp,-1e5]};
    case 1
        dt=varargin{1};
        out{1}={0,0,0,nyomas.p,[csp,-1e5]};
    case 3
        out{1}=1;
    otherwise
        error('Ismeretlen opcio');
end
