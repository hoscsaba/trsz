function [out,cso]=setport(flag,cso,varargin)

ro  = cso.tranziens_agelem_2csp.ro;
csp = cso.tranziens_agelem_2csp.csp;
Qr  = cso.tranziens_agelem_2csp.Q;
L=cso.L; D=cso.D; A=cso.A; lam=cso.lambda;

switch flag
    case -1 % inicializálás
        out{1}=1;
        
    case 0 % stacioner számítás
        out{1}={0, 0, lam*L/D/2/A^2/ro, 0, [csp(1),-1e5] [csp(2),1e5]};
        
    case 1 % tranziens számítás
%        xr = varargin{1};
	t  = varargin{2}{1}{1}(1);
	dt = varargin{2}{1}{1}(2);
	
        out{1}={L/A/dt, 0, lam*L/D/2/A^2/ro, -L/A/dt*Qr*ro, [csp(1),-1e5] [csp(2),1e5]};
   case 3
        out{1}=1;
	
    otherwise
        error('Ilyen opcio nincs; vagy stac vagy instac!');
end
