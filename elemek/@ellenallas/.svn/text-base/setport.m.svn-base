function [out,ellenallas]=setport(flag,ellenallas,varargin)

ro  = ellenallas.tranziens_agelem_2csp.ro;
csp = ellenallas.tranziens_agelem_2csp.csp;

switch flag
    case -1 % inicializálás
        out{1}=1;
        
    case 0 % stacioner számítás
        out{1}={ellenallas.K1/ro, 0, 0, ellenallas.K0, [csp(1),-1e5] [csp(2),1e5]};
        
    case 1 % tranziens számítás
        dt  = varargin{1};
        out{1}={ellenallas.K1/ro, 0, 0, ellenallas.K0, [csp(1),-1e5] [csp(2),1e5]};
   case 3
        out{1}=1;
	
    otherwise
        error('Ajjaj...');
end
