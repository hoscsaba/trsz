function [out,fojt]=setport(flag,fojt,varargin)

ro  = fojt.tranziens_agelem_2csp.ro;
csp = fojt.tranziens_agelem_2csp.csp;

switch flag
    case -1 % inicializálás
        out{1}=1;
        
    case 0 % stacioner számítás
        out{1}={0, 0, fojt.K1/ro, fojt.K0, [csp(1),-1e5] [csp(2),1e5]};
        
    case 1 % tranziens számítás
        dt  = varargin{1};
        out{1}={0, 0, fojt.K1/ro, fojt.K0, [csp(1),-1e5] [csp(2),1e5]};
   case 3
        out{1}=1;
	
    otherwise
        error('Ajjaj...');
end
