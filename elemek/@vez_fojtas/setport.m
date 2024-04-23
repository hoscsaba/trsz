function [out,vfojt]=setport(flag,vfojt,varargin)

ro  = vfojt.tranziens_agelem_2csp.ro;
csp = vfojt.tranziens_agelem_2csp.csp;

switch flag
    case -1 % inicializálás
        out{1}=1;
        
    case 0 % stacioner számítás
        epD = interp1(fojt.t,fojt.t_epD,t);
        K=interp1(fojt.K_epD,fojt.K,epD);
        KK=((1-K)/K)^2/2/fojt.A^2;
        out{1}={0, 0, fojt.K1*ro, fojt.K0, [csp(1),-1] [csp(2),1]};
        
    case 1 % tranziens számítás

        xr = varargin{1};
        t  = varargin{2}{1}{1}(1);
        dt = varargin{2}{1}{1}(2);
        
        epD = interp1(vfojt.t_vf,vfojt.epst,t);
        K   = interp1(vfojt.epsK,vfojt.K,epD);
%        if K<1e-3, K=1e-3; end
        KK=K;%((1-K)/K)^2/2/fojt.A^2;
        
%        fprintf('\n t=%5.3e  epD=%5.3e  K=%5.3e  KK=%5.3e',t,epD,K,KK);

        %out{1}={0, 0, KK*ro, 0, [csp(1),-1] [csp(2),1]};
        out{1}={0, 0, KK/ro, 0, [csp(1),-1e5] [csp(2),1e5]};
        
 case 3
   out{1}=1;
  
    otherwise
        error('Ajjaj...');
end
