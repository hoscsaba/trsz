function [out,buko]=setport(flag,buko,varargin)

ro  = buko.tranziens_agelem_2csp.ro;
csp = buko.tranziens_agelem_2csp.csp;
Qr  = buko.tranziens_agelem_2csp.Q;

switch flag
    case -1 % inicializalas
        out{1}=1;
        
    case 0 % stacioner szamitas
        error('Stacionarius szamitas!!!!!');
        
    case 1 % tranziens szamitas
        xr = varargin{1};
        t  = varargin{2}{1}{1}(1);
        dt = varargin{2}{1}{1}(2);
        pe = get(varargin{3},'p',csp(1))-1e5;
        pv = get(varargin{3},'p',csp(2))-1e5;
        he=pe/ro/9.81;
        hv=pv/ro/9.81;

        if he<buko.h0
            % nagy ellenallasu fojtas
            K=1e10;            
            out{1}={K/ro, 0, 0, 0, [csp(1),-1e5] [csp(2),1e5]};
        else
            dhr=(he-buko.h0)^(buko.kitevo-1);
            szorzo1=buko.Cd*ro*buko.B*sqrt(2*9.81)*dhr;
            szorzo=szorzo1*(1e5/ro/9.81);
            konst=szorzo1*(buko.h0+1e5/ro/9.81);
            out{1}={1, 0, 0, konst, [csp(1),-szorzo] [csp(2),0]};        
        end
   case 3
       xr = varargin{1};
       t  = varargin{2}(1);
       pe = get(varargin{3},'p',csp(1))-1e5;
       pv = get(varargin{3},'p',csp(2))-1e5;
       he=pe/ro/9.81;
       hv=pv/ro/9.81;
       
       if (buko.is_working==0) && (he>buko.h0)
           fprintf('\n\t t=%g s, a %s buko mukodesbe lepett: he=%5.3f m, hv=%5.3f m, dh=%5.3f, Q=%5.3f m3/h',t,buko.tranziens_agelem_2csp.nev,he,hv, he-buko.h0,xr*3.600);
           fprintf(' (%5.3f m3/h)',buko.Cd*ro*buko.B*sqrt(2*9.81)*(he-buko.h0)^buko.kitevo*3.6);
           buko.is_working=1;
       end
       if (buko.is_working==1) && (he<buko.h0)
           fprintf('\n\t t=%g s, a %s buko bezart',t,buko.tranziens_agelem_2csp.nev);
           buko.is_working=0;
       end
       out{1}=1;
       
    otherwise
        error('Ilyen opcio nincs; 1 (NR iteracio) vagy 3 (iteracio vegen update)');
end
