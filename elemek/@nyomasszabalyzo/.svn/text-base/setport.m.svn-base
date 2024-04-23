function [out,nysz]=setport(flag,nysz,varargin)

debug=0;
if debug>1
    fprintf('\n\n--------------------------------------------');
    fprintf('\n ---> Nyomásszab. %s',nysz.tranziens_agelem_2csp.nev);
end

g   = 9.81;
ro  = nysz.tranziens_agelem_2csp.ro;
csp = nysz.tranziens_agelem_2csp.csp;
%xr  = sziv.tranziens_agelem_2csp.Qr;

if flag==-1 % inicializálás
    out{1}=1;
    
elseif flag==3
    
    % Ide jon a belso valtozok atallitasa, pl. szivattyú kifutás
    % esetén a núj
    xr=varargin{1};
    t  = varargin{2}(1);
    dt = varargin{2}(2);
    mar= varargin{3};
    switch nysz.szab
        case 'p'
            pp=get(mar,'p',nysz.szcsp1);
        case 'dp'
            pp=get(mar,'p',nysz.szcsp1)-get(mar,'p',nysz.szcsp2);
    end
    if debug>0
        fprintf('\n pp=%+6.4e  alapjel=%+6.4e',pp,nysz.ajel);
    end
    
    if t>nysz.tbe      
        xuj=(pp/1e5-nysz.ajel)/nysz.pmax;
        xatl=(nysz.xx(4)+xuj)/2;
        nysz.ii=nysz.ii+dt*xatl;  
        epsz=nysz.szabP*xatl + nysz.szabI*nysz.ii + nysz.szabD*(xuj-nysz.xx(4))/dt;
        if abs(epsz)>dt*nysz.vmax
            epsz=dt*nysz.vmax*sign(epsz); 
            if debug>1
                fprintf('\n  --> sebességkorlát: epsz=%+5.3e',epsz);
            end
        end        
        
        d1=nysz.e+epsz-nysz.er;
        d2=nysz.er-nysz.err;
        if debug>0
            fprintf('\n t=%5.3e  d1=%+5.3e  d2=%+5.3e  e=%+5.3e  epsz=%+5.3e  e+epsz=%+5.3e',t,d1,d2,nysz.e,epsz,nysz.e+epsz);
        end
        % Alulrelaxáció a lengések elkerülése miatt
        if d1*d2<0 & 0.99*abs(d1)<abs(d2)
            nysz.e=(nysz.er+nysz.err)/2;        
            if debug>1
                fprintf('\n  --> alulrelaxacio: e=%+5.3e',nysz.e);
            end
        else
            nysz.e=nysz.e+epsz;
        end      
    else
        xuj=nysz.xx(4);
    end
    
    if nysz.e<0, nysz.e=0; end
    if nysz.e>1, nysz.e=1; end
    for i=1:3
        nysz.xx(i)=nysz.xx(i+1);
        nysz.tt(i)=nysz.tt(i+1);
    end
    nysz.xx(4)=xuj; nysz.tt(4)=t;
    nysz.err=nysz.er;
    nysz.er =nysz.e;
    
    out{1}=1;
    
else
    
    xr = varargin{1};
    t  = varargin{2}{1}{1}(1);
    dt = varargin{2}{1}{1}(2);
    if t>nysz.tbe
        xip1=interp1(nysz.tt,nysz.xx,t,'cubic','extrap');
        xatl=(nysz.xx(4)+xip1)/2;
        epsz=nysz.szabP*xatl + nysz.szabI*(nysz.ii+dt*xatl) + nysz.szabD*(xip1-nysz.xx(4))/dt;
        
        if abs(epsz)>dt*nysz.vmax, epsz=dt*nysz.vmax*sign(epsz); end
        
        ee = nysz.e+epsz;        
    else
        ee=nysz.e;
    end
    
    if ee<0, ee=0; end
    if ee>1, ee=1; end
    
    K=interp1(nysz.jge,nysz.jgk,ee);
    if K<1e-10, K=1e-10; end
    KK=((1-K)/K)^2/2/nysz.A^2;
    if debug>1
        fprintf('\n t=%5.3e  epsz=%5.3e  ee=%5.3e  X=%+5.3e  K=%5.3e  dzeta=%5.3e',t,epsz,ee,xip1,K,KK);
    end
    nysz.euj=ee;
    out{1}={0, 0, KK*ro, 0, [csp(1),-1] [csp(2),1]};
    
end
