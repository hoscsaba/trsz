function [out,legust]=setport(flag,legust,varargin)
  
ro  = legust.tranziens_agelem_1csp.ro;
csp = legust.tranziens_agelem_1csp.csp;

if flag == -1
  out{1}=1;
  
elseif flag==3
  xr = varargin{1};
  t  = varargin{2}(1);
  dt = varargin{2}(2);
  
  Q = xr;
  
  Vuj = legust.V - Q*dt;
  if Vuj>legust.A*legust.H, error('A légüst kiürült!!!');  end
  puj=legust.p + legust.n*legust.p/Vuj*Q*dt;
  
  hiba = abs(puj*Vuj^legust.n-legust.p0*legust.V0^legust.n)/(legust.p0*legust.V0^legust.n)*100;
  %  if hiba>10
  %    warning('\n A %s légüstben p*V^n értéke %4.2f %%-al eltér a kezdeti állapottól!',legust.tranziens_agelem_1csp.nev,hiba);
  %  end
  
%  if strcmp(legust.tranziens_agelem_1csp.nev,'legust1')
%    fprintf('\n --> t=%5.3e dt=%5.3e  Q=%+5.3e | V=%5.3e -> %5.3e | p=%5.3e -> %5.3e |  hiba p*V^n [%%]= %6.3f',t,dt,Q,legust.V,Vuj,legust.p/1e5,puj/1e5,hiba);
%    pause
%  end

  legust.dtee=legust.dte;  legust.dte=dt;
  legust.Vrr = legust.Vr; legust.Vr = legust.V; legust.V=Vuj;
  legust.prr = legust.pr; legust.pr = legust.p; legust.p=puj;
  out=1;
  
else
  %  xr = varargin{1};
  t  = varargin{2}{1}{1}(1);
  dt = varargin{2}{1}{1}(2);
  
  switch legust.extrap
   case 0
    pp=legust.p; VV=legust.V; 
   case 1   
    % linearis extrapolacio a kovetkezo idoszintre
    pp=legust.pr+(legust.p-legust.pr)/legust.dte*(legust.dte+dt/2);
    VV=legust.Vr+(legust.V-legust.Vr)/legust.dte*(legust.dte+dt/2);
   case 2
    pp=interp1([0, legust.dtee, legust.dtee+legust.dte],[legust.prr, legust.pr, legust.p], legust.dtee+legust.dte+dt/2,'spline','extrap');
    VV=interp1([0, legust.dtee, legust.dtee+legust.dte],[legust.Vrr, legust.Vr, legust.V], legust.dtee+legust.dte+dt/2,'spline','extrap');    
  end
  
  out{1}={dt*(legust.n*pp/VV-9.81*ro/legust.A)/ro,0,0,legust.p,[csp,-1e5]};
%  if strcmp(legust.tranziens_agelem_1csp.nev,'legust1')
%    fprintf('\n     dt=%5.3e  pp=%+5.3e  VV=%5.3e alfa=%5.3e ',dt,pp,VV,dt*(legust.n*pp/VV-9.81*ro/legust.A)/ro);
%    pause
%  end
  
end
