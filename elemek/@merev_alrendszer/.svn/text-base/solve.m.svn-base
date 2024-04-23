function mar = solve(mar,pf,varargin)

t_ide = varargin{1}(1);
t_innen = mar.t;
varargin{1}(2) = t_ide-t_innen;
dt = t_ide-t_innen;

% for i = 1:length(mar.elemek)
%     if isa(mar.elemek{i},'akna')
%         mar.elemek{i}.nev
%         mar.elemek{i}.Q
%     end
% end

%% Default hibabeallitasok
e_ag_max    =1e-3; % Agegyenletek, [bar]
e_konti_max =1e-3; % konti, [kg/s]
e_rug_max   =1e-3; % rug. csatl., [bar] 
imax=300; textout=0;

%% Egyeb beallitasok
if nargin==3
  flag=1; % tranziens futas
else
  flag=0; % stacioner futas
end

if ~(nargin==2)
    property_argin = varargin;
    while length(property_argin) >= 2,
        prop = property_argin{1};
        val  = property_argin{2};
        property_argin = property_argin(3:end);
        switch prop
            case 'emax', emax=val;
            case 'imax', imax=val;
            case 'textout', textout=val;
            otherwise, error(['Ismeretlen opcio:',prop]);
        end
    end
end

%% Adatok attoltese szamitasokhoz
n_csp=length(mar.csp); n_elem=length(mar.elemek);
for i=1:length(mar.elemek)
    temp=mar.elemek{i};  xr(i)=temp.Q*temp.ro;    
end
for i=1:length(mar.csp)
    xr(n_elem+i)=mar.csp{i}{4}/1e5;    
end
xr = xr';

%% Ha rug_cso vagy csatorna pf van, hozza kell adni uj terfogataramot is..
for i=1:length(pf)
    if length(strmatch(pf{i}{2},'cso','exact'))==1 ||...
            length(strmatch(pf{i}{2},'vizszint_&_konti','exact'))==1 || ...
            length(strmatch(pf{i}{2},'viszkcso','exact'))==1 % Itt is beleny�ltam
        xr=[xr; 0];
    end
end
xrr=xr;

%% Fogyasztasok frissitese aktualis idoszintre
mar=update_csp_fogy(mar,t_ide);

%% Es indul a madula

if textout==1, fprintf('\nSzamitas....'); end

i=1; hiba=[1e10, 1e10, 1e10];
if mar.save_level>5
  [A,B,C,D,mar] = agmatrix(mar,pf,flag,xr,varargin);
  save_abcd(mar.fnev,mar.elemek,mar.csp,A,B,C,D,'x',mar.lepes); 
end
if mar.save_level>2, save_iter(mar.fnev,1,[e_ag_max,e_konti_max,e_rug_max],imax,mar.lepes,varargin);  end
%plot_it_null=mar.plot_it(length(mar.plot_it(:,1)),1);

mar.iter=1;
while ((hiba(1)>e_ag_max) || (hiba(2)>e_konti_max) || (hiba(3)>e_rug_max)) && (i<imax+1)    
  [A,B,C,D,mar] = agmatrix(mar,pf,flag,xr,varargin);
  [hiba,xu] = nr(A,B,C,D,xr,n_csp,n_elem);   
  if mar.save_level>5, save_abcd(strcat(mar.nev,'_ABCD.out'),mar.elemek,mar.csp,A,B,C,D,'val',xr,xu,mar.lepes,mar.iter); end
  if textout==1, fprintf('\n%3d./%3d lepes    rms agegy = %7.5e   rms konti = %7.5e  rms rugcs = %7.5e',i,imax,hiba(1),hiba(2),hiba(3)); end
  if mar.save_level>2, save_iter(mar.fnev,2,i,imax,hiba); end
  %mar.plot_it=[mar.plot_it; plot_it_null+i norm(hiba)];
  
  xr=xu; i=i+1; mar.iter=mar.iter+1;
end

if i==imax+1
    uzenet='\n\n@merev_alrendszer/solve.m: A maximalis iteracioszamot elertem, valami nem stimmel...\n';    
    fprintf('\n\n %s merev alreendszer utolso lepese:\n',mar.nev);
    xu
    hiba
    error('@merev_alrendszer/solve.m: A maximalis iteracioszamot elertem!!!!');
else
    uzenet='\n\nNormal befejezes :)\n'; 
end
if textout==1, fprintf(uzenet); end
if mar.save_level>2, save_iter(mar.fnev,3,uzenet); end

% Eredmenyek visszatoltese az elemekbe
mar.t=t_ide;%mar.t+varargin{1}(2);
mar.lepes=mar.lepes+1;
for i=1:length(mar.elemek)
    mar.elemek{i}.Q = xu(i)/mar.elemek{i}.ro;

    % itt az egyeb valtozokat toltjuk vissza (pl. n_sziv)

    %   if isa(mar.elemek{i},'akna')
    %       fprintf('\n ------------------\n');
    %       mar.elemek{i}.ro
    %       mar.elemek{i}.Q
    %       pause
    %       xu
    %   end
%     if strcmp(mar.nev,'alrendsz2')
%         xu
%         pause
%     end
    [out,mar.elemek{i}] = setport(3,mar.elemek{i},mar.elemek{i}.Q,varargin{1},mar);
end

for i=1:length(mar.csp)
  mar.csp{i}{4} = xu(n_elem+i)*1e5;
  
  for j=1:length(mar.elemek)
    %        pp  =mar.elemek{j}.p*1e5;
    cspp=mar.elemek{j}.csp;
    for k=1:length(mar.elemek{j}.csp)
                 % fprintf('\n %s\t:i=%d  j=%d  cspp(%d)=%d ?=? %d',mar.elemek{j}.nev,i,j,k,cspp(k),i);
      if cspp(k)==i, pp(k)=xu(n_elem+i)*1e5; end
      %            fprintf(' -> pp(%d)=%+5.3e',k,xu(n_elem+i)*1e5);
    end            
  end
  mar.elemek{j}.p=pp;        
end

if textout==1
    fprintf('\n------------------------------------------------------------\n');
    fprintf('\nTerfogataramok es tomegaramok:');
    nQ=length(mar.elemek);
    for i=1:length(mar.elemek)
        fprintf('\nQ%2d = %+5.3e [m^3/s] = %+6.2f [l/s]        m%2d = %+5.3e [kg/s]',i,mar.elemek{i}.Q,mar.elemek{i}.Q*1e3,i,mar.elemek{i}.Q*mar.elemek{i}.ro);
    end
    
    fprintf('\n\nnyomasok:')
    for i=1:length(mar.csp)
        fprintf('\np%2d = %+5.3e [Pa] = %+6.2f [vom]',mar.csp{i}{1},mar.csp{i}{4},mar.csp{i}{4}/1e3/9.81);
    end
    fprintf('\n');
    pause
end


if mar.save_level>2, save_res(mar); end

% Eredmenyfile mentese
if mar.t > mar.dtki
    mar.dtki = mar.dtki + mar.dtkiorig;
    info(mar,3,mar.t);
end

%% Newton-Raphson lepes
function [hiba,xu]=nr(A,B,C,D,xr,n_csp,n_elem)

%fprintf('\n--------------------------------------------------------\n')

F=A*xr+B*(xr.^2)+C*(xr.*abs(xr))+D;

for i=1:length(A(:,1))
    for j=1:length(A(:,1))
        DF(i,j)=A(i,j)+2*B(i,j)*xr(j)+2*C(i,j)*abs(xr(j));
    end
end
xu = xr-inv(DF)*F;
RELAX=0.2;
xu=(1-RELAX)*xr+RELAX*xu;
hiba = [norm(F(1:n_elem))/1e5
      norm(F(n_elem+1:n_elem+n_csp));
      norm(F(n_elem+n_csp+1:length(xu)))/1e5];

% hiba_q mertekegysege [bar], hiba_p pedig [kg/s], mert ezek reziduumok.

%fprintf('\n');
%for i=1:length(F)
%    fprintf('\n i=%2d  xr=%+5.3e  xu=%+5.3e  dx=%+5.3e  F=%+5.3e',i,xr(i),xu(i),xr(i)-xu(i),F(i));
%    if i==n_elem | i==n_elem+n_csp, fprintf('\n---------------------------------------------------------'); end
%end
%pause

%% Fogyasztasok frissitese
function mar = update_csp_fogy(mar,t)

for i=1:length(mar.gorbek)
    gorbe_nevsor{i}=mar.gorbek{i}.nev;
end

for i=1:length(mar.csp)
    tmp=strcmp(mar.csp{i}{7},gorbe_nevsor);
    if sum(tmp)==0
        error(['@mar, solve.m, update_csp_fogy: nem talalom a ',mar.csp{i}{7},' nevu lefutasgorbet.']);
    else
        ezaz=1;
        while (tmp(ezaz)==0)
            ezaz=ezaz+1;
        end
        %----------------------------------------------------
        % FIXME
        % Hos Csaba, 2009. aug. 14.
        % Mi a francert kap idonkent komplex t-t? -> abs(t)
        %----------------------------------------------------
        mar.csp{i}{5}=mar.csp{i}{8}*interp1(mar.gorbek{ezaz}.x,mar.gorbek{ezaz}.y,abs(t));
    end
end




