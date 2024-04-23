function out = olvas8m(fnev,debug,wdir)

if debug>1
    fprintf('\n\n****************************************************');
    fprintf('\n*          MEREV ALRENDSZEREK OLVASASA             *');
    fprintf('\n****************************************************\n\n');
end

if debug>0
    fprintf('\n\n Feladat neve: %s, adatok olvasasa...\n',fnev);
end

fnevdat = strcat(fnev,'.mdt');
[be] = textread(fnevdat,'%s','delimiter',',');

if debug>0
    fprintf('\n Merev alrendszerek szama : %g',sum(strcmp(be,'mar')));
    fprintf('\n Rugalmas elemek:');
    fprintf('\n   nyomottvizes cso       : %g',sum(strcmp(be,'rug_cso')));
    fprintf('\n   csatorna               : %g',sum(strcmp(be,'csatorna')));
    fprintf('\n   viszkoelasztikus cso   : %g',sum(strcmp(be,'viszkcso')));
    fprintf('\n Gorbek szama             : %g',sum(strcmp(be,'gorbe')));
    fprintf('\n Opciok szama             : %g',sum(strcmp(be,'option')));
end


s_mar = struct(); % mar
s_option = struct(); % option
s_gorbe = {}; % gorbek

while ~isempty(be)
    tipus=char(be(1));
    switch tipus
        case 'mar'
            [be,s_mar] = get_mar_data(be,s_mar,debug);
            
        case 'option'
            [be,s_option] = get_option_data(be,s_option,debug);
            
        case 'gorbe'
            [be,s_gorbe] = get_gorbe_data(be,s_gorbe,debug);
            
        otherwise
            disp(be);
            error('Ismeretlen merev alrendszer elem tipus: %s (mar|option)',tipus);
    end    
end

for i=1:length(s_mar.name)
    out.mar{i}=merev_alrendszer(s_mar.name{i},s_mar.elemek{i},s_mar.cspok{i},wdir);   
    out.mar{i}.gorbek=s_gorbe;
end
out.options=s_option;

%----------------------------------------------------------------------
function [be,option_data] = get_option_data(be,option_data,debug)

% Az elem sorszama
if isfield(option_data,'name')
    no = length(option_data.name)+1;
else
    no=1;
end

if (length(be)<3)
    error('HIBA: option (%s) olvasasa, nincs eleg adat: %g van de min. %g kellene!',...
        be{2},length(be),3);
end

option_data.name{no}  = char(be(2));
option_data.value{no} = char(be(3));

if debug>0
    fprintf('\n\n OPTION:');
    fprintf('\n   nev   : %s',option_data.name{no});   
    fprintf('\n   ertek : %s',option_data.value{no});
end

be=be(4:end);

%----------------------------------------------------------------------
function [be,gorbe_data] = get_gorbe_data(be,gorbe_data,debug)
 
% Az elem sorszama
if ~isempty(gorbe_data)
    no = length(gorbe_data)+1;
    nev   = char(be(2));
    jgpsz = str2double(char(be(3)));
    x=zeros(1,jgpsz); y=x;
    for j=1:jgpsz
        x(j) = str2double(char(be(3+2*j-1)));
        y(j) = str2double(char(be(3+2*j)));
    end
    gorbe_data{no}=gorbe(nev,x,y);    
    be=be(4+2*jgpsz:end);
else
    no=1;
    nev='const';
    gorbe_data{1}=gorbe('const',[0 1e10],[1 1]);
end

if debug>0
    fprintf('\n\n GORBE:');
    fprintf('\n   %2d. gorbe (%s)',no,nev); 
end

%----------------------------------------------------------------------
function [be,mar_data] = get_mar_data(be,mar_data,debug)

% Az elem sorszama
if isfield(mar_data,'name')
    mar_no = length(mar_data.name)+1;
else
    mar_no = 1;
end

mar_data.name{mar_no} = char(be(2));
be=be(3:end);
stop = 0;
elem_no=1;
csp_no =1;

if debug>1
    fprintf('\n\n MEREV ALRENDSZER: %s',mar_data.name{mar_no});
end

while ~stop
    tipus=char(be(1));
    switch tipus
        case {'mar','gorbe','option'}
            stop=1;
           
        case 'csp'
            data = get_csp_data(be,csp_no,debug);
                        
        case 'konc_cso'
            data = get_konc_cso_data(be,elem_no,debug);
            
        case 'fojtas' 
            data = get_fojtas_data(be,elem_no,debug);
        
        case 'akna'
            data = get_akna_data(be,elem_no,debug);
            
        case 'nyomovezetek'
            data = get_nyomovezetek_data(be,elem_no,debug);
 
        case 'szivattyu'
            data = get_szivattyu_data(be,elem_no,debug);
            
        otherwise
            error('Ismeretlen merev alrendszer elem tipus: %s (neve:%s)',...
                char(be(1)),char(be(2)));
    end
    
    if ~stop
        be = data.be;
        if ~strcmp(tipus,'csp')
            cspe_nevek{elem_no} = data.cspe_nev;
            cspv_nevek{elem_no} = data.cspv_nev;
            elemek{elem_no}     = data.elem;
            elem_no             = elem_no+1;
        else
            csp{csp_no} = data.csp;
            csp_no = csp_no+1;
        end
    end
end

if debug>2
    fprintf('\n\n  Elem csompontazonositok keresese:');
end
    
cspe_num=0;
cspv_num=0;
for i=1:length(elemek)
    for j=1:length(csp)
        if strcmp(cspe_nevek{i},csp{j}{6})
            cspe_num = j;
        end

        if ~isempty(cspv_nevek{i}) 
            if strcmp(cspv_nevek{i},csp{j}{6})
                cspv_num = j;
            end
        end
    end
    
    if cspe_num==0
        error([elemek{i}.nev,' agelem eleje csp. (',cspe_nevek{i},') nincs meg!!!']);
    end
    if ~isempty(cspv_nevek{i}) && cspv_num==0
        error([elemek{i}.nev,' agelem vege csp. (',cspe_nevek{i},') nincs meg!!!']);
    end
        
    if ~isempty(cspv_nevek{i})
        elemek{i}.csp = [cspe_num cspv_num];
    else
        elemek{i}.csp = cspe_num;
    end    
end

if debug>2
    for i=1:length(elemek)
        cspi=elemek{i}.csp;
        fprintf('\n\t%s:\t %s',elemek{i}.nev,csp{cspi(1)}{6});
      if ~isempty(cspv_nevek{i})
          fprintf(' -> %s',csp{cspi(2)}{6});
      end
    end
end

mar_data.elemek{mar_no} = elemek;
mar_data.cspok{mar_no} = csp;

%-----------------------------------------------------------------------
function out = get_csp_data(be,csp_no,debug)

% Sorszam
out.csp{1} = csp_no;
% z [m]
out.csp{2} = str2double(char(be(3)));
% agazonositok
out.csp{3} = [];
% kezdeti nyomas
out.csp{4} = 1e5;
% pillanatnyi fogyasztas, ezt fogjuk szorozni a lefutas pillanatnyi ertekevel.
out.csp{5} = str2double(char(be(4)))/3600*1000;
% nev
out.csp{6} = char(be(2));
% lefutas gorbe azonositoja
out.csp{7} = char(be(5));
% Nevleges fogyasztas
out.csp{8} = str2double(char(be(4)))/3600*1000;

if debug>1
    fprintf('\n   %2d. csomopont: %s',csp_no,char(be(2)));
end

out.be = be(6:end);

%-----------------------------------------------------------------------
function out = get_konc_cso_data(be,elem_no,debug)

agnev        = char(be(2));
out.cspe_nev = char(be(3));
out.cspv_nev = char(be(4));
ro  = str2double(char(be(5))); % suruseg [kg/m3]
m0  = str2double(char(be(6))); % kezdeti tomegaram [kg/s]

D   = str2double(char(be(7))); % atmero [m]
L   = str2double(char(be(8))); % hossz [m]
lam = str2double(char(be(9))); % lambda [-]

out.elem=konc_cso(agnev,0,0,D,L,lam,ro,m0);
out.be  =be(10:end);

if debug>1
    fprintf('\n   %2d. elem: konc_cso (%s)',elem_no,agnev);
end

%-----------------------------------------------------------------------
function out = get_akna_data(be,elem_no,debug)

agnev        = char(be(2));
out.cspe_nev = char(be(3));
out.cspv_nev = '';
ro   = str2double(char(be(4))); % suruseg [kg/m3]
m0   = str2double(char(be(5))); % kezdeti tomegaram [kg/s]
            
AA   = str2double(char(be(6)));
hmin = str2double(char(be(7)));
hmax = str2double(char(be(8)));
yy0  = str2double(char(be(9)));
rajz = char(be(10));

out.elem = akna(agnev,0,AA,hmin,hmax,yy0,ro,m0,rajz);
out.be   = be(11:end);

if debug>1
    fprintf('\n   %2d. elem: akna (%s)',elem_no,agnev);
end

%-----------------------------------------------------------------------
function out = get_nyomovezetek_data(be,elem_no,debug)

agnev          = char(be(2));
out.cspe_nev = char(be(3));
out.cspv_nev = char(be(4));
ro = str2double(char(be(5))); % suruseg [kg/m3]
m0 = str2double(char(be(6))); % kezdeti tomegaram [kg/s]

D   = str2double(char(be(7)));  % atmero [m]
L   = str2double(char(be(8)));  % hossz [m]
lam = str2double(char(be(9)));  % lambda [-]
ze  = str2double(char(be(10))); % csoeleje magassag [m]
zv  = str2double(char(be(11))); % cso vege magassag [m]

out.elem = nyomovezetek(agnev,0,0,D,L,lam,ze,zv,ro,m0);
out.be = be(12:end);

if debug>1
    fprintf('\n   %2d. elem: nyomovezetek (%s)',elem_no,agnev);
end

%-----------------------------------------------------------------------
function out = get_fojtas_data(be,elem_no,debug)
agnev      = char(be(2));
out.cspe_nev = char(be(3));
out.cspv_nev = char(be(4));
ro = str2double(char(be(5))); % suruseg [kg/m3]
m0 = str2double(char(be(6))); % kezdeti tomegaram [kg/s]

K1 = str2double(char(be(7))); % ellenallas tenyezo [1/m2]

out.elem = fojtas(agnev,0,0,1e-10,K1,ro,m0);
out.be = be(8:end);

if debug>1
    fprintf('\n   %2d. elem: fojtas (%s)',elem_no,agnev);
end

%-----------------------------------------------------------------------
function out = get_szivattyu_data(be,elem_no,debug)

agnev               = char(be(2));
out.cspe_nev = char(be(3));
out.cspv_nev = char(be(4));
ro = str2double(char(be(5))); % suruseg [kg/m3]
m0 = str2double(char(be(6))); % kezdeti tomegaram [kg/s]

Ds        = str2double(char(be(7)));  % D szivocsonk [m]
Dn        = str2double(char(be(8)));  % D nyomocsonk [m]
tranziens = str2double(char(be(9)));  % tranziens tipus
jgpsz     = str2double(char(be(10))); % H(Q) jg. pontszam

be=be(11:end);

Q=zeros(jgpsz,1); H=Q; P=Q;
for j=1:jgpsz
    jj=3*j-2;
    Q=str2double(char(be(jj)));   % Q [m3/s]
    H=str2double(char(be(jj+1))); % H [m]
    P=str2double(char(be(jj+2))); % P [kW]
end

be=be(3*jgpsz+1:end);

switch tranziens
    case 0  % allando fordulatszam
        out.elem = szivattyu(agnev,0,0,Q,H,Ds,Dn,ro,m0,tranziens);
        
    case 1     % kifutas
        n    = str2double(char(be(1))); % fordulatszam [1/min]
        teta = str2double(char(be(2))); % tehetetlensegi nyomatek [kgm2]
        tki  = str2double(char(be(3))); % kifutas idopontja [s]
        out.elem = szivattyu(agnev,0,0,Q,H,Ds,Dn,ro,m0,tranziens,...
            P,n,teta,tki);
        
        be=be(4:end);
        
    case 2     % inditas, kozelito motor-jelleggorbe
        n    = str2double(char(be(1))); % fordulatszam [1/min]
        teta = str2double(char(be(2))); % tehetetlensegi nyomatek [kgm2]
        Mb   = str2double(char(be(3))); % billenonyomatek [Nm]
        nb   = str2double(char(be(4))); % billenonyomatek fordulatszama [1/min]
        nsz  = str2double(char(be(5))); % szinkronfordulatszama [1/min]
        out.elem = szivattyu(agnev,0,0,Q,H,Ds,Dn,ro,m0,tranziens,...
            P,n,teta,Mb,nb,nsz);
        
        be=be(6:end);
        
    case 3     % inditas, adott motor-jelleggorbe
        n     = str2double(char(be(1))); % fordulatszam [1/min]
        teta  = str2(doublechar(be(2))); % tehetetlensegi nyomatek [kgm2]
        jgpsz = str2double(char(be(3))); % motor jg. psz.
        be=be(4:end);
        Mm=zeros(jgpsz,1); nm=Mm;
        for j=1:jgpsz
            jj=2*j-2;
            nm = str2double(char(be(jj)));   % fordulatszam [1/min]
            Mm = str2double(char(be(jj+1))); % nyomatek [Nm]
        end
        out.elem = szivattyu(agnev,0,0,Q,H,Ds,Dn,ro,m0,tranziens,...
            P,n,teta,Mm,nm);
        
        be=be(2*jgpsz+1:end);
        
    case 4     % frekivaltos inditas, kozelito motor-jelleggorbe
        n    = str2double(char(be(1))); % fordulatszam [1/min]
        teta = str2double(char(be(2))); % tehetetlensegi nyomatek [kgm2]
        Mb   = str2double(char(be(3))); % billeno nyomatek [Nm]
        nb   = str2double(char(be(4))); % billeno-nyomatek fordulatszama [1/min]
        Mi   = str2double(char(be(5))); % indítonyomatek [Nm]
        nsz  = str2double(char(be(6))); % szinkronfordulatszam [1/min]
        tbe  = str2double(char(be(7))); % inditas időpontja [s]
        tind = str2double(char(be(8))); % felfutas idotartama [s]
        out.elem = szivattyu(agnev,0,0,Q,H,Ds,Dn,ro,m0,tranziens,...
            P,n,teta,Mb,nb,Mi,nsz,tbe,tind);
        
        be=be(9:end);
        
    case 5
        n    = str2double(char(be(1))); % fordulatszam [1/min]
        hbe  = str2double(char(be(2))); % bekapcsolasi szint [m]
        hki  = str2double(char(be(3))); % kikapcsolasi szint [m]
        uzem = str2double(char(be(4))); % szivattyu uzemel (0/1)
        out.elem = szivattyu(agnev,0,0,Q,H,Ds,Dn,ro,m0,tranziens,...
            P,n,hbe,hki,uzem);
        
        be=be(5:end);
        
    otherwise
        error(['Szivattyu ',agnev,': hibas tranziens jel: ',tranziens]);
end

if debug>0
    fprintf('\n   %2d. elem: szivattyu (%s)',elem_no,agnev);
end

out.be=be;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             fnev = char(be(1));
% ncp = str2num(char(be(2)));       %....beolvasando merev alrendszerek szama...
% eltol = 3;
% mar_szamlalo=0;
% for jc=1:ncp
%     ag_szamlalo=0;
%     mar_szamlalo=mar_szamlalo+1;
%     cspnev = char(be(eltol));               %....rugalmas-csomopont neve
%     nag = str2num(char(be(eltol+1)));       %....agak szama a csp-ben..
%     if debug>1, fprintf('\n\n%2d./%2d csomoponti alrendszer neve: %s  (%d ag)',mar_szamlalo,ncp,cspnev,nag); end
% 
%     %...az agadatok ....
%     eltol=eltol+2;
%     for k=1:nag
%         ag_szamlalo =ag_szamlalo+1;
%         tipus{k}=char(be(eltol));
%         switch tipus{k}
%                      
%             case 'vez_fojtas'
%                 agnev{k}=char(be(eltol+1));                   %...agnev..
%                 csp1{k}=char(be(eltol+2));                    %...csp elejen
%                 csp2{k}=char(be(eltol+3));                    %...csp vegen
%                 ro(k)=str2num(char(be(eltol+4)));             %...suruseg [kg/m3]
%                 A(k)=str2num(char(be(eltol+5)));              %...nevleges keresztmetszet [m2]
%                 jgpsz(k)=str2num(char(be(eltol+6)));          %...K(eps) jellegg. pontsz.
%                 
%                 jgpszt(k)=str2num(char(be(eltol+7)));        %...eps(t) fv.
%                 m0(k)=str2num(char(be(eltol+8)));            %...tomegaram [kg/s]
% 
%                 eltol=eltol+8;
%                 for j=1:jgpsz(k)
%                     jj=2*j-1;
%                     epsK(k,j)=str2num(char(be(eltol+jj)));             %...rel. helyzet [-]
%                     K(k,j)   =str2num(char(be(eltol+jj+1)));           %...ellenallas tenyezo [1/m2]
%                     %		    fprintf('\n\t %2d./%2d adat:  eps=%+5.3e  K=%5.3e',j,jgpsz(k),eps(k,j),K(k,j));
%                 end
%                 eltol=eltol+2*jgpsz(k);
% 
%                 for j=1:jgpszt(k)
%                     jj=2*j-1;
%                     t_vf(k,j)=str2num(char(be(eltol+jj)));               %...ido [s]
%                     epst(k,j)=str2num(char(be(eltol+jj+1)));             %...rel. helyzet [-]
%                     %                    fprintf('\n\t %2d./%2d adat:  t=%+5.3e  eps=%5.3e',j,jgpszt(k),t_vf(k,j),epst(k,j));
%                 end
%                 eltol=eltol+2*jgpszt(k)+1;
% 
%                 % Ez m�r k�z�s...
%                 if debug>1, fprintf('\n   %2d. elem: vez_fojtas (%s)',ag_szamlalo,agnev{k}); end
% 
%            
% 
%             case 'nyomas'
%                 agnev{k}=char(be(eltol+1));                 %...agnev..
%                 csp1{k}=char(be(eltol+2));                   %...csp elejen
%                 csp2{k}='';
%                 pk(k)=str2num(char(be(eltol+3)));            %%...adott nyomas [N/m2]
%                 ro(k)=str2num(char(be(eltol+4)));            %...adott nyomas [N/m2]
%                 eltol=eltol+5;
%                 if debug>1, fprintf('\n   %2d. elem: nyomas (%s)',ag_szamlalo,agnev{k}); end
% 
%             case 'valtozo_nyomas'
%                 agnev{k}=char(be(eltol+1));                 %...agnev..
%                 csp1{k}=char(be(eltol+2));                   %...csp elejen
%                 csp2{k}='';
%                 %                pk(k)=str2num(char(be(eltol+3)));            %%...adott nyomas [N/m2]
%                 ro(k) = str2num(char(be(eltol+3)));             %%...adott nyomas [N/m2]
%                 jgpsz = str2num(char(be(eltol+4)));             %...K(eps) jellegg. pontsz.
%                 eltol = eltol+4;
%                 for j=1:jgpsz
%                     jj=2*j-1;
%                     tt(k,j) = str2num(char(be(eltol+jj)));
%                     pp(k,j) = str2num(char(be(eltol+jj+1)));
%                 end
%                 eltol=eltol+2*jgpsz+1;
%                 if debug>1, fprintf('\n   %2d. elem: valtozo nyomas (%s)',ag_szamlalo,agnev{k}); end
% 
%             case 'valtozo_tomegaram'
%                 agnev{k}=char(be(eltol+1));                 %...agnev..
%                 csp1{k}=char(be(eltol+2));                   %...csp elejen
%                 csp2{k}='';
%                 %                pk(k)=str2num(char(be(eltol+3)));            %%...adott nyomas [N/m2]
%                 ro(k) = str2num(char(be(eltol+3)));             %%...adott suruseg [kg/m3]
%                 jgpsz = str2num(char(be(eltol+4)));             %...Q(t) jellegg. pontsz.
%                 eltol = eltol+4;
%                 for j=1:jgpsz
%                     jj=2*j-1;
%                     tt(k,j) = str2num(char(be(eltol+jj)));
%                     mm(k,j) = str2num(char(be(eltol+jj+1)));
%                 end
%                 eltol=eltol+2*jgpsz+1;
%                 if debug>1, fprintf('\n   %2d. elem: valtozo tomegaram (%s)',ag_szamlalo,agnev{k}); end                
%                 
%             case 'visszacsapo_szelep'
%                 agnev{k}= char(be(eltol+1));                 %...agnev..
%                 csp1{k} = char(be(eltol+2));                 %...csp elejen
%                 csp2{k} = char(be(eltol+3));
%                 ro(k)   = str2num(char(be(eltol+4)));
%                 m0(k)   = str2num(char(be(eltol+5)));
%                 eltol = eltol+6;
%                 if debug>1, fprintf('\n   %2d. elem: visszacsapo_szelep (%s)',ag_szamlalo,agnev{k}); end
% 
%             case 'nyomasszabalyzo'
%                 agnev{k}=char(be(eltol+1));                 %...agnev..
%                 csp1{k} =char(be(eltol+2));                   %...csp elejen
%                 csp2{k} =char(be(eltol+3));
%                 ro(k)   =str2num(char(be(eltol+4)));
%                 m0(k)   =str2num(char(be(eltol+5)));
%                 jgpsz(k)=str2num(char(be(eltol+6)));
%                 eltol=eltol+6;
%                 for j=1:jgpsz(k)
%                     jj=2*j-1;
%                     e(k,j)  = str2num(char(be(eltol+jj)));
%                     ek(k,j) = str2num(char(be(eltol+jj+1)));
%                 end
%                 eltol=eltol+2*jgpsz(k);
%                 szab{k} = char(be(eltol+1));
%                 switch szab{k}
%                     case 'p'
%                         szabcsp{k} =      char(be(eltol+2));
%                         ajel(k) = str2num(char(be(eltol+3)));
%                         pmax(k) = str2num(char(be(eltol+4)));
%                         szabP(k)= str2num(char(be(eltol+5)));
%                         szabI(k)= str2num(char(be(eltol+6)));
%                         szabD(k)= str2num(char(be(eltol+7)));
%                         e0(k)   = str2num(char(be(eltol+8)));
%                         vmax(k) = str2num(char(be(eltol+9)));
%                         tbe(k)  = str2num(char(be(eltol+10)));
%                         eltol = eltol+11;
%                     case 'dp'
%                         szabcsp1{k} =     char(be(eltol+2));
%                         szabcsp2{k} =     char(be(eltol+3));
%                         ajel(k)  = str2num(char(be(eltol+4)));
%                         pmax(k)  = str2num(char(be(eltol+5)));
%                         szabP(k) = str2num(char(be(eltol+6)));
%                         szabI(k) = str2num(char(be(eltol+7)));
%                         szabD(k) = str2num(char(be(eltol+8)));
%                         e0(k)    = str2num(char(be(eltol+9)));
%                         vmax(k)  = str2num(char(be(eltol+10)));
%                         tbe(k)   = str2num(char(be(eltol+11)));
%                         eltol = eltol+12;
%                     otherwise
%                         error('Szab�lyz�szelep csak nyom�sra  vagy nyom�sk�l�nbs�gre tud szab�lyozni!');
%                 end
% 
%                 if debug>1, fprintf('\n   %2d. elem: nyom�sszab�lyz� szelep (%s)',ag_szamlalo,agnev{k}); end
% 
%             case 'legust'
%                 agnev{k}=char(be(eltol+1));                 %...agnev..
%                 csp1{k}  =char(be(eltol+2));
%                 nn(k)   =str2num(char(be(eltol+3)));
%                 VV0(k)  =str2num(char(be(eltol+4)));
%                 pp0(k)  =str2num(char(be(eltol+5)));
%                 AA(k)   =str2num(char(be(eltol+6)));
%                 ll(k)   =str2num(char(be(eltol+7)));
%                 HH(k)   =str2num(char(be(eltol+8)));
%                 ro(k)  =str2num(char(be(eltol+9)));
%                 m0(k)  =str2num(char(be(eltol+10)));
%                 eltol=eltol+11;
%                 if debug>1, fprintf('\n   %2d. elem: legust (%s)',ag_szamlalo,agnev{k}); end
% 
%             case 'szippanto'
%                 agnev{k}=char(be(eltol+1));
%                 csp1{k}=char(be(eltol+2));
%                 csp2{k}='';
%                 pk(k)=str2num(char(be(eltol+3)));
%                 ro(k)=str2num(char(be(eltol+4)));
%                 eltol=eltol+5;
%                 if debug>1, fprintf('\n   %2d. elem: legbeszivo szelep (%s)',ag_szamlalo,agnev{k}); end
%                 
%             case 'buko'
%                 agnev{k}  = char(be(eltol+1));
%                 csp1{k}   = char(be(eltol+2));
%                 csp2{k}   = char(be(eltol+3));
%                 h0(k)     = str2num(char(be(eltol+4)));
%                 Cd(k)     = str2num(char(be(eltol+5)));
%                 B(k)      = str2num(char(be(eltol+6)));
%                 kitevo(k) = str2num(char(be(eltol+7)));
%                 ro(k)     = str2num(char(be(eltol+8)));
%                 m0(k)     = str2num(char(be(eltol+9)));
% 
%                 eltol=eltol+10;
%                 if debug>1, fprintf('\n   %2d. elem: buko (%s)',ag_szamlalo,agnev{k}); end    %              
%                 
%             otherwise
%                 error(['hibas agtipus, cspnev: %8s',cspnev]);
%         end
%     end   %...k
% 
%     %....cpnevsor a csomopontban.......
%     % Itt a nagy ciklus jc-je fel�lír�dik, ami szeintem hib�s. A belsõ
%     % jc-ket lecser�lem jcc-re.
%     cpnevsor{1}=csp1{1};
%     jcc=1;
%     for k=1:nag
%         jel='n';
%         for i=1:jcc
%             if strcmp(cpnevsor{i},csp1{k})
%                 jel='i';
%             end;
%         end;
%         if strcmp(jel,'n')&&~isempty(csp1{k})
%             jcc=jcc+1;
%             cpnevsor{jcc}=csp1{k};
%         end
%         jel='n';
%         for i=1:jcc
%             if strcmp(cpnevsor{i},csp2{k})
%                 jel='i';
%             end
%         end
%         if strcmp(jel,'n')&~isempty(csp2{k})
%             jcc=jcc+1;
%             cpnevsor{jcc}=csp2{k};
%         end
%     end
% 
%     ncp_m=jcc;  %...csomopontszam a merev alrendszerben...
% 
%     %...rangszam kiszamitasa es cpnevsor ellenorzese.............
%     % Itt szint�n lecser�ltem a jc-t jcc-re. HCs
%     rang=zeros(1,ncp_m);
%     for jcc=1:ncp_m
%         mar_volt=0;
%         for k=1:nag
%             if strcmp(cpnevsor{jcc},csp1{k}) |  strcmp(cpnevsor{jcc},csp2{k})
%                 rang(jcc)=rang(jcc)+1;
%                 mar_volt=1;
%             end;
%         end
%         if mar_volt==0
%             error(['\n cp hiba %8s',cpnevsor{jcc}]);
%         end
%     end
% 
%     %rang
%     r_max=max(rang);
%     %r_max
% 
%     %...csom osszeallitas...
%     % meg itt is jc->jcc
%     csom=zeros(ncp_m,r_max);
%     for jcc=1:ncp_m
%         j=0;
%         for k=1:nag
%             if strcmp(cpnevsor{jcc},csp1{k})
%                 j=j+1;
%                 csom(jcc,j) = -k;
%             end
% 
%             if strcmp(cpnevsor{jcc},csp2{k})
%                 j=j+1;
%                 csom(jcc,j) = k;
%             end
%         end
%     end
%     %csom        %...itt a csom vege.......
% 
%     %...csp magassag es elvetel beolvasas....
%     for jcc=1:ncp_m
%         cspnev1{jcc}=char(be(eltol));                    %...csomopontnev
%         jelzo='nincs';
%         for i=1:ncp_m
%             %fprintf('\n %8s  %8s',cspnev1{jcc},cpnevsor{i})
%             if strcmp(cpnevsor{i},cspnev1{jcc})
%                 z(i)=str2double(char(be(eltol+1)));                 %...magassag... [m]
%                 elv(i)=str2double(char(be(eltol+2)))/3600*1000;               %...elvetel... [m3/h] -> [kg/s]
%                 lefutas{i}=char(be(eltol+3));
%                 eltol=eltol+4;
%                 jelzo='van';
%             end
%         end
%         if strcmp(jelzo,'nincs')
%             error('\n\n %8s nevu csomopont nincs\n\n',cspnev1{jcc});
%         end
% 
%     end
% 
%     if debug>2
%         fprintf('\n\n    szam    nev        z[m]    fogy.[kg/s]  lefutas gorbe neve');
%         for jcc=1:ncp_m
%             fprintf('\n     %3i  %8s  %8.2f  %8.2f          %s',jcc,char(cpnevsor{jcc}),z(jcc),elv(jcc),lefutas{jcc});
%         end
%     end
% 
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % ITT NYULTAM BELE!!!!
%     % Merev agelemek letrehozasa.
% 
%     if debug>1, fprintf('\n\n    elemek epitese...'); end
% 
%     clear elemek;
%     for k=1:nag
%         if debug>1, fprintf(' %d ',k); end
%         for i=1:length(cpnevsor)
%             if strcmp(cpnevsor{i},csp1{k}), cspe_szam=i; end
%             if strcmp(cpnevsor{i},csp2{k}), cspv_szam=i; end
%         end
%         switch tipus{k}
% %             case 'konc_cso'
% %                 elemek{k}=konc_cso(agnev{k},cspe_szam,cspv_szam,D(k),L(k),lam(k),ro(k),m0(k));
% 
%             case 'szivattyu'
%                 switch tranziens(k)
%                     case 0     %...allando fordulatszam...
%                         elemek{k}=szivattyu(agnev{k},cspe_szam,cspv_szam,Q,H,Ds(k),Dn(k),ro(k),m0(k),tranziens(k));
%                     case 1     %...kifutas...
%                         elemek{k}=szivattyu(agnev{k},cspe_szam,cspv_szam,Q,H,Ds(k),Dn(k),ro(k),m0(k),tranziens(k),P,n(k),teta(k),tki(k));
%                     case 2     %...inditas, kozelito motor-jelleggorbe...
%                         elemek{k}=szivattyu(agnev{k},cspe_szam,cspv_szam,Q,H,Ds(k),Dn(k),ro(k),m0(k),tranziens(k),P,n(k),teta(k),Mb(k),nb(k),nsz(k));
%                     case 3     %...inditas, adott motor-jelleggorbe...
%                         elemek{k}=szivattyu(agnev{k},cspe_szam,cspv_szam,Q,H,Ds(k),Dn(k),ro(k),m0(k),tranziens(k),P,n(k),teta(k),Mm,nm);
%                     case 4     %...frekivaltos inditas, kozelito motor-jelleggorbe...
%                         elemek{k}=szivattyu(agnev{k},cspe_szam,cspv_szam,Q,H,Ds(k),Dn(k),ro(k),m0(k),tranziens(k),P,n(k),teta(k),Mb(k),nb(k),Mi(k),nsz(k),tbe(k),tind(k));
%                     case 5     %...atemelo szivattyu...
%                         elemek{k}=szivattyu(agnev{k},cspe_szam,cspv_szam,Q,H,Ds(k),Dn(k),ro(k),m0(k),tranziens(k),P,n(k),hbe(k),hki(k),uzem(k));
%                 end
% 
% %             case 'fojtas'
% %                 elemek{k}=fojtas(agnev{k},cspe_szam,cspv_szam,1e-10,K1(k),ro(k),m0(k));
% 
%             case 'vez_fojtas'
%                 elemek{k}=vez_fojtas(agnev{k},cspe_szam,cspv_szam,ro(k),A(k),jgpsz(k),jgpszt(k),m0(k),epsK(k,:),K(k,:),t_vf(k,:),epst(k,:));
% 
%             case 'nyomas'
%                 elemek{k}=nyomas(agnev{k},cspe_szam,pk(k),ro(k));
% 
%             case 'valtozo_nyomas'
%                 elemek{k}=valtozo_nyomas(agnev{k},cspe_szam,ro(k),tt(k,:),pp(k,:));
% 
%             case 'valtozo_tomegaram'
%                 elemek{k}=valtozo_tomegaram(agnev{k},cspe_szam,ro(k),tt(k,:),mm(k,:));
%                 
%             case 'visszacsapo_szelep'
%                 elemek{k}=visszacsapo_szelep(agnev{k},cspe_szam,cspv_szam,ro(k),m0(k));
% 
%             case 'nyomasszabalyzo'
%                 switch szab{k}
%                     case 'p'
%                         for i=1:length(cpnevsor)
%                             if strcmp(cpnevsor{i},szabcsp1{k}), szabcsp_szam=i; end
%                         end
%                         elemek{k}=nyomasszabalyzo(agnev{k},cspe_szam,cspv_szam,ro(k),m0(k),e(k,:),ek(k,:),szab,szabcsp_szam,ajel,pmax,szabP,szabI,szabD,szabcsp,e0,vmax,tbe);
% 
%                     case 'dp'
%                         for i=1:length(cpnevsor)
%                             if strcmp(cpnevsor{i},szabcsp1{k}), szabcsp_szam1=i; end
%                             if strcmp(cpnevsor{i},szabcsp2{k}), szabcsp_szam2=i; end
%                         end
%                         elemek{k}=nyomasszabalyzo(agnev{k},cspe_szam,cspv_szam,ro(k),m0(k),e(k,:),ek(k,:),szab{k},szabcsp_szam1,szabcsp_szam2,ajel(k),pmax(k),szabP(k),szabI(k),szabD(k),szabcsp1{k},szabcsp2{k},e0(k),vmax(k),tbe(k));
%                 end
% 
%             case 'legust'
%                 elemek{k}=legust(agnev{k},cspe_szam,nn(k),VV0(k),pp0(k),AA(k),ll(k),HH(k),ro(k),m0(k));
%             
%             case 'szippanto'
%                 elemek{k}=szippanto(agnev{k},cspe_szam,pk(k),ro(k));
%             
% %             case 'akna'
% %                 elemek{k}=akna(agnev{k},cspe_szam,AA(k),hmin(k),hmax(k),yy0(k),ro(k),m0(k),rajz{k});
%             
%             case 'buko'
%                 elemek{k}=buko(agnev{k},cspe_szam,cspv_szam,h0(k),Cd(k),B(k),kitevo(k),ro(k),m0(k));
%                 
%             case 'nyomovezetek'
%                 elemek{k}=nyomovezetek(agnev{k},cspe_szam,cspv_szam,D(k),L(k),lam(k),ze(k),zv(k),ro(k),m0(k));
%                 
%             otherwise
%                 error(['olvas_merev -> Ismeretlen agelem: ',tipus{k}])
%         end
%     end
%     
%     if debug>1, fprintf('   OK\n    merev alrendszer epitese...'); end;
%     clear csp;
%     maxj=length(csom(1,:));
%     for i=1:ncp_m
%         j=1; clear agak;
%         while j<=maxj && ~csom(i,j)==0
%             %	fprintf('\n csomturka: i=%d  j=%d  csom(%d,%d)=%g',i,j,i,j,csom(i,j));
%             agak(j)=csom(i,j); j=j+1;
%         end
%         csp{i}{1} = i;
%         csp{i}{2} = z(i);
%         csp{i}{3} = agak;
%         csp{i}{4} = 1e5; % nyomas
%         csp{i}{5} = elv(i); % Ez a pillanatnyi fogyasztas, ezt fogjuk szrozni a lefutas pillanatnyi ertekevel.
%         csp{i}{6} = cpnevsor{i};
%         csp{i}{7} = lefutas{i};
%         csp{i}{8} = elv(i); % Ez a nevleges fogy, ehhez nem nyulunk
%     end
% 
%     mar{jc} = merev_alrendszer(1,cspnev,elemek,csp,wdir);
%     
%     if debug>1, fprintf('   OK\n'); end;
%     % EDDIG!!!
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end   %...jc
% 
% %% A gorbek beolvasasa
% 
% no_of_gorbek=0; 
% be=be(eltol:end);
% while ~isempty(be)>0 % Van meg adat
%     if strcmp(be(1),'gorbe')
%         no_of_gorbek= no_of_gorbek+1;
%         nev   = char(be(2));
%         jgpsz = str2num(char(be(3)));
%         for j=1:jgpsz
%             x(j) = str2double(char(be(3+2*j-1)));
%             y(j) = str2double(char(be(3+2*j)));
%         end
%         gorbek{no_of_gorbek}=gorbe(nev,x,y);
%         eltol=4+2*jgpsz; 
%         be=be(eltol:end);        
%         if debug>1, fprintf('\n   %2d. gorbe (%s)\tOK',no_of_gorbek,nev); end
%     else
%         error(['Beolvastam az osszes merev alrendszert, ha van meg adat, csak __gorbe__ lehet. De nem, hanem ',nev,' .'])
%     end
% end
% gorbek{no_of_gorbek+1}=gorbe('const',[0 1e10],[1 1]);
% 
% %% Gorbek hozzadasa az osszes merev alrendszerhez
% for i=1:length(mar)
%     mar{i}.gorbek=gorbek;
% end
