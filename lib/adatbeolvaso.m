function out = adatbeolvaso(fnev,debug,wdir)

if debug>1
    fprintf('\n\n****************************************************');
    fprintf('\n*        ADATFAJL ALRENDSZEREK OLVASASA            *');
    fprintf('\n****************************************************\n\n');
end

if debug>0
    fprintf('\n Feladat neve: %s\n',fnev);
end

[be] = textread(fnev,'%s','delimiter',',');

if debug>0
    fprintf('\n Merev alrendszerek szama : %g',sum(strcmp(be,'mar')));
    fprintf('\n Rugalmas elemek:');
    fprintf('\n   nyomottvizes cso       : %g',sum(strcmp(be,'rugalmas_cso')));
    fprintf('\n   csatorna               : %g',sum(strcmp(be,'csatorna')));
    fprintf('\n   viszkoelasztikus cso   : %g',sum(strcmp(be,'viszkcso')));
    fprintf('\n Gorbek szama             : %g',sum(strcmp(be,'gorbe')));
    fprintf('\n Opciok szama             : %g',sum(strcmp(be,'option')));
end

% Adatstrukturak elokeszitese
s_mar = struct(); % mar
s_option = struct(); % option
s_gorbe = {}; % gorbek

% Rugalmas elemek
data{1}=struct(); % rugalmas cso
data{2}=struct(); % csatorna
data{3}=struct(); % viszkoelasztikus cso
amoba_data=struct();
% rugalmas elemek kozos vektorba lesznek osszefogva, tudnunk kell, az adott
% tipus elofordult-e
type_exists = zeros(1,length(data)+1); % +1 az amoba

while ~isempty(be)
    be=check_comment(be);
    tipus=char(be(1));
    switch tipus
        
        case 'option'
            [be,s_option] = get_option_data(be,s_option,debug);
            
        case 'gorbe'
            [be,s_gorbe] = get_gorbe_data(be,s_gorbe,debug);
            
        case 'mar'
            [be,s_mar] = get_mar_data(be,s_mar,debug);
            
        case 'rugalmas_cso'
            [be,data{1}] = get_rugcso_data(be,data{1},debug);
            type_exists(1)=1;
            
        case 'csatorna'
            [be,data{2}] = get_csatorna_data(be,data{2},debug);
            type_exists(2)=1;
            
        case 'viszkcso'
            [be,data{3}] = get_viszkcso_data(be,data{3},debug);
            type_exists(3)=1;
            
        case 'amoba'
            [be,amoba_data] = get_amoba_data(be,amoba_data,debug);
            type_exists(end)=1;
            
        otherwise
            disp(be);
            error('Ismeretlen elem tipus: %s (option|gorbe|mar|csatorna|rugcso|viszkcso)',tipus);
    end
end

% Merev alrendszerek osszeepitese
% A 'const' gorbet mindig hozza kell adni.
if ~isempty(s_gorbe)
    no = length(s_gorbe)+1;
else
    no=1;
end
s_gorbe{no}=gorbe('const',[0 1e10],[1 1]);

% konstruktrorok
for i=1:length(s_mar.name)
    out.mar{i}=merev_alrendszer(s_mar.name{i},s_mar.elemek{i},s_mar.cspok{i},wdir);
    out.mar{i}.gorbek=s_gorbe;
end

% Option attoltese
out.options=s_option;

% Rugalmas struktura letrehozasa
%-----------------------------------
% rugalmas csompont nevsor epitese
out.rug_struct.cpnevsor=cell(1,1);
for i=1:length(data)
    if type_exists(i)
        out.rug_struct.cpnevsor=build_cpnevsor(out.rug_struct.cpnevsor,data{i},debug);
    end
end

if debug>2
    fprintf('\n\n===================================================\n');
    fprintf('Amoba csomopontok:\n\n  jc    csp. nev');
    for jc=1:length(out.rug_struct.cpnevsor)
        fprintf('\n %3i  %8s',jc,char(out.rug_struct.cpnevsor{jc}))
    end
    fprintf('\n\n');
end

% amoba epitese
if type_exists(end)
    out.rug_struct.amoba = build_amoba(amoba_data,s_gorbe,debug);
    out.rug_struct.amoba_exists=1;
else
    out.rug_struct.amoba_exists=0;
end

% csom epitese
out.rug_struct.csom=build_csom(data,type_exists,out.rug_struct.cpnevsor,debug);        

% konstruktorok hivasa
out.rug_struct.csovek=call_constructors(data,out.rug_struct.csom,type_exists,debug);

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
else
    no=1;
end

% Hozzafuzes
nev   = char(be(2));
jgpsz = str2double(char(be(3)));
x=zeros(1,jgpsz); y=x;
for j=1:jgpsz
    x(j) = str2double(char(be(3+2*j-1)));
    y(j) = str2double(char(be(3+2*j)));
end

% Konstruktor
gorbe_data{no}=gorbe(nev,x,y);

% Egyeb
be=be(4+2*jgpsz:end);

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

if debug>0
    fprintf('\n\n MEREV ALRENDSZER: %s',mar_data.name{mar_no});
end

while ~stop            
    be=check_comment(be);
    tipus=char(be(1));
    switch tipus
        case {'mar','gorbe','option','rugalmas_cso','csatorna','viszkcso'}
            stop=1;
           
        case 'csp'
            data = get_csp_data(be,csp_no,debug);
                        
        case 'konc_cso'
            data = get_konc_cso_data(be,elem_no,debug);
            
        case 'fojtas' 
            data = get_fojtas_data(be,elem_no,debug);
        
        case 'ellenallas' 
            data = get_ellenallas_data(be,elem_no,debug);
        
        case 'akna'
            data = get_akna_data(be,elem_no,debug);
            
        case 'nyomovezetek'
            data = get_nyomovezetek_data(be,elem_no,debug);
 
        case 'szivattyu'
            data = get_szivattyu_data(be,elem_no,debug);
        
        case 'nyomas'
            data = get_nyomas_data(be,elem_no,debug);
        
        case 'visszacsapo_szelep'
            data = get_vcs_data(be,elem_no,debug);
            
        case 'valtozo_nyomas'
            data = get_vny_data(be,elem_no,debug);

        case 'valtozo_tomegaram'
            data = get_vta_data(be,elem_no,debug);            
            
        case 'vez_fojtas'
            data = get_vez_fojtas_data(be,elem_no,debug); 
        
        case 'buko'
            data = get_buko_data(be,elem_no,debug);
            
        case 'legust'
            data = get_legust_data(be,elem_no,debug);
            
        otherwise
            be
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
    fprintf('\n\n    Elem csompontazonositok keresese:');
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

if debug>0
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
function out = get_ellenallas_data(be,elem_no,debug)
agnev      = char(be(2));
out.cspe_nev = char(be(3));
out.cspv_nev = char(be(4));
ro = str2double(char(be(5))); % suruseg [kg/m3]
m0 = str2double(char(be(6))); % kezdeti tomegaram [kg/s]

K1 = str2double(char(be(7))); % ellenallas tenyezo [1/m2]

out.elem = ellenallas(agnev,0,0,1e-10,K1,ro,m0);
out.be = be(8:end);

if debug>1
    fprintf('\n   %2d. elem: fojtas (%s)',elem_no,agnev);
end

%-----------------------------------------------------------------------
function out = get_szivattyu_data(be,elem_no,debug)

agnev        = char(be(2));
out.cspe_nev = char(be(3));
out.cspv_nev = char(be(4));
ro = str2double(char(be(5))); % suruseg [kg/m3]
m0 = str2double(char(be(6))); % kezdeti tomegaram [kg/s]

Ds        = str2double(char(be(7)));  % D szivocsonk [m]
Dn        = str2double(char(be(8)));  % D nyomocsonk [m]
tranziens = str2double(char(be(9)));  % tranziens tipus
jgpsz     = str2double(char(be(10))); % H(Q) jg. pontszam

be=be(11:end);

Q=zeros(1,jgpsz); H=Q; P=Q;
for j=1:jgpsz
    jj=3*j-2;
    Q(1,j)=str2double(char(be(jj)));   % Q [m3/s]
    H(1,j)=str2double(char(be(jj+1))); % H [m]
    P(1,j)=str2double(char(be(jj+2))); % P [kW]
end

be=be(3*jgpsz+1:end);

switch tranziens
    case 0  % allando fordulatszam
        out.elem = szivattyu(agnev,0,0,Q,H,Ds,Dn,ro,m0,tranziens);
        
    case 1     % kifutas
        n    = str2double(char(be(1))); % fordulatszam [1/min]
        teta = str2double(char(be(2))); % tehetetlensegi nyomatek [kgm2]
        tki  = str2double(char(be(3))); % kifutas idopontja [s]
        out.elem = szivattyu(agnev,0,0,Q,H,Ds,Dn,ro,m0,tranziens,P,...
            n,teta,tki);
        
        be=be(4:end);
        
    case 2     % inditas, kozelito motor-jelleggorbe
        n    = str2double(char(be(1))); % fordulatszam [1/min]
        teta = str2double(char(be(2))); % tehetetlensegi nyomatek [kgm2]
        Mb   = str2double(char(be(3))); % billenonyomatek [Nm]
        nb   = str2double(char(be(4))); % billenonyomatek fordulatszama [1/min]
        nsz  = str2double(char(be(5))); % szinkronfordulatszama [1/min]
        out.elem = szivattyu(agnev,0,0,Q,H,Ds,Dn,ro,m0,tranziens,P,...
            n,teta,Mb,nb,nsz);
        
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
        out.elem = szivattyu(agnev,0,0,Q,H,Ds,Dn,ro,m0,tranziens,P,...
            n,teta,Mm,nm);
        
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
        out.elem = szivattyu(agnev,0,0,Q,H,Ds,Dn,ro,m0,tranziens,P,...
            n,teta,Mb,nb,Mi,nsz,tbe,tind);
        
        be=be(9:end);
        
    case 5 % szintkapcsolos atemelo szivattyu
        n    = str2double(char(be(1))); % fordulatszam [1/min]
        hbe  = str2double(char(be(2))); % bekapcsolasi szint [m]
        hki  = str2double(char(be(3))); % kikapcsolasi szint [m]
        uzem = str2double(char(be(4))); % szivattyu uzemel (0/1)
        out.elem = szivattyu(agnev,0,0,Q,H,Ds,Dn,ro,m0,tranziens,P,...
            n,hbe,hki,uzem);
        
        be=be(5:end);
        
    case 6 % frekvenciavaltos, nyomastarto szivattyu
        n     = str2double(char(be(1))); % nevleges fordulatszam [1/min]
        pcsp  = char(be(2));             % szabalyzott csomopont NEVE
        pp    = str2double(char(be(3))); % tartando nyomas (abs) [Pa]
        parP  = str2double(char(be(4))); % szabalyzo aranyos tag
        parI  = str2double(char(be(5))); % szabalyzo integralo tag
        nkezd = str2double(char(be(6))); % indulo fordulatszam
        out.elem = szivattyu(agnev,0,0,Q,H,Ds,Dn,ro,m0,tranziens,P,...
            n,pcsp,pp,parP,parI,nkezd);
                    
        be=be(7:end);
        
    otherwise
        error(['Szivattyu ',agnev,': hibas tranziens jel: ',tranziens]);
end

if debug>0
    fprintf('\n   %2d. elem: szivattyu (%s)',elem_no,agnev);
end

out.be=be;
 
%-----------------------------------------------------------------------
function out = get_nyomas_data(be,elem_no,debug)

agnev        = char(be(2));
out.cspe_nev = char(be(3));
out.cspv_nev = '';
ro = str2double(char(be(4))); % suruseg [kg/m3]
m0 = str2double(char(be(5))); % kezdeti tomegaram [kg/s]

pk = str2double(char(be(6))); % adott nyomas [Pa]

out.elem = nyomas(agnev,0,pk,ro);

if debug>0
    fprintf('\n   %2d. elem: allando nyomas (%s)',elem_no,agnev);
end

out.be=be(7:end);

%-----------------------------------------------------------------------
function out = get_vcs_data(be,elem_no,debug)

agnev        = char(be(2));
out.cspe_nev = char(be(3));
out.cspv_nev = char(be(4));
ro = str2double(char(be(5))); % suruseg [kg/m3]
m0 = str2double(char(be(6))); % kezdeti tomegaram [kg/s]
 
out.elem = visszacsapo_szelep(agnev,0,0,ro,m0);

if debug>0
    fprintf('\n   %2d. elem: visszacsapo_szelep (%s)',elem_no,agnev);
end

out.be = be(7:end);

%-----------------------------------------------------------------------
function out = get_vta_data(be,elem_no,debug)

agnev        = char(be(2));
out.cspe_nev = char(be(3));
out.cspv_nev = '';
ro = str2double(char(be(4)));    % suruseg [kg/m3]
m0 = str2double(char(be(5)));    % kezdeti tomegaram [kg/s]
if strcmp(char(be(6)),'file')
    jgpsz = 0.5;
    f01 = char(be(7));
    data = xlsread(f01);
    tt = data(:,1);
    mm = data(:,2);
else    
    jgpsz = str2double(char(be(6))); % m(t) jellegg. pontsz.    
    for j=1:jgpsz
        jj=2*j-1;
        tt(j) = str2double(char(be(6+jj)));
        mm(j) = str2double(char(be(6+jj+1)));
    end
end
out.elem = valtozo_tomegaram(agnev,0,0,ro,m0,tt,mm);

if debug>0
    fprintf('\n   %2d. elem: valtozo_tomegaram (%s)',elem_no,agnev);
end

out.be = be(7+2*jgpsz:end);

%-----------------------------------------------------------------------
function out = get_vny_data(be,elem_no,debug)

agnev        = char(be(2));
out.cspe_nev = char(be(3));
out.cspv_nev = '';
ro = str2double(char(be(4)));    % suruseg [kg/m3]
m0 = str2double(char(be(5)));    % kezdeti tomegaram [kg/s]

if strcmp(char(be(6)),'file')
    jgpsz = 0.5;
    f01 = char(be(7));
    data = xlsread(f01);
    tt = data(:,1);
    pp = data(:,2);
else
    jgpsz = str2double(char(be(6))); % p(t) jellegg. pontsz.
    for j=1:jgpsz
        jj=2*j-1;
        tt(j) = str2double(char(be(6+jj)));
        pp(j) = str2double(char(be(6+jj+1)));
    end
end
out.elem = valtozo_nyomas(agnev,0,0,ro,m0,tt,pp);

if debug>0
    fprintf('\n   %2d. elem: valtozo_nyomas (%s)',elem_no,agnev);
end

out.be = be(7+2*jgpsz:end);

%-----------------------------------------------------------------------
function out = get_vez_fojtas_data(be,elem_no,debug)

agnev        = char(be(2));
out.cspe_nev = char(be(3));
out.cspv_nev = char(be(4));
ro = str2double(char(be(5))); % suruseg [kg/m3]
m0 = str2double(char(be(6))); % kezdeti tomegaram [kg/s]

A      = str2double(char(be(7))); % nevleges keresztmetszet [m2]
jgpsz  = str2double(char(be(8))); % K(eps) jellegg. pontsz.
jgpszt = str2double(char(be(9))); % eps(t) fv.
be=be(10:end);

epsK=zeros(jgpsz,1); K=epsK;
for j=1:jgpsz
    jj=2*j-1;
    epsK(j) = str2double(char(be(jj)));   % rel. helyzet [-]
    K(j)    = str2double(char(be(jj+1))); % ellenallas tenyezo [1/m2]    
end
be=be(2*jgpsz+1:end);

t_vf=zeros(jgpszt,1); epst=t_vf;
for j=1:jgpszt
    jj=2*j-1;
    t_vf(j)=str2double(char(be(jj)));   % ido [s]
    epst(j)=str2double(char(be(jj+1))); % rel. helyzet [-]
end
be=be(2*jgpszt+1:end);

out.elem = vez_fojtas(agnev,0,0,ro,A,jgpsz,jgpszt,m0,epsK,K,t_vf,epst);

if debug>0
    fprintf('\n   %2d. elem: vezerelt fojtas (%s)',elem_no,agnev);
end

out.be = be;

%-----------------------------------------------------------------------
function out = get_buko_data(be,elem_no,debug)

agnev        = char(be(2));
out.cspe_nev = char(be(3));
out.cspv_nev = char(be(4));
ro = str2double(char(be(5))); % suruseg [kg/m3]
m0 = str2double(char(be(6))); % kezdeti tomegaram [kg/s]

h0     = str2double(char(be(7)));
Cd     = str2double(char(be(8)));
B      = str2double(char(be(9)));
kitevo = str2double(char(be(10)));

be=be(11:end);

out.elem = buko(agnev,0,0,h0,Cd,B,kitevo,ro,m0);

if debug>0
    fprintf('\n   %2d. elem: buko (%s)',elem_no,agnev);
end

out.be = be;

function out = get_legust_data(be,elem_no,debug)

agnev        = char(be(2));
out.cspe_nev = char(be(3));
out.cspv_nev = '';
ro = str2double(char(be(4))); % suruseg [kg/m3]
m0 = str2double(char(be(5))); % kezdeti tomegaram [kg/s]

nn = str2double(char(be(6)));
VV0 = str2double(char(be(7)));
pp0 = str2double(char(be(8)));
AA = str2double(char(be(9)));
ll = str2double(char(be(10)));
HH = str2double(char(be(11)));

% h0     = str2double(char(be(7)));
% Cd     = str2double(char(be(8)));
% B      = str2double(char(be(9)));
% kitevo = str2double(char(be(10)));

%                 agnev{k}=char(be(eltol+1));                 %...agnev..
%                 csp1{k}  =char(be(eltol+2));///
%                 nn(k)   =str2num(char(be(eltol+3)));
%                 VV0(k)  =str2num(char(be(eltol+4)));
%                 pp0(k)  =str2num(char(be(eltol+5)));
%                 AA(k)   =str2num(char(be(eltol+6)));
%                 ll(k)   =str2num(char(be(eltol+7)));
%                 HH(k)   =str2num(char(be(eltol+8)));
%                 ro(k)  =str2num(char(be(eltol+9)));///
%                 m0(k)  =str2num(char(be(eltol+10)));///
%                 eltol=eltol+11;
%                 if debug>1, fprintf('\n   %2d. elem: legust
%                 (%s)',ag_szamlalo,agnev{k}); end

be=be(12:end);

out.elem = legust(agnev,0,0,ro,m0,nn,VV0,pp0,AA,ll,HH);

if debug>0
    fprintf('\n   %2d. elem: legust (%s)',elem_no,agnev);
end

out.be = be;

%-----------------------------------------------------------------------
function be = check_comment(be)

tmp=char(be(1));
if length(tmp)>1
    if strcmp(tmp(1:2),'/*')
        be=be(2:end);
        tmp=char(be(1));
    end
end

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

%             
%             case 'valtozo_nyomas'
%                 elemek{k}=valtozo_nyomas(agnev{k},cspe_szam,ro(k),tt(k,:),pp(k,:));
% 
%             case 'valtozo_tomegaram'
%                 elemek{k}=valtozo_tomegaram(agnev{k},cspe_szam,ro(k),tt(k,:),mm(k,:));
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

%------------------------------------------------------------------------
function [be,rugcso_data] = get_rugcso_data(be,rugcso_data,debug)

% Az elem szama
if isfield(rugcso_data,'nev')
    no = length(rugcso_data.nev)+1;
else
    no=1;
end

if (length(be)<14)
    error('HIBA: rugalmas cso (neve:%s), nincs eleg adat: %g van de min. %g kellene!',...
        be{2},length(be),14);
end

agnev=char(be(2));            % agnev
cpele=char(be(3));            % cpnev eleje
cpveg=char(be(4));            % cpnev vege
ro   =str2double(char(be(5)));  % suruseg [kg/m3]
m    =str2double(char(be(6)));  % tomegaram, t=0 [m3/s]
pe   =str2double(char(be(7)));  % nyomas, t=0, cso elejen, [Pa]

d    =str2double(char(be(8)));   % atmero [m]
lam  =str2double(char(be(9)));   % lambda
del  =str2double(char(be(10)));   % falvastagsag [m]
L    =str2double(char(be(11)));   % hossz [m]
Ec   =str2double(char(be(12)));   % cso rug.mod [Pa]
Ef   =str2double(char(be(13)));  % foly. rug.mod [Pa]
opsz =str2double(char(be(14)));  % osztaspontok szama
op_h =char(be(15));  % auto vagy user

Er =(1/Ef+d/del/Ec)^(-1);
aa =sqrt(Er/ro);
dx =L/(opsz-1);
dt =dx/aa;

if debug>1
    fprintf('\n\n===================================================\n tipus: rugalmas cso (%d.)\n',no);
    fprintf('\n         agnev     cp_1  ->  cp_2    d[m]     del[m]   L[m]     E[Pa]    opsz  m0[kg/s]  p0[bar]  ro[kg/m^3]   Ef[Pa]   a[m/s]   dt[s]');
    fprintf('\n--------------------------------------------------------------------------------------------------------------------------------------------');
    fprintf('\n %2d.  %8s  %8s  %8s  %6.4f   %6.4f   %6.1f   %6.2e  %3d   %6.2f    %5.2f      %4.1f    %6.2e  %7.2f  %3.3f',no,agnev,cpele,cpveg,d,del,L,Ec,opsz,m,pe/1e5,ro,Ef,aa,dt);
end

% osztasponti adatok
Ds=L/(opsz-1);
A=d^2*pi/4;
seb=m/ro/A;
Dpv=lam*Ds/d*ro/2*seb*abs(seb);

switch op_h
    case 'auto'
        ze=str2double(char(be(16)));
        zv=str2double(char(be(17)));
        h=interp1([0 L],[ze zv],linspace(0,L,opsz),'linear');
        be=be(17+1:end);
    case 'user'
        for i=1:opsz
            h(i)=str2double(char(be(15+i)));
        end
        be=be(15+opsz+1:end);
    otherwise
        error('Ismeretlen a %d. rugalmas cso osztaspont magassaganak feltoltese: %s, megengedett: user,auto',no,op_h);
end


for i=1:opsz    
    if i==1
        p(i)=pe;
    else
        p(i)=p(i-1)-ro*9.81*(h(i)-h(i-1))-Dpv; 
    end
end

if debug>2
    fprintf('\n\n   opsz.    h          p        v');
    fprintf('\n  -------------------------------------');
    for i=1:opsz
        fprintf('\n   %3i  %7.3f  %10.0f  %6.3f',i,h(i),p(i),seb);
    end
    fprintf('\n  -------------------------------------');
end

% Adatok beolvasasa kesz, utolagos muveletek

rugcso_data.nev{no}   = agnev;
rugcso_data.cpele{no} = cpele;           
rugcso_data.cpveg{no} = cpveg;
rugcso_data.d(no)     = d;
rugcso_data.lam(no)   = lam;
rugcso_data.del(no)   = del;
rugcso_data.L(no)     = L;
rugcso_data.Ec(no)    = Ec;
rugcso_data.opsz(no)  = opsz;
rugcso_data.ro(no)    = ro;
rugcso_data.Ef(no)    = Ef;
rugcso_data.seb(no)   = seb;
rugcso_data.p{no}     = p;
rugcso_data.h{no}     = h;

% Ilyen a rugcso konstruktora:
% agnev,cpele_szam,cpveg_szam,d,lam,del,L,Ec,opsz,ro,Ef,seb,p,h

if debug>2
    fprintf('\n\n %s nevu rugalmas cso adatainak beolvasasa kesz.',agnev);
end

%------------------------------------------------------------------------
function [be,data] = get_csatorna_data(be,data,debug)

% Az elem szama
if isfield(data,'nev')
    no = length(data.nev)+1;
else
    no=1;
end

if (length(be)<15)
    error('HIBA: csatorna (neve:%s), nincs eleg adat: %g van de min. %g kellene!',...
        be{2},length(be),14);
end

agnev = char(be(2));               % agnev
cpele = char(be(3));               % cpnev eleje
cpveg = char(be(4));               % cpnev vege
ro    = str2double(char(be(5)));  % suruseg [kg/m3]
m     = str2double(char(be(6)));  % tomegaram, t=0 [m3/s]
tipus = char(be(7));               % tipus: kor,teglalap
switch tipus
    case 'kor'
        dvB = str2double(char(be(8)));   % atmero [m]
    case 'teglalap'
        dvB = str2double(char(be(8)));   % szelesseg [m]
    otherwise
        error('Ismeretlen a %d. csatorna tipusa: %s, megengedett: kor,teglalap',no,tipus);
end
L     = str2double(char(be(9)));   % hossz [m]
ze    = str2double(char(be(10)));   % eleje fenekszint [m]
zv    = str2double(char(be(11)));   % vege fenekszint [m]
n     = str2double(char(be(12)));  % Manning-fele n
y0    = str2double(char(be(13)));  % vizszint t=0, [m]
opsz  = str2double(char(be(14)));  % osztaspontok szama
op_h  = char(be(15));              % osztaspontok magassaganak feltoltese: 'auto', 'user'
rajz  = char(be(16));              % rajzoljon? 'rajz' vagy barmi mas
dttype= char(be(17));              % 'auto' vagy egy szam [s]

v = zeros(1,opsz);
y = ones(1,opsz)*y0;

switch op_h
    case 'auto'
        h=interp1([0 L],[ze zv],linspace(0,L,opsz),'linear');
        opsz_eltol = 0;
    case 'user'
        for i=1:opsz
            h(i)=str2double(char(be(15+i)));
        end
    otherwise
        error('Ismeretlen a %d. csatorna osztaspont magassaganak feltoltese: %s, megengedett: user,auto',no,op_h);
end

if debug>1
    if debug>2
        fprintf('\n\n===================================================\n tipus: csatorna (%d.)\n',no);
    end
    fprintf('\n         agnev     cp_1  ->  cp_2   dvB[m]    n[-]     L[m]     ze[m]    zv[m]  m0[kg/s]   y0[m]  ro[kg/m^3]  opsz');
    fprintf('\n----------------------------------------------------------------------------------------------------------------------');
    fprintf('\n %2d.  %8s  %8s  %8s  %6.4f   %6.4f   %6.1f   %6.2f  %6.2f    %5.2f      %4.1f    %6.2e  %d\n',no,agnev,cpele,cpveg,dvB,n,L,ze,zv,m,y0,ro,opsz);
end

if debug>2
    fprintf('\n   opsz.    z          y        v');
    fprintf('\n  -------------------------------------');
    for i=1:opsz
        fprintf('\n   %3i  %7.3f  %7.3f  %6.3f',i,h(i),y(i),v(i));
    end
    fprintf('\n  -------------------------------------');
end

% Adatok beolvasasa kesz, utolagos muveletek

be=be(17+opsz_eltol+1:end);

data.nev{no}   = agnev;
data.cpele{no} = cpele;           
data.cpveg{no} = cpveg;
data.tipus{no} = tipus;
data.dvB(no)   = dvB;
data.n(no)     = n;
data.L(no)     = L;
data.ze(no)    = ze;
data.zv(no)    = zv;
data.y0(no)    = y0;
data.ro(no)    = ro;
data.y{no}     = y;
data.v{no}     = v;
%data.h{no}     = h;
data.rajz{no}  = rajz;
data.dttype{no}=dttype;

% Ilyen a csatorna konstruktora:

if debug>2
    fprintf('\n\n %s nevu csatorna adatainak beolvasasa kesz.',agnev);
end

%------------------------------------------------------------------------
function [be,viszkcso_data] = get_viszkcso_data(be,viszkcso_data,debug)

% Az elem szama
if isfield(viszkcso_data,'nev')
    no = length(viszkcso_data.nev)+1;
else
    no=1;
end

if (length(be)<16)
    error('HIBA: viszkoelasztikus cso (neve:%s), nincs eleg adat: %g van de min. %g kellene!',...
        be{2},length(be),16);
end

%viszk_cso,vcso1,csp1,csp2,0.006,9e-6,0.001,1.86,2e2,30,0.00,1.1e5,900,2e2,9e6,9e5,150000

agnev=char(be(2));            % agnev
cpele=char(be(3));            % cpnev eleje
cpveg=char(be(4));            % cpnev vege
d0    =str2double(char(be(5)));   % atmero [m]
nu  =str2double(char(be(6)));    % kin. viszkozitas [m^2/s]
del  =str2double(char(be(7)));   % falvastagsag [m]
L    =str2double(char(be(8)));   % hossz [m]
%Ec   =str2double(char(be(9)));   % cso rug.mod [Pa]
%Ef   =str2double(char(be(10)));  % foly. rug.mod [Pa]
m    =str2double(char(be(9)));   % tomegaram, t=0 [m3/s]
pe   =str2double(char(be(10)));  % nyomas, t=0, cso elejen, [Pa]
ro   =str2double(char(be(11)));  % suruseg [kg/m3]
opsz =str2double(char(be(12)));  % osztaspontok szama
op_h =char(be(13));  % osztaspontok magassaganak feltoltese: 'auto', 'user'
E1   =str2double(char(be(14)));  % E1 rug. mod. (Stuart modell)
E2   =str2double(char(be(15)));  % E2 rug. mod. (Stuart modell)
eta2 =str2double(char(be(16)));  % eta2 csillap�t�s (Stuart modell)

%Er =(1/Ef+d/del/Ec)^(-1);
%aa =sqrt(Er/ro);

% Kezdeti hull�msebess�g be�ll�t�sa

alp = d0/del;
gam = (1 + alp*pe/E1);   % gamma --> A
%bet = sqrt(E1/(alp*ro*gam));     % beta --> lambda     
a0 = sqrt(E1*gam/(alp*ro));
dx =L/(opsz-1);
dt0 =dx/a0;

if debug>1
    %fprintf('\n\n===================================================\n tipus: viszkoelasztikus cso (%d.)\n',no);
    %fprintf('\n         agnev     cp_1  ->  cp_2    d[m]     del[m]   L[m]     E[Pa]    opsz  m0[kg/s]  p0[bar]  ro[kg/m^3]   Ef[Pa]   a[m/s]   dt[s]');
    %fprintf('\n--------------------------------------------------------------------------------------------------------------------------------------------');
    %fprintf('\n %2d.  %8s  %8s  %8s  %6.4f   %6.4f   %6.1f   %6.2e  %3d   %6.2f    %5.2f      %4.1f    %6.2e  %7.2f  %3.3f',no,agnev,cpele,cpveg,d0,del,L,Ec,opsz,m,pe/1e5,ro,Ef,aa,dt);
end

% osztasponti adatok
Ds = L/(opsz-1);
A = d0^2*pi/4;
seb = m/ro/A;
%Dpv = lam*Ds/d*ro/2*seb*abs(seb);
Dpv = 32*nu/(d0^2)*abs(seb)*ro*Ds; % Nyom�svesztes�g

switch op_h
    case 'auto'
        ze=str2double(char(be(17)));
        zv=str2double(char(be(18)));
        h=interp1([0 L],[ze zv],linspace(0,L,opsz),'linear');
        %opsz = 0;
        opsz_dummy = 2;
    case 'user'
        for i=1:opsz
            h(i)=str2double(char(be(16+i)));
        end
        opsz_dummy = opsz;
    otherwise
        error('Ismeretlen a %d. viszkcso osztaspont magassaganak feltoltese: %s, megengedett: user,auto',no,op_h);
end

for i=1:opsz
    %h(i) = str2double(char(be(15+i)));
    a(i) = a0;
    d(i) = d0;
    if i==1
        p(i)=pe;
    else
        p(i)=p(i-1)-ro*9.81*(h(i)-h(i-1))-Dpv; 
    end
end

if debug>2
    fprintf('\n\n   opsz.    h          p        v         a');
    fprintf('\n  -------------------------------------');
    for i=1:opsz
        fprintf('\n   %3i  %7.3f  %10.0f  %6.3f %6.3f',i,h(i),p(i),seb,a(i));
    end
    fprintf('\n  -------------------------------------');
end

% Adatok beolvasasa kesz, utolagos muveletek

be = be(16+opsz_dummy+1:end);

viszkcso_data.nev{no}   = agnev;
viszkcso_data.cpele{no} = cpele;           
viszkcso_data.cpveg{no} = cpveg;
viszkcso_data.d0(no)     = d0;
viszkcso_data.nu(no)   = nu;
viszkcso_data.del(no)   = del;
viszkcso_data.L(no)     = L;
%viszkcso_data.Ec(no)    = Ec;
viszkcso_data.opsz(no)  = opsz;
viszkcso_data.ro(no)    = ro;
%viszkcso_data.Ef(no)    = Ef;
viszkcso_data.seb(no)   = seb;
viszkcso_data.p{no}     = p;
viszkcso_data.h{no}     = h;
viszkcso_data.E1(no)     = E1;
viszkcso_data.E2(no)     = E2;
viszkcso_data.eta2(no)     = eta2;

% Ilyen a rugcso konstruktora:
% agnev,cpele_szam,cpveg_szam,d,lam,del,L,Ec,opsz,ro,Ef,seb,p,h

if debug>2
    fprintf('\n\n %s nevu viszkoelasztikus cso adatainak beolvasasa kesz.',agnev);
end

%------------------------------------------------------------------------
function [be,amoba_data] = get_amoba_data(be,amoba_data,debug)

% Az elem szama
if isfield(amoba_data,'nev')
    no = length(amoba_data.nev)+1;
else
    no=1;
end

if (length(be)<4)
    error('HIBA: amoba csp. (neve:%s), nincs eleg adat: %g van de min. %g kellene!',...
        be{2},length(be),4);
end

amoba_data.nev{no}  = char(be(2));
amoba_data.h(no)    = str2double(char(be(3)));
amoba_data.fogy(no) = str2double(char(be(4)))/3600*1000;
amoba_data.gorbe{no} = char(be(5));
amoba_data.nevl(no) = str2double(char(be(4)))/3600*1000;

be=be(6:end);

if debug>1
    fprintf('\n\n===================================================\n tipus: amoba csomopont (%d.)\n',no);
    fprintf('\n nev : %s',amoba_data.nev{no});   
    fprintf('\n h   : %g [m]',amoba_data.h(no));
    fprintf('\n fogy: %g [kg/s]',amoba_data.fogy(no));
    fprintf('\n gorbe: %s',amoba_data.gorbe{no});
end

if debug>2
    fprintf('\n\n %s nevu amoba csompont adatainak beolvasasa kesz.',amoba_data.nev{no});
end

%------------------------------------------------------------------------
function cpnevsor = build_cpnevsor(cpnevsor,data,debug)

offset=length(cpnevsor);
if offset==1
    offset=0;
end

for i=1:length(data.nev)
    cpnevsor{offset+i}=data.cpele{i};
end
for i=1:length(data.nev)
    cpnevsor{offset+length(data.nev)+i}=data.cpveg{i};
end
cpnevsor=unique(cpnevsor);

%------------------------------------------------------------------------
function csom = build_csom(data,type_exists,cpnevsor,debug)
% csom matrix:
% sorok: csomopont szam
% oszlopol: erkezo (-) ill. tavozo (+) csovek szama elojelesen

for i=1:length(cpnevsor)
    j=0;
    for type_no=1:length(type_exists)-1 % Az utolso type az amoba csp.
        if type_exists(type_no)
            for k=1:length(data{type_no}.nev)
                if strcmp(cpnevsor{i},data{type_no}.cpele{k})
                    j=j+1;
                    csom(i,j)=-k;
                end
                if strcmp(cpnevsor{i},data{type_no}.cpveg{k})
                    j=j+1;
                    csom(i,j)=k;
                end
            end
        end
    end
end

%------------------------------------------------------------------------
function csovek = call_constructors(data,csom,type_exists,debug)

% Kikeressuk az eleje es vege csompont szamat es meghivjuk a konstruktort
db=1;
for type=1:length(data)
    if type_exists(type)
        for k=1:length(data{type}.nev)
            data{type}.cpele_szam(k)=0;
            data{type}.cpveg_szam(k)=0;
            for i=1:length(csom(:,1))
                for j=1:length(csom(1,:))
                    if csom(i,j)==-k
                        data{type}.cpele_szam(k)=i;
                    end
                    if csom(i,j)==k
                        data{type}.cpveg_szam(k)=i;
                    end
                end
            end

            if data{type}.cpele_szam(k)==0
                error('Nem talalom az eleje csomopontot: %s elem',data{type}.nev{k});
            elseif data{type}.cpveg_szam(k)==0
                error('Nem talalom a vege csomopontot: %s elem',data{type}.nev{k});
            else
                % Rugalmas cso
                if type==1
                    if debug>1
                        fprintf('\n%2d./%2d rugalmas cso letrehozasa... ',...
                            k,length(data{type}.nev));
                    end
                    tmp=data{type};
                    csovek{db}=cso(...
                        tmp.nev{k},...
                        tmp.cpele_szam(k),...
                        tmp.cpveg_szam(k),...
                        tmp.d(k),...
                        tmp.lam(k),...
                        tmp.del(k),...
                        tmp.L(k),...
                        tmp.Ec(k),...
                        tmp.opsz(k),...
                        tmp.ro(k),...
                        tmp.Ef(k),...
                        tmp.seb(k),...
                        tmp.p{k},...
                        tmp.h{k});
                    if debug>1, fprintf(' OK'); end
                end
                
                % Csatorna
                
                if type==2
                    if debug>1
                        fprintf('\n%2d./%2d csatorna letrehozasa... ',...
                            k,length(data{type}.nev));
                    end
                    tmp=data{type};
                    csovek{db}=csatorna(...
                        tmp.nev{k},...
                        tmp.cpele_szam(k),...
                        tmp.cpveg_szam(k),...
                        tmp.tipus{k},...
                        tmp.dvB(k),...
                        tmp.L(k),...
                        tmp.ze(k),...
                        tmp.zv(k),...
                        tmp.n(k),...
                        tmp.ro(k),...
                        tmp.y{k},...
                        tmp.v{k},...
                        tmp.rajz{k},...
                        tmp.dttype{k});
                    if debug>1, fprintf(' OK'); end
                end
                
                % Viszkoelasztikus cso
                
                if type == 3
                    if debug>1
                        fprintf('\n%2d./%2d viszkoelasztikus cso letrehozasa... ',...
                            k,length(data{type}.nev));
                    end
                    tmp = data{type};
                    csovek{db}=viszkcso(...
                        tmp.nev{k},...
                        tmp.cpele_szam(k),...
                        tmp.cpveg_szam(k),...
                        tmp.d0(k),...
                        tmp.nu(k),...
                        tmp.del(k),...
                        tmp.L(k),...
                        tmp.opsz(k),...
                        tmp.ro(k),...
                        tmp.seb(k),...
                        tmp.p{k},...
                        tmp.h{k},...
                        tmp.E1(k),...
                        tmp.E2(k),...
                        tmp.eta2(k));

                    if debug>1, fprintf(' OK'); end                    
                end
                
            end
            db=db+1;
        end
    end
end

%------------------------------------------------------------------------
function amoba = build_amoba(amoba_data,s_gorbe,debug)

amoba{1}{6} = s_gorbe;

for i=1:length(amoba_data.nev)
    amoba{i}{1}=amoba_data.nev{i}; 
    amoba{i}{2}=amoba_data.fogy(i);
    amoba{i}{3}=amoba_data.h(i);
    amoba{i}{4}=amoba_data.gorbe{i};
    amoba{i}{5}=amoba_data.nevl(i);
end
