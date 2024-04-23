function trsz(df,tmax,varargin)

global dtmin

close all
fclose all;
fprintf('\nTranziens szimulator');
fprintf('\n------------------------------\n\n');

%% Elokeszuletek
lepesmax=10; cput=cputime; t=0; lepes=0; tt=0;
debug=0; dtmin_set='auto';
options=varargin;

while length(options)>1
    mit=options{1};
    mire=options{2};
    options=options(3:end);
    switch mit
        case 'debug', debug=mire;
        case 'dtmin'
            dtmin=mire;
            dtmin_set='user';
        case 't0', t=mire;
        otherwise,      fprintf('\n\n HIBA: Ismeretlen opcio. Lehetseges ertekek: debug(0|1|2|3), rajz(0|1), dtmin(xx s)\n');
    end
end

dtout = (tmax-t)/lepesmax;

%% Esetleges maradvanyfajlok torlese. 
fprintf('Konyvtar tisztatasa');

resfiles=dir('*.res');
if ~isempty(resfiles)
    for i=1:length(resfiles)
        if debug>1
            fprintf('\n\t%s torlese...',resfiles(i).name);
        end
        delete(resfiles(i).name);
    end
end

outfiles=dir('*.out');
if ~isempty(outfiles)
    for i=1:length(outfiles)
        if debug>1
            fprintf('\n\t%s torlese... ',outfiles(i).name);
        end
        delete(outfiles(i).name);
    end
end

fprintf('\t\t\t\tOK');

%% Esetleges abrak torlese
figs_open=get(0,'Children');
for i=1:length(figs_open)
    figure(figs_open(i));
    clf;
end

%% Working directory megallapitasa
wdir = what;

%% Adatbeolvasas
fprintf('\nAdatfajl olvasasa');
out=adatbeolvaso(df,debug,wdir.path);
mar        = out.mar;
options    = out.options;
rug_struct = out.rug_struct;
csovek     = rug_struct.csovek;
csp_nevsor = rug_struct.cpnevsor;
csom       = rug_struct.csom;
if rug_struct.amoba_exists
    amoba  = rug_struct.amoba;
end

%% Rendszer epitese
%-----------------------------------------------
%% Merev alrendszer

if debug>2
    fprintf('\n\n\n***********************************************');
    fprintf('\n*     Merev alrendszer csomopont Biblia       *');
    fprintf('\n***********************************************\n');
    fprintf('\n  mar.nev    mar.ssz.    csp.nev     csp.szam     h[m]     fogy.[kg/s]');
    fprintf('\n---------------------------------------------------------------------------');
end

k=1;
for i=1:length(mar)
    for j=1:length(mar{i}.csp)
        
        csp=mar{i}.csp;
       
        if debug>2
            fprintf('\n%9s      %2d      %9s       %2d      %7.2f   %7.2f',mar{i}.nev,i,csp{j}{6},j,csp{j}{2},csp{j}{5});
        end
        cspb{k}{1}=i;
        cspb{k}{2}=j;
        cspb{k}{3}=csp{j}{6};
        k=k+1;
    end
    if debug>2
        fprintf('\n---------------------------------------------------------------------------');
    end
end
if debug>2
    pause;
end
if debug<1, fprintf('\t\t\t\tOK'); end

% Rugalmas csovek
%-------------------------------------
if debug<1, fprintf('\nRugalmas csovek epitese\t'); end

% dtmin beallatasa
%-------------------------------------
% Ha ennal kisebbet lepne ket cso, egyszerre
% leptetjuk oket.
switch dtmin_set
    case 'auto'
        dtmin=1e5;
        for i=1:length(csovek), dtmin=min(dtmin,csovek{i}.dt); end
        dtmin=dtmin/1000;
end

% Tejcsarda (tjcs) matrix osszeallitasa, azaz a
% rugalmas as a merev alrendszerek csatlakozasa, az utolso
% sor amoba csomopont eseten a fogyasztas
%----------------------------------------------------
n_cso=length(csovek); n_mar=length(mar);
tjcs=zeros(length(csp_nevsor),n_cso+n_mar+1);

if debug>1
    fprintf('\n\n****************************************************');
    fprintf('\n*   RUGALMAS ES MEREV ALRENDSZEREK CSATLAKOZaSA    *');
    fprintf('\n****************************************************\n');
end
if debug<1, fprintf('\nRugalmas csovek es merev alrendszerek csatlakozasa\t'); end

% vegigmegyunk a csoveken....
tjcs_sor = [];
tjcs_oszlop = [];
amoba_no = [];
for i=1:length(csovek)
    temp=csovek{i}.csp;
    cspe_szam=temp(1); cspv_szam=temp(2);
    cspe=csp_nevsor{cspe_szam}; cspv=csp_nevsor{cspv_szam};
    
    % vegigmegyunk a csomopont Biblian,
    % es megkeressuk az elejen es vegen levo
    % csomapontot a neve alapjan.
    emegvan=0; vmegvan=0;
    for j=1:length(cspb)
        if strcmp(cspe,cspb{j}{3})
            if debug>2
                fprintf('\ncso:%5s  %8s -> %8s | %8s | ',csovek{i}.nev,cspe,cspv,cspb{j}{3});
                fprintf(' eleje megvan: %d.mar %d.csp',cspb{j}{1},cspb{j}{2});
            end
            k=1;
            while k<=length(csp_nevsor) && ~strcmp(cspe,csp_nevsor{k}), k=k+1; end
            tjcs(k,i)=-1; tjcs(k,n_cso+cspb{j}{1})=cspb{j}{2};
            gcsp_tipus(k) =1;
            emegvan=1;
        end
        
        if strcmp(cspv,cspb{j}{3})
            if debug>2
                fprintf('\ncso:%5s  %8s -> %8s | %8s | ',csovek{i}.nev,cspe,cspv,cspb{j}{3});
                fprintf(' vege  megvan: %d.mar %d.csp',cspb{j}{1},cspb{j}{2});
            end
            k=1;
            while k<=length(csp_nevsor) && ~strcmp(cspv,csp_nevsor{k}), k=k+1; end
            tjcs(k,i)=1; tjcs(k,n_cso+cspb{j}{1})=cspb{j}{2};
            gcsp_tipus(k) =1;
            vmegvan=1;
        end
    end
    
    % Ha nincs meg a csomopont, amobarol van szo...
    if emegvan==0
        if debug>2
            fprintf('\ncso:%5s  %8s -> %8s | %8s | ',csovek{i}.nev,cspe,cspv,cspb{j}{3});
            fprintf(' eleje nincs meg -> amoba');
        end
        k=1;
        while k<=length(csp_nevsor) && ~strcmp(cspe,csp_nevsor{k}), k=k+1; end
        tjcs(k,i)=-1; gcsp_tipus(k) =2;
        iii=1;
        % jav�tva iii<=length(amoba) -r�l; BG; 2010.01.29.
        while (iii<length(amoba) && ~strcmp(cspe,amoba{iii}{1})), iii=iii+1; end
        tjcs(k,n_cso+n_mar+1)=amoba{iii}{2};
        tjcs_sor = [tjcs_sor k];
        tjcs_oszlop = [tjcs_oszlop n_cso+n_mar+1];
        amoba_no = [amoba_no iii];
    end
    if vmegvan==0
        if debug>2
            fprintf('\ncso:%5s  %8s -> %8s | %8s | ',csovek{i}.nev,cspe,cspv,cspb{j}{3});
            fprintf(' vege nincs meg -> amoba');
        end
        k=1;
        while ((k<=length(csp_nevsor)) && (~strcmp(cspv,csp_nevsor{k}))), k=k+1; end
        tjcs(k,i)=1;  gcsp_tipus(k) =2;
        iii=1;
        while (iii<=length(amoba) && (~strcmp(cspv,amoba{iii}{1}))), iii=iii+1; end
        tjcs(k,n_cso+n_mar+1)=amoba{iii}{2};
        tjcs_sor = [tjcs_sor k];
        tjcs_oszlop = [tjcs_oszlop n_cso+n_mar+1];
        amoba_no = [amoba_no iii];
    end
    
    if debug>2
        fprintf('\n---------------------------------------------------');
    end
end

if debug<1
    fprintf('OK');
end

if debug>2
    fprintf('\n\nA tjcs matrix:\n\n');
    fprintf('   csp.nev   tipus    |');
    for i=1:length(csovek), fprintf(' cso%2d',i); end
    fprintf(' | '); for i=1:length(mar), fprintf(' mar%2d ',i); end
    fprintf('| amoba fogy.|');
    fprintf('\n----------------------+');
    for i=1:length(csovek), fprintf('------'); end
    fprintf('-+-'); for i=1:length(mar), fprintf('-------'); end
    fprintf('+------------+');
    for i=1:length(csp_nevsor)
        fprintf('\n%9s ',csp_nevsor{i});
        if gcsp_tipus(i)==1, fprintf('  rug<->mar |'); end
        if gcsp_tipus(i)==2, fprintf('    amoba   |'); end
        for j=1:n_cso
            if ~tjcs(i,j)==0, fprintf('  %+2d  ',tjcs(i,j));
            else  fprintf('      '); end
        end
        fprintf(' | ');
        for j=1:n_mar
            if ~tjcs(i,n_cso+j)==0, fprintf('   %2d  ',tjcs(i,n_cso+j));
            else  fprintf('       '); end
        end
        fprintf('|');
        if gcsp_tipus(i)==2, fprintf(' %5.3e  |',tjcs(i,n_cso+n_mar+1));
        else fprintf('            |'); end
    end
    fprintf('\n----------------------+');
    for i=1:length(csovek), fprintf('------'); end
    fprintf('-+-'); for i=1:length(mar), fprintf('-------'); end
    fprintf('+------------+');
end

for i=1:length(csovek), info(csovek{i},3,wdir.path); end
for i=1:length(mar),    info(mar{i},3,t);  end

% A cso_mar matrix j-edik sora megmondja, hogy ahhoz, hogy a j-edik
% rugalmas csovel lepjek, melyik merev alrendszerektol
% kell info. Az elso ket szam mar es csp, a 3. az
% esetleges amoba csp sorszama.
for i=1:n_cso
    cso_mar{i}{1}=[0,0,0]; cso_mar{i}{2}=[0,0,0];
end

for i=1:n_cso
    for j=1:length(gcsp_tipus)
        if ~(tjcs(j,i)==0)
            for k=1:n_mar
                if ~(tjcs(j,n_cso+k)==0)
                    if tjcs(j,i)<0
                        cso_mar{i}{1}=[k, tjcs(j,n_cso+k)];
                    else
                        cso_mar{i}{2}=[k, tjcs(j,n_cso+k)]; 
                    end
                end
            end
        end
    end
end

% A mar_cso matrix j-edik sora megmondja, hogy ahhoz, hogy a j-edik
% merev alrendszerrel lepjek, honnan kell info.
mar_cso=[];
for i=1:n_mar
    l=1;
    for j=1:length(gcsp_tipus)
        if ~(tjcs(j,n_cso+i)==0)
            mar_cso{i}{l}{1}=tjcs(j,n_cso+i);
            melyikcso=[];
            for k=1:n_cso
                if ~(tjcs(j,k)==0), melyikcso=[melyikcso, k*tjcs(j,k)]; end
            end
            mar_cso{i}{l}{2}=melyikcso;
            l=l+1;
        end
    end
end

if debug>1;
    %        disp(tjcs);
    fprintf('\n\nA rendszer a rugalmas csovek szemszogebol (cso_mar):');
    for i=1:n_cso
        fprintf('\n cso:%d,   kapcsolodo mar.ek:',i);
        if ~cso_mar{i}{1}(1)==0,  fprintf('  eleje:  %d/%d  ',cso_mar{i}{1}(1),cso_mar{i}{1}(2));
        else fprintf('  eleje: amoba '); end
        if ~cso_mar{i}{2}(1)==0,  fprintf(' vege:   %d/%d',cso_mar{i}{2}(1),cso_mar{i}{2}(2));
        else fprintf(' age:  amoba'); end
    end
    
    fprintf('\n\nA rendszer a merev alrendszerek szemszogebol (mar_cso):');
    for i=1:length(mar_cso)
        fprintf('\n mar:%d,  kapcsolodo csomopontok szama: %d',i,length(mar_cso{i}))
        for j=1:length(mar_cso{i})
            fprintf('\n\t   csp: %d    csovek:',mar_cso{i}{j}{1});
            for k=1:length(mar_cso{i}{j}{2}), fprintf(' %d,', mar_cso{i}{j}{2}(k));end
        end
    end
end

% Inicializalas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kiiratas idokoze

if ~isfield(options,'name')
    options.name = 'dt_save';
    options.value = 'auto';
end

for i = 1:length(options)
    if strcmp(options(i).name,'dt_save')
        dtsave = options(i).value;
    else
        dtsave = 'auto';
    end
end

if ~strcmp(dtsave,'auto')
    dtki = str2double(dtsave);
else
    dtki = 0.0;
end

% merev alrendszer

for i = 1:length(mar)
    mar{i}.dtki = dtki;
    mar{i}.dtkiorig = dtki;
end
% rugalmas alrenszer
for i = 1:length(csovek)
    csovek{i}.dtki = dtki;
end

% Szamitas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if debug>1
    fprintf('\n\n****************************************************');
    fprintf('\n*                      SZAMITAS                    *');
    fprintf('\n****************************************************\n');
end
if debug<1, fprintf('\n\nSzamitas\n'); end

while t<tmax
    [dt,icso]=dtuj(csovek,t);
    
    if debug>1
        fprintf('\n\n\n     *****************************************************');
        fprintf('\n     *   %6d.lepes: t=%7.4e dt=%7.4e         *',lepes,t,dt);
        fprintf('\n     *****************************************************\n');
        fprintf('\nEzekkel a csovekkel kell lepnunk: '); disp(icso);
    end
    
    % Amoba fogyasztasok frissitese
    
    if length(amoba_no) > 0
        amoba = update_amoba_fogy(amoba,t);
    end
    
    for i = 1:length(tjcs_sor)
        kk = tjcs_sor(i);
        mm = tjcs_oszlop(i);
        tjcs(kk,mm) = amoba{amoba_no(i)}{2};
    end
     
    mars_updated={};
    for i=1:length(icso)
        acso=icso(i);
        %temp=csovek{acso}.csp;
        
        if debug>1, fprintf('\n\nrugalmas elem: %d\n-----------------\n',acso); end
        
        % Csatlakozo merev alrendszerek frissitese:
        %-------------------------------------------------
        if debug>1, fprintf('\n  A kapcsolodo merev alrendszerek frissitese:\n'); end
        
        % cso eleje:
        amare = cso_mar{acso}{1}(1);
        acspe = cso_mar{acso}{1}(2);
        
        if ~(amare==0)
            % csak ha meg nem leptunk ezzel a merev alrendszerrel:
            if sum(strcmp(mars_updated,num2str(amare)))==0
                mar = update_mar(amare,mar_cso,mar,csovek,'e',t,dt,dtmin,dtki,debug);
                mars_updated{end+1}=num2str(amare);
            else
                if debug>1, fprintf('\n\tEbben a korben mar leptem a(z) %s merev alrendszerrel.\n',mar{amare}.nev); end
            end
        else
            if debug>1, fprintf('\n\tAmoba csomopont a rugalmas cso elejen, nem kell frissiteni.\n'); end
        end
        
        % cso vege:
        amarv = cso_mar{acso}{2}(1);
        acspv = cso_mar{acso}{2}(2);
        
        if ~(amarv==0)
            % csak ha meg nem leptunk ezzel a merev alrendszerrel:
            if sum(strcmp(mars_updated,num2str(amarv)))==0
                mar = update_mar(amarv,mar_cso,mar,csovek,'v',t,dt,dtmin,dtki,debug);
                mars_updated{end+1}=num2str(amarv);
            else
                if debug>1, fprintf('\n\tEbben a korben mar leptem a(z) %s. merev alrendszerrel.\n',mar{amarv}.nev); end
            end
        else
            if debug>1, fprintf('\n\tAmoba csomopont a rugalmas cso vegen, nem kell frissiteni.\n'); end
        end
        
%         if t > 0.3
%             debug = 5;
%         end
        
        % Rugalmas elem leptetese:
        %-------------------------------------------------
        
        if debug>1, fprintf('\n Es most johet a %s (%d.) cso frissitese:',csovek{acso}.nev,acso);end
        temp=csovek{acso}.csp;
        
        pf{1} = update_rug_elem('eleje',gcsp_tipus(temp(1)),mar,amare,acspe,cspb,csovek,csp_nevsor,acso,tjcs,t,dt,dtmin,debug);
        pf{2} = update_rug_elem('vege',gcsp_tipus(temp(2)), mar,amarv,acspv,cspb,csovek,csp_nevsor,acso,tjcs,t,dt,dtmin,debug);
               
        csovek{acso} = solve(csovek{acso},pf);
       
        if csovek{acso}.t > csovek{acso}.dtki
            csovek{acso}.dtki = csovek{acso}.dtki + dtki;
            info(csovek{acso},3,wdir.path);
        end

    end
    
    if t>tt
        fprintf('\n%6d.lepes: t=%7.3e / %7.3e  dt=%5.3e',lepes,t,tmax,dt);
        tt=tt+dtout;
    end
    
    t=t+dt; lepes=lepes+1;

    if debug>1, pause; end
    fclose('all');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\nAz utolso lepes eredmenyei:');
for i=1:length(mar)
    fprintf('\n\n  => merev alrendszer: %s',mar{i}.nev);
    for j=1:length(mar{i}.csp)
        temp=mar{i}.csp;
        pp=get(mar{i},'p',temp{j}{1});
        fprintf('\n         %8s (%3d): p=%+6.3f [bar] = %5.2f vom',temp{j}{6},j,pp/1e5,pp/1000/9.81);
    end
    fprintf('\n');
    for j=1:length(mar{i}.elemek)
        temp=mar{i}.elemek;
        Q=get(mar{i},'Q',j);
        m=get(mar{i},'m',j);
        dp=get(mar{i},'dp',j);
        fprintf('\n         %8s (%3d): Q=%+5.3f[m3/s]   m=%+5.2f[kg/s]   dp=%+5.3f[bar]  dH=%+6.2f[m]',temp{j}.nev,j,Q,m,dp/1e5,dp/9.81/1000);
    end
end


fprintf('\n\n => rugalmas csovek:\n');
fprintf('\n    Nev      |        eleje         |        vege          |          ');
fprintf('\n             |   p[bar]    Q[m3/s]  |   p[bar]    Q[m3/s]  |  m[kg/s] ');
fprintf('\n-------------+----------------------+----------------------+-----------');

for i=1:length(csovek)
    ppp=csovek{i}.p;
    qqq=(csovek{i}.v).*(csovek{i}.A);
    nnn=csovek{i}.N;
    mmm=csovek{i}.ro*qqq;
    fprintf('\n   %7s  |   %6.2f    %+6.4f  |   %6.2f    %+6.4f  |  %+6.2f ',csovek{i}.nev,ppp(1)/1e5,qqq(1),ppp(nnn+1)/1e5,qqq(nnn+1),mmm(nnn+1));
end
% CPU ido
cput=cputime-cput; th=floor(cput/3600); cput=cput-3600*th; tm=floor(cput/60); cput=cput-60*tm;
fprintf('\n\nszamitasi ido:   %gh  %gmin  %gs\n\n',th,tm,cput);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fp=fopen('eredmenyek.out','w');
c = clock;
fprintf(fp,'Tranziens szimulator\n------------------------------------');

fprintf(fp,'\n A szamitas vege: %d.%d.%d  %d:%d',c(1),c(2),c(3),c(4),c(5));
fprintf(fp,'\n Adatfile: %s',df);
fprintf(fp,'\n\n Az utolso lepes eredmenyei:');
for i=1:length(mar)
    fprintf(fp,'\n\n  => merev alrendszer: %s',mar{i}.nev);
    for j=1:length(mar{i}.csp)
        temp=mar{i}.csp;
        pp=get(mar{i},'p',temp{j}{1});
        fprintf(fp,'\n         %8s (%3d): p=%+6.3f [bar]',temp{j}{6},j,pp/1e5);
    end
    fprintf('\n');
    for j=1:length(mar{i}.elemek)
        temp=mar{i}.elemek;
        Q=get(mar{i},'Q',j);
        m=get(mar{i},'m',j);
        dp=get(mar{i},'dp',j);
        fprintf(fp,'\n         %8s (%3d): Q=%+5.3f [m3/s]\t m=%+5.2f [kg/s]\t dp=%+5.3f [bar]\t dH=%+6.2f [m]',temp{j}.nev,j,Q,m,dp/1e5,dp/9.81/1000);
    end
end

fprintf(fp,'\n\n => rugalmas csovek:\n');
fprintf(fp,'\n    Nev      |        eleje         |        vege          |          ');
fprintf(fp,'\n             |   p[bar]    Q[m3/s]  |   p[bar]    Q[m3/s]  |  m[kg/s] ');
fprintf(fp,'\n-------------+----------------------+----------------------+-----------');
for i=1:length(csovek)
    ppp=csovek{i}.p;
    qqq=csovek{i}.v.*csovek{i}.A;   
    mmm=csovek{i}.ro*qqq;
    fprintf(fp,'\n   %7s  |   %6.2f    %+6.4f  |   %6.2f    %+6.4f  |  %+6.2f ',csovek{i}.nev,ppp(1)/1e5,qqq(1),ppp(end)/1e5,qqq(end),mmm(end));
end
fprintf(fp,'\n\nszamitasi ido:   %gh  %gmin  %gs\n\n',th,tm,cput);
fclose(fp);

%----------------------------------------------------------------------

function [dt,icso]=dtuj(csovek,t)

global dtmin

dt=1e5;
flag = 0;
for i=1:length(csovek)
    if isa(csovek{i},'viszkcso')
        dtm = dtviszko(csovek{i});
        dtt = csovek{i}.t+dtm-t;
        flag = 1;
    else
        dtt = csovek{i}.t+csovek{i}.dt-t;
    end
    temp(i,1:2) = [i dtt];
    if dtt<dt
        dt = dtt; 
        icso = i; 
    end
end

%Ha valamelyik cso ehhez nagyon kozel van, azzal is lepunk (flag: viszkcso
%eseten ezt kihagyjuk!)
if flag == 0
    for i=[1:icso-1,icso+1:length(csovek)]
        dtt=csovek{i}.t+csovek{i}.dt-t;
        if abs(dtt-dt) < dtmin
            icso=[icso,i];
        end
    end
elseif flag == 1
    for i=[1:icso-1,icso+1:length(csovek)]
    if isa(csovek{i},'viszkcso')
        dtm = dtviszko(csovek{i});
        dtt = csovek{i}.t+dtm-t;
    else
        dtt = csovek{i}.t+csovek{i}.dt-t;
    end
        if abs(dtt-dt) < dtmin
            icso=[icso,i];
        end
    end
end

%[Y I] = sort(temp);
% mask = find(temp(:,2) < dtmin & temp(:,1) ~= icso & temp(:,2) < temp(icso,2)*2);
% temp2 = I(mask,2).';
% icso=[icso,temp2];

%fprintf('\nt=%7.5e  dt1=%7.5e  dt2=%7.5e  t1=%7.5e  t2=%7.5e',t,csovek{1}.t+csovek{1}.dt-t,csovek{2}.t+csovek{2}.dt-t,csovek{1}.t,csovek{2}.t);
%fprintf(' -> dt=%7.5e, csovek: %d ==>',dt,length(icso));
%for i=1:length(icso), fprintf('  %g  ',icso(i)); end

%----------------------------------------------------------------------

function out = dtviszko(viszkcso)

% data = viszkcso.all;
% x = data{1}.';
% v = data{2}.';
% a = data{3}.';
% t = data{4}.';
% N = data{5}.';
% L = data{6}.';
% 
% for i = 2:N
%     xL = x(i-1); xR = x(i+1);
%     vL = v(i-1); vR = v(i+1);
%     aL = a(i-1); aR = a(i+1);
% 
%     tP = (xR - xL - t*(vR - aR) + t*(vL + aL))/((vL + aL) - (vR - aR));
%     dt(i) = tP - t;
% end
% 
% tP1 = t - x(2)/(v(2) - a(2));
% dt(1) = tP1 - t;
% tPN1 = t + (L - x(N))/(v(N) + a(N));
% dt(N+1) = tPN1 - t;
% dtm = min(dt);
dtm = viszkcso.dtuj;

out = dtm;
    
%----------------------------------------------------------------------

function out = update_mar(amar,mar_cso,mar,csovek,eleje_v_vege,t,dt,dtmin,dtki,debug)

% amar : a merev alrendszer szama, amit frissiteni akarunk

if ~(amar==0) % ha van egyaltalan merev alr.
    
    % Elso lepeskent be kell kerni a amar-hoz kapcsolodo rugalmas elemek
    % peremfelteteleit: pf-be gyujtjuk
    l=1; pf={};
    for j=1:length(mar_cso{amar})
        for k=1:length(mar_cso{amar}{j}{2})
            ezacso=mar_cso{amar}{j}{2}(k);
            %temp=csovek{abs(ezacso)}.csp;
            
            if eleje_v_vege=='e'
                pf_tipus=csovek{abs(ezacso)}.pf_eleje_tipus;
            else
                pf_tipus=csovek{abs(ezacso)}.pf_vege_tipus;
            end
            
            switch pf_tipus
                case 'karakterisztika'
                    if ezacso<0
                        % a rugalmas cso ELEJEn van a mar
                        ezacso=abs(ezacso);
                        %pf_tipus = csovek{ezacso}.pf_eleje_tipus;
                        if isa(csovek{ezacso}, 'viszkcso')
                            pf{l}={mar_cso{amar}{j}{1},'viszkcso',karm(csovek{ezacso},t+dt-csovek{ezacso}.t,dtmin)};
                            
                         
                        elseif isa(csovek{ezacso}, 'cso')
                            pf{l}={mar_cso{amar}{j}{1},'cso',karm(csovek{ezacso},t+dt-csovek{ezacso}.t,dtmin)};
                        end
                        kar_jel='-';
                    else
                        % a rugalmas cso VEGEn van a mar
                        %pf_tipus = csovek{ezacso}.pf_vege_tipus;
                        if isa(csovek{ezacso}, 'viszkcso')
                            pf{l}={mar_cso{amar}{j}{1},'viszkcso',karp(csovek{ezacso},t+dt-csovek{ezacso}.t,dtmin)};
                        elseif isa(csovek{ezacso}, 'cso')
                            pf{l}={mar_cso{amar}{j}{1},'cso',karp(csovek{ezacso},t+dt-csovek{ezacso}.t,dtmin)};
                        end
                        kar_jel='+';
                    end
                    if debug>1
                        fprintf('\n\t %d. pf -> csp:%d  cso: %s (%g.) C%s kar., dt=%+5.3e',...
                            l,mar_cso{amar}{j}{1},csovek{ezacso}.nev,ezacso,kar_jel,t+dt-csovek{ezacso}.t);
                    end
                    
                case 'vizszint_&_konti'
                    if ezacso<0
                        % a csatorna ELEJEn van a mar
                        ezacso=abs(ezacso);
                        %pf_tipus = csovek{ezacso}.pf_eleje_tipus;
                        pf{l}={mar_cso{amar}{j}{1},'vizszint_&_konti',karm(csovek{ezacso},t+dt-csovek{ezacso}.t,dtmin)};

                    else
                        % a csatorna VEGEn van a mar
                        %pf_tipus = csovek{ezacso}.pf_vege_tipus;
                        pf{l}={mar_cso{amar}{j}{1},'vizszint_&_konti',karp(csovek{ezacso},t+dt-csovek{ezacso}.t,dtmin)};
                    end
                    
                otherwise
                    error('Hibas peremfeltetel: %s',pf_tipus);
            end
            
            % pf{.} leptetese
            l=l+1;
        end
    end
    
    % Frissites:
    if debug>1
        fprintf('\n -->Lepek egyet a %d. merev alrendszerrel: %5.3e [s] -> ',...
            amar,mar{amar}.t);
    end
    
    mar{amar} = solve(mar{amar},pf,t+dt);
    
    out = mar;

    if debug>1, fprintf(' %5.3e [s]\n',mar{amar}.t); end
    
end

%----------------------------------------------------------------------

function pf = update_rug_elem(flag,gcsp_tipus,mar,amar,acsp,cspb,csovek,csp_nevsor,acso,tjcs,t,dt,dtmin,debug)

n_cso=length(csovek);
n_mar=length(mar);

if strcmp(flag,'eleje')
    ev_string = 'Eleje';
else
    ev_string = 'Vege';
end

switch gcsp_tipus
    case 1 % Merev alrendszer kapcsolodik      
        pf={'p',get(mar{amar},'p',acsp)}; % abszolut nyomas, nullszinthez viszonyitva!!!
        if debug>2
            tmp=mar{amar}.csp;
            fprintf('\n\t %s: %s (%d.) merev alrendszer %s (%d.) csomopontja',ev_string,mar{amar}.nev,amar,tmp{acsp}{6},acsp);
        end
        
    case 2 % Csak rugalmas csovek talalkoznak
        if debug>1, fprintf('\n\n\t A cso %sn amoba csomopont van',lower(ev_string)); end
        temp=csovek{acso}.csp;
        if strcmp(flag,'eleje')
            acsp=temp(1);
        else
            acsp=temp(2);
        end
        if debug>2, fprintf('\n\t\t\t ezek a rugalmas csovek talalkoznak a %s (%d) csomopontban:',csp_nevsor{acsp},acsp); end
        
        % csovek kivalogatasa
        gcso=[];       
        for j=1:n_cso
            if ~(tjcs(acsp,j)==0), gcso=[gcso,j*tjcs(acsp,j)]; end
        end
        if debug>1, disp(gcso); end
        
        % Csomoponti konti+kar. egyenletek megoldasa
        if debug>2, fprintf('\t Amoba egyenletek megoldasa:'); end
        szum1=0; szum2=0;
        for j=1:length(gcso)
            if gcso(j)>0
                ezacso=gcso(j);
                out = karp(csovek{ezacso},t+dt-csovek{ezacso}.t,dtmin);
                if debug>2, fprintf('\n\t\t   %d. rugalmas elem interpolaljon C+ menten ide: %7.5f[s]  innen: %7.5f[s]  azaz dt: %7.5f[s], val=%7.5f, p=%7.5f, v=%7.5f',ezacso,t+dt,csovek{ezacso}.t,t+dt-csovek{ezacso}.t,out(5)/1e5,out(6)/1e5,out(7)); end
            else
                ezacso=abs(gcso(j));
                out = karm(csovek{ezacso},t+dt-csovek{ezacso}.t,dtmin);
                if debug>2, fprintf('\n\t\t   %d. rugalmas elem interpolaljon C- menten ide: %7.5f[s]  innen: %7.5f[s]  azaz dt: %7.5f[s], val=%7.5f, p=%7.5f, v=%7.5f',ezacso,t+dt,csovek{ezacso}.t,t+dt-csovek{ezacso}.t,out(5)/1e5,out(6)/1e5,out(7)); end
            end
            szum1=szum1+out(4)*out(3)/out(2); szum2=szum2+out(4)/out(2);
        end
        
        Qbe = tjcs(acsp,n_mar+n_cso+1);
        p = (szum1-Qbe)/szum2;
        
        if debug>1, fprintf('\n\t\tAmoba csomoponti nyomas:  %g [bar],  Qbe: %g [m^3/s]',p/1e5,Qbe); end
        pf={'p',p};
        
    case 3 % Adott nyomasu pont
        pf={'p',gcsp_p(temp(1))};
        
    otherwise
        error('Megszivtad:  gcsp_tipus=%s', gcsp_tipus);
end

%% Amoba fogyasztasok frissitese
function amoba = update_amoba_fogy(amoba,t)

temp = amoba{1}{6};
for i=1:length(temp)
    gorbe_nevsor{i}=temp{i}.nev;
end

for i=1:length(amoba)
    tmp=strcmp(amoba{i}{4},gorbe_nevsor);
    if sum(tmp)==0
        error(['trsz.m, update_amoba_fogy: nem talalom a ',amoba{i}{4},' nevu lefutasgorbet.']);
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
        amoba{i}{2}=amoba{i}{5}*interp1(temp{ezaz}.x,temp{ezaz}.y,abs(t));
    end
end
