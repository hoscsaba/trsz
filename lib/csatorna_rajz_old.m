function csatorna_rajz(varargin)

clc

mdt  = varargin{1};
rug  = varargin{2};
teleje = varargin{3};
tvege = varargin{4};
dt  = varargin{5};
wdir = what;

j = 0;
for i = 6:nargin
    j = j+1;
    elemnev{j} = varargin{i};
end

element_found=zeros(size(elemnev));

numelem = j;
debug = 0;
mar = olvas_merev(mdt,debug,wdir.path);

rug_struct = olvas_rugalmas(rug,debug);
csovek     = rug_struct.csovek;
csp_nevsor = rug_struct.cpnevsor;
csom       = rug_struct.csom;
if rug_struct.amoba_exists
    amoba  = rug_struct.amoba;
end

%% Aknak keresese

marszam = 0;
maxakna = 0;
for jj = 1:numelem
    megvan = 0;
    for i = 1:length(mar)
        for j = 1:length(mar{i}.elemek)
            if strcmp(mar{i}.elemek{j}.nev,elemnev{jj})
                fprintf('%s megvan!\n',mar{i}.elemek{j}.nev);
                mar_num = i;
                akna_num = j;
                megvan = 1;
            end
        end
    end
    
    if megvan == 0
        continue
    else
        marszam = marszam+1;
        element_found(jj)=1;
    end
    
    data = load(strcat(wdir.path,'\',mar{mar_num}.nev,'.res'));
    n_ag = length(mar{mar_num}.elemek);
    n_csp = length(mar{mar_num}.csp);
    akna_talp(marszam) = mar{mar_num}.elemek{akna_num}.hmin;
    
    % Meg kell szamolni a szivattyukat es a nyomasszabalyzokat
    temp=mar{mar_num}.elemek;
    n_sziv = 0;
    n_nysz = 0;
    for j=1:length(temp)
        if isa(temp{j},'szivattyu'), n_sziv = n_sziv+1; end
        if isa(temp{j},'nyomasszabalyzo'), n_nysz = n_nysz+1; end
    end
    
    eltol = 0;
    temp = mar{mar_num}.elemek;
    kk = 0;
    
    for i=1:n_ag
        if isa(temp{i},'akna')
            kk=kk+1;
            if strcmp(temp{i}.nev,elemnev{jj})
                eltol=kk;
            end
        end
    end
    
    adat_hely = 1+2*n_ag+n_csp+n_sziv+n_nysz+eltol;
    
    takt = teleje;
    k = 1;
    while takt < tvege+dt
        j = 1;
        while data(j,1) <= takt
            elotte = j;
            utana = j+1;
            j = j+1;
        end
        t_elotte = data(elotte,1);
        t_utana = data(utana,1);
        y_elotte = data(elotte,adat_hely);
        y_utana = data(utana,adat_hely);
        
        akna_y{marszam}(k) = interp1([t_elotte t_utana],[y_elotte; y_utana],takt);
        
        takt = takt + dt;
        k = k+1;
    end
end

%% Rugalmas elemek keresese

x = [];
yy = [];
zf = [];
zt = [];
rugszam = 0;
csatl(1) = 0;
for jj = 1:numelem
    megvan = 0;
    for i = 1:length(csovek)
        if strcmp(csovek{i}.nev,elemnev{jj})
            fprintf('%s megvan!\n',csovek{i}.nev);
            csat_num = i;
            megvan = 1;
        end
    end
    
    if megvan == 0
        continue
    else
        rugszam = rugszam+1;
        element_found(jj)=1;
    end
    
    data = load(strcat(wdir.path,'\',csovek{csat_num}.nev,'.res'));
    csat = csovek{csat_num};
    
    takt = teleje;
    k = 1;
    maxy(rugszam) = 0;
    while takt < tvege+dt
        j = 1;
        while data(j,1) < takt
            elotte = j;
            utana = j+1;
            j = j+1;
        end
        t_elotte = data(elotte,1);
        t_utana = data(utana,1);
        y_elotte = data(elotte,10:10+length(csat.y)-1);
        y_utana = data(utana,10:10+length(csat.y)-1);
        if max(y_elotte) > maxy(rugszam)
            maxy(rugszam) = max(y_elotte);
        end
        
        y_intp{rugszam}{k} = interp1([t_elotte t_utana],[y_elotte; y_utana],takt);
        if rugszam > 1
            yy(k) = {[yy{k} y_intp{rugszam}{k}]};
        else
            yy{k} = y_intp{rugszam}{k};
        end
        t{rugszam}(k) = takt;
        takt = takt + dt;
        k = k+1;
    end
    
    maxk = length(t{rugszam});
    xx = linspace(0,csat.L,(csat.N+1));
    zzf = linspace(csat.ze,csat.zv,csat.N+1);
    zzt = zzf+csat.dvB;
    
    if rugszam > 1
        xx = xx + max(x);
    end
    
    x = [x xx];
    zf = [zf zzf]; 
    zt = [zt zzt];
    csatl(rugszam+1) = max(x);
end
fclose all;

for i=1:length(element_found)
    if element_found(i)==0
        error([elemnev{jj},' csatornat vagy aknat nem talalom!']);
    end
end

%% Rajz indul

for k = 1:marszam
    akna_y{k} = akna_y{k} - min(akna_talp);
    if max(akna_y{k}(:))> maxakna
        maxakna = max(akna_y{k}(:));
    end
end

zf = zf - min(akna_talp);
zt = zt - min(akna_talp);

f1 = figure(1);
screen_size = get(0, 'ScreenSize');
set(f1, 'Position', [50 50 0.8*screen_size(3) 0.8*screen_size(4) ] );
for i = 1:maxk

    %h2 = area(x,[z.' y_intp{1}{i}.']);
    
    %h2 = area(x,[yy{i}.' z.']);
    %h2 = area(x, [z.' yy{i}.']);
    %h2 = area(x, z.');
    
    % fenek
    hf = plot(x,zf,'k-'); hold on
    % teteje
    ht = plot(x,zt,'k-'); hold on
    % vizszint
    hv = plot(x,(yy{i}.'+zf.').','b-'); hold on   
    % akna
    for j = 1:marszam
        ha = plot(csatl(j),akna_y{j}(i),'r*');
        %set(h3,'MarkerFaceColor','r')
        %set(h3,'MarkerSize',8)
        %line([csatl(j) csatl(j)],[0 1e5],'Color','r')
    end
    hold off, grid on
    
    set(hf,'LineWidth',2);
    set(ht,'LineWidth',2);
    set(hv,'LineWidth',2);
    
    %set(h2(1),'FaceColor',[0.2 0.4 0.3])
    %set(h2(2),'FaceColor',[0.0 0.0 1.0])
    %set(h2(1),'LineWidth',2)
    %set(h2(2),'LineWidth',2)
    
    ttl =title(['Csatornaszintek a ' , csovek{csat_num}.nev , ' csatornaban, t = ' , num2str(t{1}(i)/60, '%5.2f') , ' min']);
    set(ttl,'FontSize',16)
    xlabel('Hossz [m]','FontSize',14)
    ylabel('Szint [m]','FontSize',14)
    %maxylim = max( max(maxy)+max(z) , maxakna );
    %ylim([0 1.1*maxylim]);
    %F(i) = getframe(gcf);
    pause
end

% fileout = strcat(csovek{csat_num}.nev,'.avi');
% movie2avi(F,fileout);
