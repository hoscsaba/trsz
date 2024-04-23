function mar=olvas5m(fnev,debug)
%...merev-adatfile olvasasa..............
%global cpnevsor ncp

grav=9.806;

if debug>1
    fprintf('\n\n****************************************************');
    fprintf('\n*          MEREV ALRENDSZEREK OLVAS�SA             *');   
    fprintf('\n****************************************************\n\n');  
end

if debug>0
    fprintf('\n  Feladat neve: %s',fnev);
    fprintf('\n  Adatfile megnyit�sa...')
end

%...file-nev es kiterjesztes...
kiterjeszt='.mdt';
fnevdat=strcat(fnev,kiterjeszt);
if debug>0,  fprintf('\n  Adatok olvas�sa...');end

% be string vektorba beolvassuk az eg�sz file-t.
[be]=textread(fnevdat,'%s','delimiter',',');

%...k: cso-index...j: jelleggorbe-index...jc: csomopont index...
%..feladat neve
fnev=char(be(1));
ncp=str2num(char(be(2)));       %....beolvasand� merev alrendszerek sz�ma...   
eltol=3;

mar_szamlalo=0; 
for jc=1:ncp
    ag_szamlalo=0;
    mar_szamlalo=mar_szamlalo+1;
    cspnev=char(be(eltol));             %....rugalmas-csomopont neve
    nag=str2num(char(be(eltol+1)));       %....agak szama a csp-ben..   
    if debug>1, fprintf('\n\n%2d./%2d csomoponti alrendszer neve: %s  (%d �g)',mar_szamlalo,ncp,cspnev,nag); end
    
    %...az agadatok ....
    eltol=eltol+2;
    for k=1:nag
        ag_szamlalo =ag_szamlalo+1;	
        tipus{k}=char(be(eltol));
        switch tipus{k}
            case 'konc_cso'
                agnev{k}=char(be(eltol+1));                 %...agnev..
                csp1{k}=char(be(eltol+2));                  %...csp elejen
                csp2{k}=char(be(eltol+3));                  %...csp vegen
                D(k)=str2num(char(be(eltol+4)));            %...atmero [m]
                L(k)=str2num(char(be(eltol+5)));            %...aghossz [m]
                lam(k)=str2num(char(be(eltol+6)));          %...lambda
                ro(k)=str2num(char(be(eltol+7)));           %...suruseg [kg/m3]
                m0(k)=str2num(char(be(eltol+8)));            %...tomegaram [kg/s]
                eltol=eltol+9;
                if debug>1, fprintf('\n   %2d. elem: konc_cso (%s)',ag_szamlalo,agnev{k}); end

            case 'fojtas'
                agnev{k}=char(be(eltol+1));                 %...agnev..
                csp1{k}=char(be(eltol+2));                  %...csp elejen
                csp2{k}=char(be(eltol+3));                  %...csp vegen
                K1(k)=str2num(char(be(eltol+4)));           %...ellenallas tenyezo [1/m2]
                ro(k)=str2num(char(be(eltol+5)));           %...suruseg [kg/m3]
                m0(k)=str2num(char(be(eltol+6)));            %...tomegaram [kg/s]
                eltol=eltol+7;
                if debug>1, fprintf('\n   %2d. elem: fojtas (%s)',ag_szamlalo,agnev{k}); end

            case 'vez_fojtas'
                agnev{k}=char(be(eltol+1));                   %...agnev..
                csp1{k}=char(be(eltol+2));                    %...csp elejen
                csp2{k}=char(be(eltol+3));                    %...csp vegen
                ro(k)=str2num(char(be(eltol+4)));             %...suruseg [kg/m3]
                A(k)=str2num(char(be(eltol+5)));              %...nevleges keresztmetszet [m2]
                jgpsz(k)=str2num(char(be(eltol+6)));          %...K(eps) jellegg. pontsz.
							      
		% Ez a r�gi					      
                %m0(k)=str2num(char(be(eltol+7)));            %...tomegaram [kg/s]
                %eltol=eltol+7;
                %for j=1:jgpsz(k)
                %    jj=2*j-1;
                %    eps(k,j)=str2num(char(be(eltol+jj)));             %...rel. helyzet [-]
                %    K(k,j)=str2num(char(be(eltol+jj+1)));             %...ellenallas tenyezo [1/m2]
                %end
                %eltol=eltol+2*jgpsz(k)+1;
		
		% Ez az uj:
		jgpszt(k)=str2num(char(be(eltol+7)));        %...eps(t) fv.
		m0(k)=str2num(char(be(eltol+8)));            %...tomegaram [kg/s]
		
                eltol=eltol+8;
                for j=1:jgpsz(k)
                    jj=2*j-1;
                    epsK(k,j)=str2num(char(be(eltol+jj)));             %...rel. helyzet [-]
                    K(k,j)   =str2num(char(be(eltol+jj+1)));           %...ellenallas tenyezo [1/m2]
%		    fprintf('\n\t %2d./%2d adat:  eps=%+5.3e  K=%5.3e',j,jgpsz(k),eps(k,j),K(k,j));
                end
                eltol=eltol+2*jgpsz(k);
		
		for j=1:jgpszt(k)
                    jj=2*j-1;
                    t_vf(k,j)=str2num(char(be(eltol+jj)));               %...ido [s]
                    epst(k,j)=str2num(char(be(eltol+jj+1)));             %...rel. helyzet [-]
%                    fprintf('\n\t %2d./%2d adat:  t=%+5.3e  eps=%5.3e',j,jgpszt(k),t_vf(k,j),epst(k,j));
		end
                eltol=eltol+2*jgpszt(k)+1;
		
		% Ez m�r k�z�s...
                if debug>1, fprintf('\n   %2d. elem: vez_fojtas (%s)',ag_szamlalo,agnev{k}); end
		
            case 'szivattyu'
                agnev{k}=char(be(eltol+1));                  %...agnev..
                csp1{k}=char(be(eltol+2));                   %...csp elejen
                csp2{k}=char(be(eltol+3));                   %...csp vegen
                ro(k)=str2num(char(be(eltol+4)));            %...suruseg [kg/m3]
                Ds(k)=str2num(char(be(eltol+5)));            %...D szivocsonk [m]
                Dn(k)=str2num(char(be(eltol+6)));            %...D nyomocsonk [m]
                tranziens(k)=str2num(char(be(eltol+7)));     %...tranziens jellemzo..[-]
                jgpsz(k)=str2num(char(be(eltol+8)));         %...K(eps) jellegg. pontsz.
                m0(k)=str2num(char(be(eltol+9)));            %...tomegaram [kg/s]
                eltol=eltol+9;
                for j=1:jgpsz(k)
                    jj=3*j-2;
                    Q(k,j)=str2num(char(be(eltol+jj)));               %...terfogataram[m3/s]
                    H(k,j)=str2num(char(be(eltol+jj+1))); %...sz�llitomagassag [m]
                    P(k,j)=str2num(char(be(eltol+jj+2)));             %...felvett teljesitmeny [kW]
                    %                    if tranziens(k)>0
                    %                        fprintf('\n      Q: %7.3f  H: %7.3f  P: %7.3f',Q(k,j),H(k,j),P(k,j));
                    %                    else
                    %                        fprintf('\n      Q: %7.3f  H: %7.3f',Q(k,j),H(k,j));
                    %                    end
                end
                eltol=eltol+3*jgpsz(k);
                switch tranziens(k)
                    case 0     %...allando fordulatszam...
                        eltol=eltol+1;
                    case 1     %...kifutas...
                        n(k)=str2num(char(be(eltol+1)));            %...fordulatszam [1/min]
                        teta(k)=str2num(char(be(eltol+2)));         %...forgoreszek tehetetlensegi nyomateka [kgm2]
                        tki(k)=str2num(char(be(eltol+3)));          %...kifutas idopontja [sec]
                        eltol=eltol+4;
                    case 2     %...inditas, kozelito motor-jelleggorbe...
                        n(k)=str2num(char(be(eltol+1)));            %...fordulatszam [1/min]
                        teta(k)=str2num(char(be(eltol+2)));         %...forgoreszek tehetetlensegi nyomateka [kgm2]
                        Mb(k)=str2num(char(be(eltol+3)));           %...motor billeno-nyomatek [Nm]
                        nb(k)=str2num(char(be(eltol+4)));           %...motor billeno-nyomatek fordulatszama [1/min]
                        nsz(k)=str2num(char(be(eltol+5)));          %...motor szinkron-fordulatszama [1/min]
                        eltol=eltol+6;
                    case 3     %...inditas, adott motor-jelleggorbe...
                        n(k)=str2num(char(be(eltol+1)));            %...fordulatszam [1/min]
                        teta(k)=str2num(char(be(eltol+2)));         %...forgoreszek tehetetlensegi nyomateka [kgm2]
                        jgpsz(k)=str2num(char(be(eltol+3)));        %...motor billeno-nyomatek [Nm]
                        eltol=eltol+3;
                        for j=1:jgpsz(k)
                            jj=2*j-2;
                            Mm(k,j)=str2num(char(be(eltol+jj)));          %...motor-nyomatek [Nm]
                            nm(k,j)=str2num(char(be(eltol+jj+1)));        %...motor fordulatszam [1/min]
                        end
                        eltol=eltol+2*jgpsz(k)+1;
                    otherwise
                        fprintf('\nhibas tranziens jel, cspnev: %8s agnev: %8s ',cspnev,agnev{k});
                end
                if debug>1, fprintf('\n   %2d. elem: szivattyu (%s)',ag_szamlalo,agnev{k}); end
            case 'nyomas'
                agnev{k}=char(be(eltol+1));                 %...agnev..
                csp1{k}=char(be(eltol+2));                   %...csp elejen
                csp2{k}='';
                pk(k)=str2num(char(be(eltol+3)));            %...adott nyomas [N/m2]
                eltol=eltol+4;
                if debug>1, fprintf('\n   %2d. elem: nyomas (%s)',ag_szamlalo,agnev{k}); end
                
            case 'visszacsapo_szelep'
                agnev{k}=char(be(eltol+1));                 %...agnev..
                csp1{k} =char(be(eltol+2));                   %...csp elejen
                csp2{k} =char(be(eltol+3));
                ro(k)   =str2num(char(be(eltol+4)));
                m0(k)   =str2num(char(be(eltol+5)));
                eltol=eltol+6;
                if debug>1, fprintf('\n   %2d. elem: visszacsapo_szelep (%s)',ag_szamlalo,agnev{k}); end
                
            case 'nyomasszabalyzo'
                agnev{k}=char(be(eltol+1));                 %...agnev..
                csp1{k} =char(be(eltol+2));                   %...csp elejen
                csp2{k} =char(be(eltol+3));
                ro(k)   =str2num(char(be(eltol+4)));
                m0(k)   =str2num(char(be(eltol+5)));
                jgpsz(k)=str2num(char(be(eltol+6)));
                eltol=eltol+6;
                for j=1:jgpsz(k)
                    jj=2*j-1;
                    e(k,j) =str2num(char(be(eltol+jj)));
                    ek(k,j)=str2num(char(be(eltol+jj+1)));
                end
                eltol=eltol+2*jgpsz(k);
                szab{k} = char(be(eltol+1));
                switch szab{k}
                    case 'p'
                        szabcsp{k}=       char(be(eltol+2));
                        ajel(k) = str2num(char(be(eltol+3)));
                        pmax(k) = str2num(char(be(eltol+4)));
                        szabP(k)= str2num(char(be(eltol+5)));
                        szabI(k)= str2num(char(be(eltol+6)));
                        szabD(k)= str2num(char(be(eltol+7)));
                        e0(k)   = str2num(char(be(eltol+8)));
                        vmax(k) = str2num(char(be(eltol+9)));
                        tbe(k)  = str2num(char(be(eltol+10)));
                        eltol=eltol+11;  
                    case 'dp'
                        szabcsp1{k}=       char(be(eltol+2));
                        szabcsp2{k}=       char(be(eltol+3));
                        ajel(k) = str2num(char(be(eltol+4)));
                        pmax(k) = str2num(char(be(eltol+5)));
                        szabP(k)= str2num(char(be(eltol+6)));
                        szabI(k)= str2num(char(be(eltol+7)));
                        szabD(k)= str2num(char(be(eltol+8)));
                        e0(k)   = str2num(char(be(eltol+9)));
                        vmax(k) = str2num(char(be(eltol+10)));
                        tbe(k)  = str2num(char(be(eltol+11)));
                        eltol=eltol+12;  
                    otherwise
                        error('Szab�lyz�szelep csak nyom�sra  vagy nyom�sk�l�nbs�gre tud szab�lyozni!');
                end
                
                if debug>1, fprintf('\n   %2d. elem: nyom�sszab�lyz� szelep (%s)',ag_szamlalo,agnev{k}); end 
                
            otherwise
                error(['\nhibas agtipus, cspnev: %8s',cspnev]);
            end
        end   %...k
        
        %....cpnevsor a csomopontban.......
        % Itt a nagy ciklus jc-je fel�l�r�dik, ami szeintem hib�s. A bels�
        % jc-ket lecser�lem jcc-re.
        cpnevsor{1}=csp1{1};
        jcc=1;
        for k=1:nag
            jel='n';
            for i=1:jcc
                if strcmp(cpnevsor{i},csp1{k})
                    jel='i';
                end;
            end;
            if strcmp(jel,'n')&~isempty(csp1{k})
                jcc=jcc+1;
                cpnevsor{jcc}=csp1{k};
            end
            jel='n';
            for i=1:jcc
                if strcmp(cpnevsor{i},csp2{k})
                    jel='i';
                end
            end
            if strcmp(jel,'n')&~isempty(csp2{k})
                jcc=jcc+1;
                cpnevsor{jcc}=csp2{k};
            end
        end
        
        ncp_m=jcc;  %...csomopontszam a merev alrendszerben...
        
        %...rangszam kiszamitasa es cpnevsor ellenorzese............. 
        % Itt szint�n lecser�ltem a jc-t jcc-re. HCs
        rang=zeros(1,ncp_m);
        for jcc=1:ncp_m 
            mar_volt=0; 
            for k=1:nag 
                if strcmp(cpnevsor{jcc},csp1{k}) |  strcmp(cpnevsor{jcc},csp2{k})   
                    rang(jcc)=rang(jcc)+1; 
                    mar_volt=1; 
                end; 
            end 
            if mar_volt==0 
                error(['\n cp hiba %8s',cpnevsor{jcc}]); 
            end 
        end 
        
        %rang 
        r_max=max(rang);
        %r_max
        
        %...csom osszeallitas... 
        % meg itt is jc->jcc
        csom=zeros(ncp_m,r_max); 
        for jcc=1:ncp_m 
            j=0; 
            for k=1:nag 
                if strcmp(cpnevsor{jcc},csp1{k})
                    j=j+1; 
                    csom(jcc,j)=-k; 
                end 
                
                if strcmp(cpnevsor{jcc},csp2{k})
                    j=j+1; 
                    csom(jcc,j)=k; 
                end 
            end 
        end           
        %csom        %...itt a csom vege.......
        
        %...csp magassag es elvetel beolvasas....
        for jcc=1:ncp_m
            cspnev1{jcc}=char(be(eltol));                    %...csomopontnev 
            jelzo='nincs';
            for i=1:ncp_m
                %fprintf('\n %8s  %8s',cspnev1{jcc},cpnevsor{i})
                if strcmp(cpnevsor{i},cspnev1{jcc})
                    z(i)=str2num(char(be(eltol+1)));                 %...magassag... [m]
                    elv(i)=str2num(char(be(eltol+2)));               %...elvetel... [kg/s]
                    eltol=eltol+3;
                    jelzo='van';
                end
            end
            if strcmp(jelzo,'nincs')
                fprintf('\n\n %8s nev� csom�pont nincs\n\n',cspnev1{jcc}); pause
            end
            
        end
        
        if debug>1
            fprintf('\n\n    csp. sz�m:  ');
            for jcc=1:ncp_m, fprintf('   %d   ',jcc); end
            fprintf('\n    csp. n�v: ');
            for jcc=1:ncp_m, fprintf(' %5s ',cpnevsor{jcc}); end
        end
        
        if debug>2
            fprintf('\n\n    jc  cpnevsor   z           elvetel');
            for jcc=1:ncp_m
                fprintf('\n   %3i  %8s  %8.3f  %8.3f',jcc,char(cpnevsor{jcc}),z(jcc),elv(jcc))
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ITT NY�LTAM BELE!!!!
        % Merev �gelemek l�trehoz�sa.
        
        if debug>1, fprintf('\n\n    elemek �p�t�se...'); end
        
        clear elemek;
        for k=1:nag 
            if debug>1, fprintf(' %d ',k); end
            for i=1:length(cpnevsor)
                if strcmp(cpnevsor{i},csp1{k}), cspe_szam=i; end
                if strcmp(cpnevsor{i},csp2{k}), cspv_szam=i; end
            end
            switch tipus{k}                        
                case 'konc_cso'
                    elemek{k}=konc_cso(agnev{k},cspe_szam,cspv_szam,D(k),L(k),lam(k),ro(k),m0(k));
                    
                case 'szivattyu'
                    switch tranziens(k)
                        case 0     %...allando fordulatszam...
                            elemek{k}=szivattyu(agnev{k},cspe_szam,cspv_szam,Q(k,:),H(k,:),Ds(k),Dn(k),ro(k),m0(k),tranziens(k));                        
                        case 1     %...kifutas...
                            elemek{k}=szivattyu(agnev{k},cspe_szam,cspv_szam,Q(k,:),H(k,:),Ds(k),Dn(k),ro(k),m0(k),tranziens(k),P(k,:),n(k),teta(k),tki(k));                        
                        case 2     %...inditas, kozelito motor-jelleggorbe...
                            elemek{k}=szivattyu(agnev{k},cspe_szam,cspv_szam,Q(k,:),H(k,:),Ds(k),Dn(k),ro(k),m0(k),tranziens(k),P(k,:),n(k),teta(k),Mb(k),nb(k),nsz(k));
                        case 3     %...inditas, adott motor-jelleggorbe...
                            elemek{k}=szivattyu(agnev{k},cspe_szam,cspv_szam,Q(k,:),H(k,:),Ds(k),Dn(k),ro(k),m0(k),tranziens(k),P(k,:),n(k),teta(k),Mm(k,:),nm(k,:));
                    end
                    
                case 'fojtas'                
                    elemek{k}=fojtas(agnev{k},cspe_szam,cspv_szam,1e-10,K1(k),ro(k),m0(k));
	     
	        case 'vez_fojtas' 
                    elemek{k}=vez_fojtas(agnev{k},cspe_szam,cspv_szam,ro(k),A(k),jgpsz(k),jgpszt(k),m0(k),epsK(k,:),K(k,:),t_vf(k,:),epst(k,:));
                    
                case 'nyomas'
                    elemek{k}=nyomas(agnev{k},cspe_szam,pk(k));
                    
                case 'visszacsapo_szelep'
                    elemek{k}=visszacsapo_szelep(agnev{k},cspe_szam,cspv_szam,ro(k),m0(k));
                    
                case 'nyomasszabalyzo'
                    switch szab{k}
                        case 'p'
                            for i=1:length(cpnevsor)
                                if strcmp(cpnevsor{i},szabcsp1{k}), szabcsp_szam=i; end
                            end
                            elemek{k}=nyomasszabalyzo(agnev{k},cspe_szam,cspv_szam,ro(k),m0(k),e(k,:),ek(k,:),szab,szabcsp_szam,ajel,pmax,szabP,szabI,szabD,szabcsp,e0,vmax,tbe);   
                            
                        case 'dp'
                            for i=1:length(cpnevsor)
                                if strcmp(cpnevsor{i},szabcsp1{k}), szabcsp_szam1=i; end
                                if strcmp(cpnevsor{i},szabcsp2{k}), szabcsp_szam2=i; end                        
                            end
                            elemek{k}=nyomasszabalyzo(agnev{k},cspe_szam,cspv_szam,ro(k),m0(k),e(k,:),ek(k,:),szab{k},szabcsp_szam1,szabcsp_szam2,ajel(k),pmax(k),szabP(k),szabI(k),szabD(k),szabcsp1{k},szabcsp2{k},e0(k),vmax(k),tbe(k));   
                    end    
                otherwise
                    error(['olvas5m -> Ismeretlen �gelem: ',tipus{k}])
            end
        end
        
        if debug>1, fprintf('   OK\n    merev alrendszer �p�t�se...'); end;
        clear csp;
        maxj=length(csom(1,:));
        for i=1:ncp_m
            j=1; clear agak;
            while j<=maxj & ~csom(i,j)==0
                %	fprintf('\n csomturka: i=%d  j=%d  csom(%d,%d)=%g',i,j,i,j,csom(i,j));
                agak(j)=csom(i,j); j=j+1;
            end      
            csp{i}{1}=i;
            csp{i}{2}=z(i);
            csp{i}{3}=agak;
            csp{i}{4}=1e5; % nyomas
            csp{i}{5}=elv(i);
            csp{i}{6}=cpnevsor{i};
        end
        
        mar{jc} = merev_alrendszer(1,cspnev,elemek,csp);
        
        % EDDIG!!!    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end   %...jc
    
    
