function [out,sziv]=setport(flag,sziv,varargin)

global g ro

g   = 9.81;
ro  = sziv.tranziens_agelem_2csp.ro;
csp = sziv.tranziens_agelem_2csp.csp;
cspe = csp(1);
cspv = csp(2);
nev = sziv.tranziens_agelem_2csp.nev;

if flag==-1 % inicializálás
    out{1}=1;
elseif flag==3
    % Adatok frissítése iteráció után
    % FONTOS!!!! Ez az xr TÉRFOGATáram!!!!
    xr = max(0,varargin{1});
    t  = varargin{2}(1);
    % dt = varargin{2}(2);
    
    % fprintf('\n %5.3e, %5.3e',varargin{1}/ro,sziv.Q);
            
    out{1}=1;
    sziv.n =sziv.nn;
    sziv.M =sziv.MM;
    sziv.Mm=sziv.MMm;
    
    if sziv.tranziens == 4
        % Frekvenciaváltós szivattyú indítás esetén lépésenkénti jelleggörbe
        % rajzolás.
        % Motor és szivattyú nyomaték jelleggörbe
        [nsz,nb] = motor_ford(sziv.nsz,sziv.nb,sziv.tbe,sziv.tind,t);
        
        p_n_sziv=linspace(0,sziv.nsz);
        p_n_mot =linspace(0,nsz*1.1);
        for i=1:length(p_n_sziv)
            p_Msz(i)= interp1(sziv.qq,sziv.mm,xr/sziv.n,'linear','extrap')*p_n_sziv(i)^2;
            p_Mm(i) =  motor_jg(nsz,nb,sziv.Mb,sziv.Mi,p_n_mot(i));
        end
        
        % Aktuális pontok
        p_Msza  = interp1(sziv.qq,sziv.mm,xr/sziv.n,'linear','extrap')*sziv.n^2;
        p_Mma   = motor_jg(nsz,nb,sziv.Mb,sziv.Mi,sziv.n);
        
        
        %         figure(1)
        %         clf
        %         plot(p_n_sziv*60,p_Msz/1000,'r-',p_n_mot*60,p_Mm/1000,'b-',sziv.n*60,p_Msza/1000,'r*',sziv.n*60,p_Mma/1000,'b*')
        %         xlabel('n [1/min]'), ylabel('M [kNm]'), grid on
        %         %hold on
        %         legend('szivattyu','motor','munkapont sziv.','munkapont motor',4)
        %         title(['t=',num2str(t),' s, n=',num2str(sziv.n*60),'/perc']);
        
        %fprintf('\n->       t=%5.3e t/tbe=%5.3e | Q=%+5.3e | n=%5.1f nsz=%5.1f nb=%5.1f | s=%+5.3f sb=%+5.3f Mm=%+5.2f Msz=%+5.2f | dn=%+5.2f nuj=%5.1f\n',t,t/sziv.tbe,xr/sziv.n,sziv.n*60,nsz*60,nb*60,1-sziv.n/nsz,1-nb/nsz,p_Mma,p_Msza,0,sziv.nn*60);
    end
else
    % iteracio
    % FONTOS!!!! Ez az xr TÖMEGáram!!!!
    xr = varargin{1};
    
    t  = varargin{2}{1}{1}(1);
    %    if strcmp(sziv.tranziens_agelem_2csp.nev,'edh1'), fprintf('\n\t\t setport:%7.5f',t); end
    dt = varargin{2}{1}{1}(2);
    mar = varargin{3};
    

    switch sziv.tranziens
        case 0 % stacioner számítás
            [a,b] = sz_eh(xr/ro,sziv.n,sziv.qq,sziv.a,sziv.b);
            out{1} = {-a/ro, 8/pi^2/9.81*(1/sziv.Dn^4-1/sziv.Ds^4)/ro^2, 0, -b, [csp(1),-1e5], [csp(2),1e5] };            
            
        case 1 % szivattyú kiesés
            if t>sziv.tki
                if xr/sziv.n<sziv.qq(1)
                    M=sziv.mm(1)*sziv.n^2;
                else
                    M=interp1(sziv.qq,sziv.mm,xr/sziv.n/ro)*sziv.n^2;
                end
                dn  = -dt*M/2/pi/sziv.teta;
                sziv.nn = max(sziv.n+dn,1/60);
                %if sziv.n<1/60, sziv.n=1/60; end
            else
                M=0; dn=0;
            end            

            [a,b] = sz_eh(xr/ro,sziv.n,sziv.qq,sziv.a,sziv.b);
            out{1} = {-a/ro, 8/pi^2/ro*(1/sziv.Dn^4-1/sziv.Ds^4), 0, -b, [csp(1),-1e5], [csp(2),1e5] };
            
        case 2 % direkt szivattyú indítás közelítő motor jelleggörbével
            if xr/ro/sziv.n<sziv.qq(1)
                M = sziv.mm(1)*sziv.n^2;
            elseif xr/ro/sziv.n>sziv.qq(length(sziv.qq))
                M = sziv.mm(length(sziv.mm))*sziv.n^2;
            else
                M = interp1(sziv.qq,sziv.mm,xr/sziv.n/ro)*sziv.n^2;
            end
            
            s   = 1-sziv.n/sziv.nsz;
            sb  = 1-sziv.nb/sziv.nsz;
            Mm  = 2*sziv.Mb/(s/sb+sb/s);
            dn  = dt/2/pi/sziv.teta*(Mm-M);
            sziv.nn = max(dn + sziv.nn,1/60);
            
            [a,b] = sz_eh(xr/ro,n,sziv.qq,sziv.a,sziv.b);
            out{1} = {-a/ro, 8/pi^2/9.81*(1/sziv.Dn^4-1/sziv.Ds^4)/ro^2, 0, -b, [csp(1),-1e5], [csp(2),1e5] };
            sziv.res=[sziv.res; t M Mm sziv.n];
            
        case 4
            % frekivaltos szivattyu inditas kozelito motor
            % jelleggorbevel
            
            % Szivattyu nyomatekigeny
            M = interp1(sziv.qq,sziv.mm,xr/sziv.n/ro,'linear','extrap')*sziv.n^2;
            
            % Motor nyomatek
            [nsz,nb] =motor_ford(sziv.nsz,sziv.nb,sziv.tbe,sziv.tind,t);
            Mm=motor_jg(nsz,nb,sziv.Mb,sziv.Mi,sziv.n);
            
            % Uj frodulatszam
            dn  = dt/2/pi/sziv.teta*(Mm-M);
            sziv.nn = dn + sziv.n;
            
            % Aktualis uzemi fordulatszam
            ni=linspace(0,nsz);
            for i=1:length(ni)
                Mmi(i) =  motor_jg(nsz,nb,sziv.Mb,sziv.Mi,ni(i));
                Msz(i) = interp1(sziv.qq,sziv.mm,xr/sziv.n/ro,'linear','extrap')*ni(i)^2;
            end
            nuz = interp1(Mmi-Msz,ni,0);
            
            %fprintf('\n\t t=%5.3e t/tbe=%5.3e | Q=%+5.3e | nr=%5.1f nsz=%5.1f nb=%5.1f | s=%+5.3f sb=%+5.3f Mm=%+5.1f Msz=%+5.1f | dn=%+5.2f nuj=%5.1f nuz=%5.1f',t,t/sziv.tbe,xr/sziv.n/ro,sziv.n*60,nsz*60,nb*60,1-sziv.nn/nsz,1-nb/nsz,Mm,M,dn*60,sziv.nn*60,nuz*60);
            
            % Ha a fordulatszam az uzemi fole szaladna, visszafogjuk.
            if sziv.nn>nuz
                %fprintf('\n\n\t A %s szivattyút hajtó motor fordulatszáma meghaladja az aktuális üzemi fordulatszámot!',nev);
                Mm  = interp1(ni,Mmi,nuz);
                M   = interp1(sziv.qq,sziv.mm,xr/sziv.n/ro,'linear','extrap')*nuz^2;
                %fprintf('\n\t Felulirom az aktualis uzemi fordulatszammal: t=%5.2f, n=%g/perc -> %g/perc (nb=%g/perc, nsz=%g/perc)\n',t,sziv.nn*60,nuz*60,nb*60,nsz*60);
                %fprintf('\n\t Ellenőrzés: ezen a forulatszámon Mm=%g Nm, Msz=%g Nm',Mm,M);
                sziv.nn=nuz;
            end
            
            sziv.MM=M; sziv.MMm=Mm;
            
            [a,b] = sz_eh(xr/ro,sziv.nn,sziv.qq,sziv.a,sziv.b);
            out{1} = {-a/ro, 8/pi^2/9.81*(1/sziv.Dn^4-1/sziv.Ds^4)/ro^2, 0, -b, [csp(1),-1e5], [csp(2),1e5] };
            
        case 5
            
            if sziv.init == 0
                for i = 1:length(mar.elemek)
                    if isa(mar.elemek{i},'akna')
                        csp_akna = mar.elemek{i}.csp;
                        if csp_akna(1) == cspe
                            hakt = mar.elemek{i}.y;
                            sziv.aknae = i;
                        end
                    end
                end
                sziv.init = 1;
            else
                hakt = mar.elemek{sziv.aknae}.y;
            end
           
            if sziv.uzem == 0
                sziv.nn = 1/60;
                if hakt >= sziv.hbe
                    sziv.uzem = 1;
                end
            elseif sziv.uzem == 1
                sziv.nn = sziv.n0;
                if hakt <= sziv.hki
                    sziv.uzem = 0;
                end
            end
            
            [a,b] = sz_eh(xr/ro,sziv.nn,sziv.qq,sziv.a,sziv.b);
            
            if sziv.uzem == 0 
                %a = 1e-6;
                %b = 1e-6;
                alfa=1e10;
                out{1} = {alfa/ro, 0, 0, 0, [csp(1),-1e5], [csp(2),1e5] };
            else
                out{1} = {-a/ro, 8/pi^2/ro*(1/sziv.Dn^4-1/sziv.Ds^4), 0, -b, [csp(1),-1e5], [csp(2),1e5] };
            end
            
            %fprintf('sziv: %s Q=%5.3e hakt: %8.3f a: %8.3f b: %8.3f sziv.uzem: %5d \n',nev,xr,hakt,a,b,sziv.uzem);

            case 6
            % frekivaltos szivattyu szab�lyoz�s nyom�sra
            % nyom�s megvizsg�l�sa:
              %pp = get(mar,'p',csp(2));
            %  deltap = (sziv.pny-pp); % vez�rl� jel
            % �j adat felv�tele
            %sziv.veg=sziv.veg+1;
            %sziv.adat(i) = deltap * dt;
            %sziv.ido(i) = t;
            %while (sziv.ido(sziv.kezd)<t-sziv.intdelay) 
            %        sziv.kezd = sziv.kezd+1;
            %end
            %Ide j�n az integr�l�s r�sze  !!!!!!!!!
            %szumma = 0;
            %for i=sziv.kezd:sziv.veg;
            %    szumma = szumma+sziv.adat(i);
            %end
            % integr�l n�vel�se
            if (t>0.01)
                
            if  (t>sziv.told)
               sziv.pp=get(mar,'pnev',sziv.pcsp); 
               sziv.deltap =(sziv.pny-sziv.pp);
               sziv.intszumma = sziv.intszumma + sziv.deltap*dt;           
               sziv.Hregi = sziv.Huj;
               sziv.told = t;
            end
            % Uj frodulatszam
            % elt�r�s becsl�se affinit�s alapj�n
            %Hregi = (pp - get(mar,'p',csp(1)))/ro/g;

            deltap = sziv.deltap; % vez�rl� jel
            
            if (t<1.5*dt) Hregi=interp1(sziv.qq,sziv.hh,xr/sziv.n/ro)*sziv.n^2; end;
            
            %puj  = sziv.pp + sziv.ParP*deltap + sziv.ParI*sziv.intszumma;
            sziv.Huj = sziv.Hregi + (sziv.ParP*deltap + sziv.ParI*sziv.intszumma)/ro/g;
            if sziv.Huj<0, sziv.Huj = 5; end;
            %sziv.Huj = (puj - get(mar,'p',csp(1)))/ro/g;
            %Quj = xr/ro;
%             if (Quj<1e-6) 
%                 nuj = 0;
%                 kk=0;
%                 Qhely = 0;
%             else
%                 kk = Huj/Quj^2;
%                 
%                 Qhely = interp1(sziv.kk,sziv.qq,kk,'linear','extrap');
%                 nuj=Quj/Qhely;
%             end;
%             %disp(sziv.qq);
            %disp(sziv.kk);
            
            
            nuj = (sziv.n^2*(sziv.Huj/sziv.Hregi))^0.5;
            
            if (nuj>sziv.nmax), nuj = sziv.nmax;  end
            if (nuj<0.5), nuj = 0.5; sziv.Huj = 5; end
            %fprintf('\nvez. szivatty� pakt=%5.2f [Pa]-> Hiba = %5.2f [Pa] -> puj = %5.2f [Pa]\n-> Qakt=%5.2e >> kk = %5.2e -> Qhely = %5.2e, nuj=%5.2f', sziv.pp, deltap, puj,Quj, kk, Qhely, 60*nuj); 
            %fprintf('\nvez. szivatty� t= %5.3f [s]: pakt=%5.2f [Pa]-> Hiba = %5.2f [Pa] ->\nHregi = %5.2e [m] Huj = %5.2e [m] -> nregi=%5.2f [1/min]   nuj=%5.2f [1/min]', t, sziv.pp, deltap, sziv.Hregi, sziv.Huj,60*sziv.n, 60*nuj); 
            sziv.nn = nuj;
            %pause;
            else 
               sziv.nn = sziv.n;
            end
            
            [a,b] = sz_eh(xr/ro,sziv.nn,sziv.qq,sziv.a,sziv.b);
            out{1} = {-a/ro, 8/pi^2/9.81*(1/sziv.Dn^4-1/sziv.Ds^4)/ro^2, 0, -b, [csp(1),-1e5], [csp(2),1e5] };
           
            
        otherwise
            error('Meg nincs kesz...');
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [a,b] = sz_eh(Q,n,qq,aa,bb)

global g ro

if Q/n<qq(1)
    b=bb(1)*n^2*ro*g;
    a=aa(1)*n*ro*g;
elseif Q/n>qq(length(qq))
    b=bb(length(qq))*n^2*ro*g;
    a=aa(length(qq))*n*ro*g;
else
    %b=interp1(qq,bb,Q/n,'linear')*n^2*ro*g;
    %a=interp1(qq,aa,Q/n,'linear')*n*ro*g;
    
    m1=find( qq<=Q/n ); m2=find( qq>Q/n );
    if Q/n > qq(end) || Q/n < qq(1)
        warning('extrap!')
    end
    b = ((bb(m2(1))-bb(m1(end)))*(Q/n-qq(m1(end)))/(qq(m2(1))-qq(m1(end)))+bb(m1(end)))*n^2*ro*g;
    a = ((aa(m2(1))-aa(m1(end)))*(Q/n-qq(m1(end)))/(qq(m2(1))-qq(m1(end)))+aa(m1(end)))*n*ro*g;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M = motor_jg(nsz,nb,Mb,Mi,n)
s   = 1- n/nsz;
sb  = 1-nb/nsz;

if n<nb
    M=interp1([0 nb*0.9 nb nsz],[Mi Mi*1.1 Mb 0],n,'cubic');
else
    M=2*Mb/(s^2+sb^2)*s*sb;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ansz,anb] = motor_ford(nsz,nb,tbe,tind,t)

mini=0.01;

if t<tbe
    ansz=nsz*mini;
    anb =nb*mini;
else
    ansz=nsz*min(1,mini+(1-mini)*(t-tbe)/tind);
    anb =nb *min(1,mini+(1-mini)*(t-tbe)/tind);
end

