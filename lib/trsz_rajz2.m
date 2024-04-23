function trsz_rajz2(varargin)

  switch nargin
   case 5
    mdt  = varargin{1};
    rug  = varargin{2};
    hol1 = varargin{3};
    hol2 = varargin{4};
    mit  = varargin{5};
    van_df=0;
   case 7
    mdt  = varargin{1};
    rug  = varargin{2};
    hol1 = varargin{3};
    hol2 = varargin{4};
    mit  = varargin{5};
    dfnev = varargin{6};
    dfyy = varargin{7};
    fdata=load(dfnev);% fdata=fdata';
    dfx = fdata(:,1); dfy = fdata(:,dfyy);
    van_df=1;
   case 8
    mdt  = varargin{1};
    rug  = varargin{2};
    hol1 = varargin{3};
    hol2 = varargin{4};
    mit  = varargin{5};
    dfnev = varargin{6};
    dfyy = varargin{7};
    teltol = varargin{8};    
    fdata=load(dfnev);% fdata=fdata';
    dfx = fdata(:,1)-teltol; dfy = fdata(:,dfyy);
    van_df=1;
   case 10
    mdt  = varargin{1};
    rug  = varargin{2};
    hol1 = varargin{3};
    hol2 = varargin{4};
    mit  = varargin{5};
    dfnev = varargin{6};
    dfyy = varargin{7};
    teltol = varargin{8};    
    dataeltol = varargin{9};
    dataszorzo= varargin{10};
    fdata=load(dfnev);
    % fdata=fdata';
    dfx = fdata(:,1)+teltol; dfy = (fdata(:,dfyy)+dataeltol)*dataszorzo;
    van_df=1;
   otherwise
    error('balbal')
  end
  
debug=0;
restipus='nuku';

% Merev alrendszerr�l vagy rugalmas cs�r�l van sz�?

mar=olvas6m(mdt,debug);
[csovek,csp_nevsor,csom]=olvas7r(rug,debug);

n_sziv=0; n_nysz=0;
for i=1:length(mar)
    %  fprintf('\n Keresett merev alrendszer neve �s sz�ma: %s ?  %s, %d\n',hol1,mar{i}.nev,i);    
    if strcmp(mar{i}.nev,hol1)
        
        melyik1=i;
        n_ag=length(mar{i}.elemek);
        n_csp=length(mar{i}.csp);
        restipus='mar';
        % Meg kell sz�molni a szivattyukat es a nyomasszabalyzokat
        temp=mar{i}.elemek;
        for j=1:length(temp)   
            if isa(temp{j},'szivattyu'), n_sziv=n_sziv+1; end
            if isa(temp{j},'nyomasszabalyzo'), n_nysz=n_nysz+1; end
        end
%        fprintf('\n Keresett merev alrendszer neve �s sz�ma: %s, %d',mar{i}.nev,i);
%        fprintf('\n n_ag=%d  n_csp=%d  n_sziv=%d  n_nysz=%d\n',n_ag,n_csp,n_sziv,n_nysz);
    end
end

for i=1:length(csovek)
  if strcmp(csovek{i}.nev,hol1)
%    fprintf('\n Keresett rugalmas cs� neve �s sz�ma: %s, %d\n',csovek{i}.nev,i);
    melyik1=i;
    restipus='rug';
  end
end

if strcmp(restipus,'nuku')
  error(['Nincs ilyen merev alrendszer vagy rugalmas cso a rendszerben!! -> ',hol1]);
end

% Akkor most rajta!

switch restipus 
 case 'mar'   
  data=load(strcat(mar{melyik1}.nev,'.res'));
  switch mit
    % nyom�s csom�pontban
   case 'p'
     temp=mar{melyik1}.csp; 
     eltol=0;  
     for i=1:n_csp
       if strcmp(temp{i}{6},hol2)
	 eltol=i;	 
       end   
     end
     if eltol==0
       error(['Nincs ',hol2,' csomopont a ',mar{melyik1}.nev,' merev alrendszerben!'])  
     end

     if van_df==0
         figure
         plot(data(:,1),data(:,1+2*n_ag+eltol)/1e5);
         xlabel('t [s]'), ylabel('p [bar]'), grid on
         title(['Nyomas a ',mar{melyik1}.nev,' merev alrendszer ',temp{eltol}{6},' csomopontjaban']);
     else
         figure
         plot(data(:,1),data(:,1+2*n_ag+eltol)/1e5,dfx,dfy);
         xlabel('t [s]'), ylabel('p [bar]'), grid on
         title(['Nyomas a ',mar{melyik1}.nev,' merev alrendszer ',temp{eltol}{6},' csomopontjaban']);
         legend('Szimulacio',strcat('Adatfile: ',dfnev));
     end

     % t�rfogat�ram �gban
    case 'Q'
     eltol=0;
     temp=mar{melyik1}.elemek; 
     for i=1:n_ag
       if strcmp(temp{i}.nev,hol2)
	 eltol=i;	 
       end   
     end
     if eltol==0
       error(['Nincs ',hol2,' ág a ',mar{melyik1}.nev,' merev alrendszerben!'])  
     end
     
     if van_df==0
       figure
       plot(data(:,1),data(:,1+eltol));
       xlabel('t [s]'), ylabel('Q [m^3/s]'), grid on
       title(['Terfogataram a ',mar{melyik1}.nev,' merev alrendszer ',hol2,' agaban']);
     else
       figure
       plot(data(:,1),data(:,1+eltol),dfx,dfy);
       xlabel('t [s]'), ylabel('Q [m^3/s]'), grid on
       title(['Terfogataram a ',mar{melyik1}.nev,' merev alrendszer ',hol2,' agaban']); 
       legend('Szimulacio',strcat('Adatfile: ',dfnev));    
     end
     
     % t�meg�ram �gban
   case 'm'
    temp=mar{melyik1}.elemek; 
     eltol=0;  
     for i=1:n_ag
       if strcmp(temp{i}.nev,hol2)
	 eltol=i;	 
       end   
     end
     if eltol==0
       error(['Nincs ',hol2,' ag a ',mar{melyik1}.nev,' merev alrendszerben!'])  
     end

     if van_df==0
       figure
       plot(data(:,1),data(:,1+n_ag+eltol));
       xlabel('t [s]'), ylabel('m [kg/s]'), grid on
       title(['Tomegaram a ',mar{melyik1}.nev,' merev alrendszer ',hol2,' agaban']);
     else
       figure
       plot(data(:,1),data(:,1+n_ag+eltol),dfx,dfy);
       xlabel('t [s]'), ylabel('m [kg/s]'), grid on
       title(['Tomegaram a ',mar{melyik1}.nev,' merev alrendszer ',hol2,' agaban']); 
       legend('Szimulacio',strcat('Adatfile: ',dfnev));    
     end

     % szivatty� fordulatsz�m
    case 'n'
     eltol=0;  
     temp=mar{melyik1}.elemek;
     jj=0;
     for i=1:n_ag
       if isa(temp{i},'szivattyu')
	 jj=jj+1;
	 if strcmp(temp{i}.nev,hol2)
	   eltol=jj;
	 end
       end   
     end
     if eltol==0
       error(['Nincs ',hol2,' nevu szivattyu a ',mar{melyik1}.nev,' merev alrendszerben!'])  
     end

     % szivatty�
if van_df==0
     figure
     plot(data(:,1),data(:,1+2*n_ag+n_csp+eltol)*60);
     xlabel('t [s]'), ylabel('n [1/perc]'), grid on
     title(['A(z) ',hol2,' szivattyu fordulatszama (',mar{melyik1}.nev,' merev alrendszer) '])
     else
       figure
       plot(data(:,1),data(:,1+2*n_ag+n_csp+eltol)*60,dfx,dfy);
     xlabel('t [s]'), ylabel('n [1/perc]'), grid on
     title(['A(z) ',hol2,' szivattyu fordulatszama (',mar{melyik1}.nev,' merev alrendszer) '])
       legend('Szimulacio',strcat('Adatfile: ',dfnev));    
     end
     
    case 'x'
     eltol=0;  
     temp=mar{melyik1}.elemek;
     jj=0;
     for i=1:n_ag
       if isa(temp{i},'nyomasszabalyzo')
	 jj=jj+1;
	 if strcmp(temp{i}.nev,hol2)
	   eltol=jj;
	 end
       end   
     end
     if eltol==0
       error(['Nincs ',hol2,' nevű nyomásszabályzó a ',mar{melyik1}.nev,' merev alrendszerben!'])  
     end
     figure
     subplot(2,1,1)
     plot(data(:,1),data(:,1+2*n_ag+n_csp+n_sziv+2*eltol-1));
     ylabel('alapjel x [-]'), grid on
     title(['A(z) ',hol2,' nyomászabályzó jellemzői (',mar{melyik1}.nev,' merev alrendszer) '])
     subplot(2,1,2)
     plot(data(:,1),data(:,1+2*n_ag+n_csp+n_sziv+2*eltol));
     xlabel('t [s]'), ylabel('elmozduás e [-]'), grid on
     title(['A(z) ',hol2,' nyomászabályzó jellemzői (',mar{melyik1}.nev,' merev alrendszer) '])

    % vez_fojt jellegg�rbe
    case 'jg'
     eltol=0;  
     temp=mar{melyik1}.elemek;
     jj=0;
     for i=1:n_ag
%       class(temp{i})
       if isa(temp{i},'vez_fojtas')
	 jj=jj+1;
	 if strcmp(temp{i}.nev,hol2)
	   eltol=jj;
	   ez_az_elem=temp{i};
	 end
       end   
     end
     if eltol==0
       error(['Nincs ',hol2,' nev� vez�relt fojt�s a ',mar{melyik1}.nev,' merev alrendszerben!'])  
       end
       
     figure
     subplot(2,1,1)
     plot(ez_az_elem.epsK,ez_az_elem.K,'x-');
     xlabel('e/D [-]'), ylabel('K [-]'), grid on
     title(['A(z) ',hol2,' vez�relt fojt�s jellegg�rb�i (',mar{melyik1}.nev,' merev alrendszer) '])

     subplot(2,1,2)
     plot(ez_az_elem.t_vf,ez_az_elem.epst,'x-');
     ylabel('e/D [-]'), xlabel('t [s]'), grid on
     
   otherwise
    error('Ismeretlen opci� merev alrendszerre (p,Q,m,n)')
  
  end 
  
 case 'rug'
  data=load(strcat(csovek{melyik1}.nev,'.res'));
  switch mit
   case 'p'
    switch hol2
      case 'e'
       figure
       if van_df==0
	 plot(data(:,1),data(:,2)/1e5)
       else	 
	 plot(data(:,1),data(:,2)/1e5,dfx,dfy)
	 legend('Számítás',strcat('Adatfile: ',dfnev));
       end
       xlabel('t [s]'), ylabel('p [bar]'), grid on
       title(['Nyomás a ',csovek{melyik1}.nev,' cső elején']);
       
     case 'v'
       figure
       if van_df==0
	 plot(data(:,1),data(:,3)/1e5)
       else	 
	 plot(data(:,1),data(:,3)/1e5,dfx,dfy)
	 legend('Számítás',strcat('Adatfile: ',dfnev));
       end
       xlabel('t [s]'), ylabel('p [bar]'), grid on
       title(['Nyomás a ',csovek{melyik1}.nev,' cső végén']);
       
     case 'ev'
      figure
      if van_df==0
	plot(data(:,1),data(:,2)/1e5,data(:,1),data(:,3)/1e5)
	legend('eleje','vége');
      else
	plot(data(:,1),data(:,2)/1e5,data(:,1),data(:,3)/1e5,dfx,dfy,'r')
	legend('eleje','vége',strcat('Adatfile: ',dfnev));
      end
      xlabel('t [s]'), ylabel('p [bar]'), grid on
      title(['Nyomás a ',csovek{melyik1}.nev,' cső elején és végén']);
    end
    
   case 'v'
    switch hol2
      case 'e'
       figure      
       if van_df==0
	 plot(data(:,1),data(:,4));
       else	 
	 plot(data(:,1),data(:,4),dfx,dfy);
	 legend('Számítás',strcat('Adatfile: ',dfnev));
       end
       xlabel('t [s]'), ylabel('v [m/s]'), grid on
       title(['Sebesség a ',csovek{melyik1}.nev,' cső elején']);
       
     case 'v'
      figure
      if van_df==0
	plot(data(:,1),data(:,5));
      else	 
	plot(data(:,1),data(:,5),dfx,dfy)
	legend('Számítás',strcat('Adatfile: ',dfnev));
      end
      xlabel('t [s]'), ylabel('v [m/s]'), grid on
      title(['Sebesség a ',csovek{melyik1}.nev,' cső végén']);
      
     case 'ev'
      figure
      if van_df==0
	plot(data(:,1),data(:,4),data(:,1),data(:,5)), hold on       
	legend('eleje','vége')	
      else
	plot(data(:,1),data(:,4),data(:,1),data(:,5),dfx,dfy)
	legend('eleje','vége',strcat('Adatfile: ',dfnev));
      end      
      xlabel('t [s]'), ylabel('v [m/s]'), grid on
      title(['Sebesség a ',csovek{melyik1}.nev,' cső elején és végén']);end
   
   case 'Q'
    switch hol2
      case 'e'
       figure
       plot(data(:,1),data(:,6));
       xlabel('t [s]'), ylabel('Q [m^3/s]'), grid on
       title(['Térfogatáram a ',csovek{melyik1}.nev,' cső elején']);
     case 'v'
       figure
       plot(data(:,1),data(:,7));
       xlabel('t [s]'), ylabel('Q [m^3/s]'), grid on
       title(['Térfogatáram a ',csovek{melyik1}.nev,' cső végén']);
     case 'ev'
       figure
       plot(data(:,1),data(:,6),data(:,1),data(:,7));
       xlabel('t [s]'), ylabel('Q [m^3/s]'), grid on
       title(['Térfogatáram a ',csovek{melyik1}.nev,' cső elején és végén']);
       legend('eleje','v�ge')
    end
    
    case 'm'
    switch hol2
      case 'e'
       figure
       plot(data(:,1),data(:,8));
       xlabel('t [s]'), ylabel('m [kg/s]'), grid on
       title(['Tomegaram a ',csovek{melyik1}.nev,' cso elejen']);
     case 'v'
       figure
       plot(data(:,1),data(:,9));
       xlabel('t [s]'), ylabel('m [kg/s]'), grid on
       title(['Tomegaram a ',csovek{melyik1}.nev,' cso vegen']);
     case 'ev'
       figure
       plot(data(:,1),data(:,8),data(:,1),data(:,9));
       xlabel('t [s]'), ylabel('m [kg/s]'), grid on
       title(['Tomegaram a ',csovek{melyik1}.nev,' cso elejen es vegen']);
       legend('eleje','vege')
    end
   otherwise
    error('Ismeretlen opció rugalmas csőre (p,v,Q,m)')
  end
end

