function info(trag2,flag,varargin)

% flag=0 -> semmi
% flag=1 -> kepernyo
% flag=2 -> file, ami a varargin{1}-ben van
% flag=3 -> futas kozbeni informacio

switch flag
    case 1
        fprintf('\nElem: %s\n   %s:',trag2.nev,class(trag2));
        fprintf('\n\tcsomopont:       %d -> %d',trag2.csp(1),trag2.csp(2));
        fprintf('\n\tnyomasok:        p1=%g    p2=%g [Pa]',trag2.p(1),trag2.p(2));
        fprintf('\n\tterfogataram:    %g [m^3/s]',trag2.Q);
        fprintf('\n\ttomegaram:       %g [kg/s]',trag2.Q*trag2.ro);
        fprintf('\n\tfolyadek:        %s',trag2.folynev);
        fprintf('\n\t  suruseg:       %g [kg/m^3]',trag2.ro);
        fprintf('\n\t  rug.mod:       %g [Pa]',  trag2.B);
       
        
    case 2
        fp=fopen(char(varargin{1}),'a');
        fprintf(fp,'\nNev: %s\n   %s:',trag2.nev,class(trag2));
        fprintf(fp,'\n\tcsomopont:       %d -> %d',trag2.csp(1),trag2.csp(2));
        fprintf(fp,'\n\tnyomasok:        p1=%g    p2=%g [Pa]',trag2.p(1),trag2.p(2));
        fprintf(fp,'\n\tterfogataram:    %g [m^3/s]',trag2.Q);
        fprintf(fp,'\n\ttomegaram:       %g [kg/s]',trag2.Q*trag2.ro);
        fprintf(fp,'\n\tfolyadek:        %s',trag2.folynev);
        fprintf(fp,'\n\t  suruseg:       %g [kg/m^3]',trag2.ro);
        fprintf(fp,'\n\t  rug.mod:       %g [Pa]',  trag2.B);        
        fclose(fp);
        
    case 3
        t=varargin{1};
        fnev=char(strcat(trag2.tranziens_agelem.nev,'.res'));
        if fopen(fnev,'r')==-1
            fp = fopen(fnev,'w');
            fprintf(fp,'     t [s]             p1[Pa]            p2[Pa]            Q[m3/s]\n');
        else
            fp = fopen(fnev,'a'); 
        end   
        fprintf(fp,'   %5.3e       %+5.3e       %+5.3e       %+5.3e\n',t,trag2.pp(1),trag2.pp(2),trag2.tranziens_agelem.Q);
        fclose(fp);
        
    case 4
        fnev=char(strcat(trag2.tranziens_agelem.nev,'.res'));
        if fopen(fnev,'r')==-1
            error(['Nem tal�lom az eredm�nyfile-t!  -> ',fnev]);
        else
            fp = fopen(fnev,'r'); 
        end   
        temp=fgetl(fp);
        data=fscanf(fp,'%g %g %g %g',[4 inf]); data=data';
        fclose(fp);
        
        figure
        subplot(2,1,1),  plot(data(:,1),data(:,2),'+-',data(:,1),data(:,3),'o-');
        ylabel('p [Pa]'), grid on, title(['Eredm�nyek: ',trag2.tranziens_agelem.nev]);
        subplot(2,1,2),  plot(data(:,1),data(:,4),'+-');
        xlabel('t [s]'), ylabel('Q [m3/s]'), grid on
        
    otherwise
        fprintf('\n');
end