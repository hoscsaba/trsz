function info(viszkcso,flag,varargin)

switch flag
    case 1
        info(viszkcso.tranziens_agelem_2csp,flag);
        fprintf('\n   %s:',class(viszkcso));
        fprintf('\n\tátmérõ:          %g [m]',viszkcso.D);
        fprintf('\n\tfalvastagság:    %g [m]',viszkcso.s);
        fprintf('\n\thossz:           %g [m]',viszkcso.L);
        fprintf('\n\tlambda:          %g\n',viszkcso.lambda);

    case 2
        info(viszkcso.tranziens_agelem_2csp,flag,varargin);

        fp=fopen(char(varargin{1}),'a');
        fprintf(fp,'\n   %s:',class(viszkcso));
        fprintf(fp,'\n\tátmérõ:          %g [m]',viszkcso.D);
        fprintf(fp,'\n\tfalvastagság:    %g [m]',viszkcso.s);
        fprintf(fp,'\n\thossz:           %g [m]',viszkcso.L);
        fprintf(fp,'\n\tlambda:          %g',viszkcso.lambda);
        fprintf(fp,'\n-----------------------------------------------------------------\n');
        fclose(fp);

    case 3
        %    t=varargin{1};
        workdir = varargin{1};
        temp = viszkcso.tranziens_agelem_2csp;
        %fnev = char(strcat(temp.nev,'.res'));
        fnev=fullfile(workdir,strcat(temp.nev,'.res'));
        if fopen(fnev,'r')==-1
            fp = fopen(fnev,'w');
            %        fprintf(fp,'t[s],p1[Pa],pN[Pa],v1[m/s],vN[m/s],m1[kg/s],mN[kg/s]\n');
        else
            fp = fopen(fnev,'a');
        end
        Q1 = viszkcso.v(1)*viszkcso.A(1);
        QN = viszkcso.v(viszkcso.N+1)*viszkcso.A(viszkcso.N+1);
        m1 = viszkcso.v(1)*viszkcso.A(1)*temp.ro;
        mN = viszkcso.v(viszkcso.N+1)*viszkcso.A(viszkcso.N+1)*temp.ro;
        %dseb1 = viszkcso.dseb(1);
        %dsebN = viszkcso.dseb(viszkcso.N+1);
        d1 = viszkcso.d(1);
        dN = viszkcso.d(viszkcso.N+1);
        epsz1 = viszkcso.epsz(1);
        epszN = viszkcso.epsz(viszkcso.N+1);
        epsz21 = viszkcso.epsz2(1);
        epsz2N = viszkcso.epsz2(viszkcso.N+1);
        a1 = viszkcso.a(1);
        aN = viszkcso.a(viszkcso.N+1);

        fprintf(fp,'%7.5e %+7.5e %+7.5e %+7.5e %+7.5e %+7.5e %+7.5e %+7.5e %+7.5e %+7.5e %+7.5e %+7.5e %+7.5e %+7.5e %+7.5e %+7.5e %+7.5e\n' ... 
                    ,viszkcso.t,viszkcso.p(1),viszkcso.p(viszkcso.N+1),viszkcso.v(1),viszkcso.v(viszkcso.N+1),Q1,QN,m1,mN,d1,dN,...
                    epsz1,epszN,epsz21,epsz2N,a1,aN);
        fclose(fp);

    case 5
        temp1=viszkcso.tranziens_agelem_2csp;
        %fnev=char(strcat(temp1.nev,'.res'));
        fnev=fullfile(workdir,strcat(temp.nev,'.res'));
        if fopen(fnev,'r')==-1
            error(['Nem találom az eredményfile-t!  -> ',fnev]);
        else
            fp = fopen(fnev,'r');
        end
        data=fscanf(fp,'%g %g %g %g %g %g %g %g %g',[9 inf]); data=data';
        fclose(fp);

        figure
        subplot(3,1,1),  plot(data(:,1),data(:,2),data(:,1),data(:,3));
        ylabel('p [Pa]'), legend('p_1','p_N'), grid on, title(['Eredmények: ',temp1.nev]);
        subplot(3,1,2),  plot(data(:,1),data(:,4),data(:,1),data(:,5));
        xlabel('t [s]'), ylabel('v [m/s]'), legend('v_1','v_N'), grid on
        subplot(3,1,3),  plot(data(:,1),data(:,8),data(:,1),data(:,9));
        xlabel('t [s]'), ylabel('m [kg/s]'), legend('m_1','m_N'), grid on

    case 4

        %  t=varargin{1};
        temp=viszkcso.tranziens_agelem_2csp;
        fnev=char(strcat(temp.nev,'_full.res'));
        if fopen(fnev,'r')==-1
            fp = fopen(fnev,'w');
            fprintf(fp,'t[s],p[Pa],v[m/s]\n');
        else
            fp = fopen(fnev,'a');
        end
        fprintf(fp,' %+5.3e ',viszkcso.t);
        for i=1:viszkcso.N+1, fprintf(fp,' %+5.3e ',viszkcso.p(i)); end
        for i=1:viszkcso.N+1, fprintf(fp,' %+5.3e ',viszkcso.v(i)); end
        fprintf(fp,'\n');
        fclose(fp);

    otherwise
        fprintf('\n');
end
