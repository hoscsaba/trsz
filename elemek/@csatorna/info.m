function info(csatorna,flag,varargin)

switch flag
    case 1
        info(csatorna.tranziens_agelem_2csp,flag);
        fprintf('\n   %s:',class(csatorna));
        fprintf('\n\tdvB:            %g [m]',csatorna.dvB);
        fprintf('\n\thossz:           %g [m]',csatorna.L);
        fprintf('\n\tn:          %g\n',csatorna.n);
        fprintf('\n\t  eleje pf. tipus:    %s',trag2.pf_eleje_tipus);
        fprintf('\n\t  vege pf. tipus :    %s\n',  trag2.pf_vege_tipus);
    case 2
        info(csatorna.tranziens_agelem_2csp,flag,varargin);

        fp=fopen(char(varargin{1}),'a');
        fprintf(fp,'\n   %s:',class(csatorna));
        fprintf('\n\tdvB:            %g [m]',csatorna.dvB);
        fprintf('\n\thossz:           %g [m]',csatorna.L);
        fprintf('\n\tn:          %g\n',csatorna.n);
        fprintf(fp,'\n-----------------------------------------------------------------\n');
        fclose(fp);
        fprintf('\n\t  eleje pf. tipus:    %s',trag2.pf_eleje_tipus);
        fprintf('\n\t  vege pf. tipus :    %s\n',  trag2.pf_vege_tipus);
    case 3
        %    t=varargin{1};
        workdir = varargin{1};

        temp=csatorna.tranziens_agelem_2csp;
        %fnev=char(strcat(temp.nev,'.res'));
        fnev=fullfile(workdir,strcat(temp.nev,'.res'));
        %fnev=char(strcat(workdir,'\',temp.nev,'.res'));

        Ae= get_A(csatorna,csatorna.y(1));
        %Av= get_A(csatorna,csatorna.y(end));
        Av= get_A(csatorna,csatorna.y(end-1));
        Q1=csatorna.v(1)*Ae;
        %QN=csatorna.v(end)*Av;
        QN=csatorna.v(end-1)*Av;
        m1=csatorna.v(1)*Ae*temp.ro;
        %mN=csatorna.v(end)*Av*temp.ro;
        mN=csatorna.v(end-1)*Av*temp.ro;

        %if csatorna.kiir == 1
            if fopen(fnev,'r')==-1
                fp = fopen(fnev,'w');
                %        fprintf(fp,'t[s],p1[Pa],pN[Pa],v1[m/s],vN[m/s],m1[kg/s],mN[kg/s]\n');
            else
                fp = fopen(fnev,'a');
            end

            fprintf(fp,'\n %7.5e %+7.5e  %+7.5e  %+7.5e  %+7.5e  %+7.5e  %+7.5e  %+7.5e  %+7.5e',csatorna.t,csatorna.y(1),csatorna.y(end),csatorna.v(1),csatorna.v(end),Q1,QN,m1,mN);
            for i = 1:(csatorna.N+1)
                fprintf(fp,'  %+7.5e',csatorna.y(i));
            end

            fprintf(fp,'  %+7.5e',csatorna.dt);
            fclose(fp);


    case 5
        temp1=csatorna.tranziens_agelem_2csp;
        fnev=char(strcat(temp1.nev,'.res'));
        if fopen(fnev,'r')==-1
            error(['Nem talalom az eredmenyfile-t!  -> ',fnev]);
        else
            fp = fopen(fnev,'r');
        end
        data=fscanf(fp,'%g %g %g %g %g %g %g %g %g',[9 inf]); data=data';
        fclose(fp);

        figure
        subplot(3,1,1),  plot(data(:,1),data(:,2),data(:,1),data(:,3));
        ylabel('p [Pa]'), legend('p_1','p_N'), grid on, title(['Eredmenyek: ',temp1.nev]);
        subplot(3,1,2),  plot(data(:,1),data(:,4),data(:,1),data(:,5));
        xlabel('t [s]'), ylabel('v [m/s]'), legend('v_1','v_N'), grid on
        subplot(3,1,3),  plot(data(:,1),data(:,8),data(:,1),data(:,9));
        xlabel('t [s]'), ylabel('m [kg/s]'), legend('m_1','m_N'), grid on

    case 4

        %  t=varargin{1};
        temp=csatorna.tranziens_agelem_2csp;
        fnev=char(strcat(temp.nev,'_full.res'));
        if fopen(fnev,'r')==-1
            fp = fopen(fnev,'w');
            fprintf(fp,'t[s],p[Pa],v[m/s]\n');
        else
            fp = fopen(fnev,'a');
        end
        fprintf(fp,' %+5.3e ',csatorna.t);
        for i=1:csatorna.N+1, fprintf(fp,' %+5.3e ',csatorna.p(i)); end
        for i=1:csatorna.N+1, fprintf(fp,' %+5.3e ',csatorna.v(i)); end
        fprintf(fp,'\n');
        fclose(fp);

    otherwise
        fprintf('\n');
end
