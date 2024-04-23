function info(cso,flag,varargin)

switch flag
    case 1
        info(cso.tranziens_agelem_2csp,flag);
        fprintf('\n   %s:',class(cso));
        fprintf('\n\tatmero:          %g [m]',cso.D);
        fprintf('\n\tfalvastagsag:    %g [m]',cso.s);
        fprintf('\n\thossz:           %g [m]',cso.L);
        fprintf('\n\tlambda:          %g\n',cso.lambda);
        fprintf('\n\t  eleje pf. tipus:    %s',cso.pf_eleje_tipus);
        fprintf('\n\t  vege pf. tipus :    %s\n',  cso.pf_vege_tipus);
    case 2
        info(cso.tranziens_agelem_2csp,flag,varargin);
        
        fp=fopen(char(varargin{1}),'a');
        fprintf(fp,'\n   %s:',class(cso));
        fprintf(fp,'\n\tatmero:          %g [m]',cso.D);
        fprintf(fp,'\n\tfalvastagsag:    %g [m]',cso.s);
        fprintf(fp,'\n\thossz:           %g [m]',cso.L);
        fprintf(fp,'\n\tlambda:          %g',cso.lambda);
        fprintf(fp,'\n-----------------------------------------------------------------\n');
        fclose(fp);
        fprintf('\n\t  eleje pf. tipus:    %s',cso.pf_eleje_tipus);
        fprintf('\n\t  vege pf. tipus :    %s\n',  cso.pf_vege_tipus);
    case 3
        %    t=varargin{1};
        temp=cso.tranziens_agelem_2csp;
        fnev=char(strcat(temp.nev,'.res'));
        if fopen(fnev,'r')==-1
            fp = fopen(fnev,'w');
            %        fprintf(fp,'t[s],p1[Pa],pN[Pa],v1[m/s],vN[m/s],m1[kg/s],mN[kg/s]\n');
        else
            fp = fopen(fnev,'a');
        end
        Q1=cso.v(1)*cso.A;
        QN=cso.v(cso.N+1)*cso.A;
        m1=cso.v(1)*cso.A*temp.ro;
        mN=cso.v(cso.N+1)*cso.A*temp.ro;
        fprintf(fp,'%7.5e  %+7.5e  %+7.5e  %+7.5e  %+7.5e  %+7.5e  %+7.5e  %+7.5e  %+7.5e\n',cso.t,cso.p(1),cso.p(cso.N+1),cso.v(1),cso.v(cso.N+1),Q1,QN,m1,mN);
        fclose(fp);

    case 5
        temp1=cso.tranziens_agelem_2csp;
        fnev=char(strcat(temp1.nev,'.res'));
        if fopen(fnev,'r')==-1
            error(['Nem tal�lom az eredm�nyfile-t!  -> ',fnev]);
        else
            fp = fopen(fnev,'r');
        end
        data=fscanf(fp,'%g %g %g %g %g %g %g %g %g',[9 inf]); data=data';
        fclose(fp);

        figure
        subplot(3,1,1),  plot(data(:,1),data(:,2),data(:,1),data(:,3));
        ylabel('p [Pa]'), legend('p_1','p_N'), grid on, title(['Eredm�nyek: ',temp1.nev]);
        subplot(3,1,2),  plot(data(:,1),data(:,4),data(:,1),data(:,5));
        xlabel('t [s]'), ylabel('v [m/s]'), legend('v_1','v_N'), grid on
        subplot(3,1,3),  plot(data(:,1),data(:,8),data(:,1),data(:,9));
        xlabel('t [s]'), ylabel('m [kg/s]'), legend('m_1','m_N'), grid on

    case 4

        %  t=varargin{1};
        temp=cso.tranziens_agelem_2csp;
        fnev=char(strcat(temp.nev,'_full.res'));
        if fopen(fnev,'r')==-1
            fp = fopen(fnev,'w');
            fprintf(fp,'t[s],p[Pa],v[m/s]\n');
        else
            fp = fopen(fnev,'a');
        end
        fprintf(fp,' %+5.3e ',cso.t);
        for i=1:cso.N+1, fprintf(fp,' %+5.3e ',cso.p(i)); end
        for i=1:cso.N+1, fprintf(fp,' %+5.3e ',cso.v(i)); end
        fprintf(fp,'\n');
        fclose(fp);

    otherwise
        fprintf('\n');
end
