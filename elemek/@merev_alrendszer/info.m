function info(mar,flag,varargin)

switch flag
    case 1
        for i=1:mar.n_elem
            temp=mar.elemek;
            info(temp{i},flag);
        end

    case 2
        for i=1:mar.n_elem
            temp=mar.elemek;
            info(temp{i},flag);
        end

    case 3
        t=varargin{1};
        workdir = mar.wdir;
        fnev=fullfile(workdir,strcat(mar.nev,'.res'));

        if fopen(fnev,'r')==-1
            fp = fopen(fnev,'w');
            %    fprintf(fp,'t[s],');
            %    for i=1:mar.n_elem, fprintf(fp,'Q%d,',i); end
            %    for i=1:mar.n_elem, fprintf(fp,'m%d,',i); end
            %    for i=1:mar.n_csp,  fprintf(fp,'p%d,',mar.csp{i}{1}); end
        else, fp = fopen(fnev,'a'); end

        fprintf(fp,'\n %8.6e ',t);
        for i=1:mar.n_elem
            temp=mar.elemek; fprintf(fp,' %+8.6e ',temp{i}.Q);
        end
        for i=1:mar.n_elem
            temp=mar.elemek;
            fprintf(fp,' %+8.6e ',temp{i}.Q*temp{i}.ro);
        end

        for i=1:mar.n_csp
            temp=mar.csp;
            fprintf(fp,' %+8.6e ',temp{i}{4});%-temp{i}{2}*1000*9.81);
        end

        for i=1:mar.n_elem
            temp=mar.elemek;
            if isa(temp{i},'szivattyu')
                fprintf(fp,' %+8.6e ',temp{i}.n);
            end
        end

        for i=1:mar.n_elem
            temp=mar.elemek;
            if isa(temp{i},'nyomasszabalyzo')
                fprintf(fp,' %+8.6e ',temp{i}.x);
                fprintf(fp,' %+8.6e ',temp{i}.e);
            end
        end
        
        for i=1:mar.n_elem
            temp=mar.elemek;
            if isa(temp{i},'akna')
                fprintf(fp,' %+8.6e ',temp{i}.y);
            end
        end   
        
        fclose(fp);

    case 5
        fnev=char(strcat(mar.nev,'.res'));
        data=load(fnev);
        plot_el  = varargin{1};
        plot_csp = varargin{2};
        temp=mar.elemek;

        figure
        for i=1:length(plot_el)
            subplot(length(plot_el),1,i)
            plot(data(:,1),data(:,plot_el(i)+1)); hold on
            legend(num2str(plot_el(i)));
            grid on, ylabel('Q [m3/s]')
            if i==1, title(['T�rfogat�ramok a "',mar.nev,'" alrendszer egyes elemeiben']); end
            if i==length(plot_el), xlabel('t [s]'); end
        end

        temp=mar.csp;
        nn=1+2*length(mar.elemek);
        figure
        for i=1:length(plot_csp)
            ize=1; while ~(temp{ize}{1}==plot_csp(i)), ize=ize+1; end
            subplot(length(plot_csp),1,i)
            plot(data(:,1),data(:,nn+ize)); hold on
            legend(num2str(plot_csp(i)));
            grid on, ylabel('p [Pa]')
            if i==1, title(['Nyom�slefut�sok a "',mar.nev,'" alrendszer egyes csom�pontjaiban']); end
            if i==length(plot_csp), xlabel('t [s]'); end
        end

    otherwise
        fprintf('\n');
end
