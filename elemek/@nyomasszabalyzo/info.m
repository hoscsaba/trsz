function info(nysz,flag,varargin)

info(nysz.tranziens_agelem_2csp,flag,varargin)

switch flag
    case 1
        fprintf('\n   %s:  ',class(nysz));
        fprintf('\n\tjelleggörbe:')
        fprintf('\n\te [-] : ');
        for i=1:length(nysz.jge), fprintf(' %7.4e ',nysz.jge(i)); end
        fprintf('\n\tK [-] : ');
        for i=1:length(nysz.jgk), fprintf(' %7.4e ',nysz.jgk(i)); end
        
        switch nysz.szab
            case 'p'
                fprintf('\n\tszab. tipus:     nyomás');
                fprintf('\n\tszab.csp neve:   %s',nysz.szcspnev1);
                fprintf('\n\tszab.csp szama:  %d',nysz.szcsp1);
            case 'dp'                     
                fprintf('\n\tszab. tipus:     nyomáskülönbség');
                fprintf('\n\tszab.cspok neve:   %s - %s',nysz.szcspnev1,nysz.szcspnev2);
                fprintf('\n\tszab.cspok szama:  %d - %d',nysz.szcsp1,nysz.szcsp2);
        end
        fprintf('\n\talapjel:         %g [bar]',nysz.ajel);
        fprintf('\n\tpmax:            %g [bar]',nysz.pmax);
        fprintf('\n\tarányos tag:     %g',nysz.szabP);
        fprintf('\n\tint. tag:        %g',nysz.szabI);
        fprintf('\n\tdiff. tag:       %g',nysz.szabD);
        fprintf('\n\tkezdeti állás:   %g',nysz.e);
        fprintf('\n\tbekapcsolás:     %g [s]',nysz.tbe);
        
    case 2
        fp=fopen(char(varargin{1}),'a');
        fprintf(fp,'\n   %s:  ',class(nysz));
        fprintf(fp,'\n\tjelleggörbe:');
        fprintf(fp,'\n\te [-] : ');
        for i=1:length(nysz.jge), fprintf(fp,' %7.4e ',nysz.jge(i)); end
        fprintf(fp,'\n\tK [-] : ');
        for i=1:length(nysz.jgk), fprintf(fp,' %7.4e ',nysz.jgk(i)); end
        switch nysz.szab
            case 'p'
                fprintf(fp,'\n\tszab. tipus:     nyomás');
                fprintf(fp,'\n\tszab.csp neve:   %s',nysz.szcspnev1);
                fprintf(fp,'\n\tszab.csp szama:  %d',nysz.szcsp1);
            case 'dp'                     
                fprintf(fp,'\n\tszab. tipus:     nyomáskülönbség');
                fprintf(fp,'\n\tszab.cspok neve:   %s - %s',nysz.szcspnev1,nysz.szcspnev2);
                fprintf(fp,'\n\tszab.cspok szama:  %d - %d',nysz.szcsp1,nysz.szcsp2);
        end
        fprintf(fp,'\n\talapjel:         %g [bar]',nysz.ajel);
        fprintf(fp,'\n\tpmax:            %g [bar]',nysz.pmax);
        fprintf(fp,'\n\tarányos tag:     %g',nysz.szabP);
        fprintf(fp,'\n\tint. tag:        %g',nysz.szabI);
        fprintf(fp,'\n\tdiff. tag:       %g',nysz.szabD);
        fprintf(fp,'\n\tkezdeti állás:   %g',nysz.e);
        fprintf(fp,'\n\tbekapcsolás:     %g [s]',nysz.tbe);
        
        fprintf(fp,'\n-----------------------------------------------------------------\n');
        fclose(fp);
        
    otherwise
        error('Ismeretlen opcio!');
    end
    
    