function info(vcssz,flag,varargin)

info(vcssz.tranziens_agelem_2csp,flag,varargin)

switch flag
    case 1
        fprintf('\n   %s:  ',class(vcssz));        
        
    case 2
        fp=fopen(char(varargin{1}),'a');
        fprintf(fp,'\n   %s:  ',class(vcssz));
	fclose(fp);

    otherwise
        error('Ismeretlen opcio!');
end
    
    