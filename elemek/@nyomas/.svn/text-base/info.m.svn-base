function info(nyomas,flag,varargin)

info(nyomas.tranziens_agelem_1csp,flag,varargin)

switch flag
    case 1
        fprintf('\n   %s:  ',class(nyomas));
        
    case 2
        fp=fopen(char(varargin{1}),'a');
        fprintf(fp,'\n   %s:  ',class(nyomas));
	fprintf(fp,'\n-----------------------------------------------------------------\n');
	fclose(fp);

    otherwise
        error('Ismeretlen opcio!');
end
    
    