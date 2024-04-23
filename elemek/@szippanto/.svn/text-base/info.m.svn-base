function info(szip,flag,varargin)

info(szip.tranziens_agelem_1csp,flag,varargin);

switch flag
    case 1
        fprintf('\n   %s:  ',class(szip));
	fprintf('\n\tnyitónyomás:          %g [-]',szip.pnyit);
 
 case 2
        fp=fopen(char(varargin{1}),'a');
        fprintf(fp,'\n   %s:  ',class(szip));
	fprintf(fp,'\n\tnyitónyomás:          %g [-]',szip.pnyit);
	fprintf(fp,'\n-----------------------------------------------------------------\n');
	fclose(fp);

    otherwise
        error('Ismeretlen opcio!');
end
    
    