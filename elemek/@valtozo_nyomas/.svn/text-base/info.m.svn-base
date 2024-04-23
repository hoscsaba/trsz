function info(nyomas,flag,varargin)

info(nyomas.tranziens_agelem_1csp,flag,varargin)

switch flag
    case 1
        fprintf('\n   %s:  ',class(nyomas));
        
    case 2
        fp=fopen(char(varargin{1}),'a');
        fprintf(fp,'\n   %s:  ',class(nyomas));
	fprintf(fp,'\n jelleggörbe:');
	fprintf(fp,'\n\tt [s]   : ');  for i=1:length(nyomas.tt), fprintf(fp,' %5.3e ',nyomas.tt); end
        fprintf(fp,'\n\tp [bar] : ');  for i=1:length(nyomas.pp), fprintf(fp,' %5.3e ',nyomas.pp); end
	fprintf(fp,'\n-----------------------------------------------------------------\n');
	fclose(fp);

    otherwise
        error('Ismeretlen opcio!');
end
    
    