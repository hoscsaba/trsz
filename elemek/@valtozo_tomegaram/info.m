function info(tomegaram,flag,varargin)

info(tomegaram.tranziens_agelem_1csp,flag,varargin)

switch flag
    case 1
        fprintf('\n   %s:  ',class(tomegaram));
        
    case 2
        fp=fopen(char(varargin{1}),'a');
        fprintf(fp,'\n   %s:  ',class(tomegaram));
	fprintf(fp,'\n jelleggörbe:');
	fprintf(fp,'\n\tt [s]   : ');  for i=1:length(tomegaram.tt), fprintf(fp,' %5.3e ',tomegaram.tt); end
    %fprintf(fp,'\n\tp [bar] : ');  for i=1:length(tomegaram.pp), fprintf(fp,' %5.3e ',tomegaram.pp); end
	fprintf(fp,'\n-----------------------------------------------------------------\n');
	fclose(fp);

    otherwise
        error('Ismeretlen opcio!');
end
    
    