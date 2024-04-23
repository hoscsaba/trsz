function info(legust,flag,varargin)

info(legust.tranziens_agelem_1csp,flag,varargin);

switch flag
    case 1
        fprintf('\n   %s:  ',class(legust));
	fprintf('\n\tpolitrpokus kitevo:          %g [-]',legust.n);
        fprintf('\n\tkezdeti légtérfogat:         %g [m^3]',legust.V0);        
        fprintf('\n\tkezdeti légnyomás:           %g [bar]',legust.p0/1e5);
	fprintf('\n\tfelület:                     %g [m^2]',legust.A);
	fprintf('\n\tmagassag:                    %g [m]',legust.H);	
	fprintf('\n\ttalppont és bemenet közti magasságkülönbség: %g [m]\n\n',legust.l);
        
    case 2
        fp=fopen(char(varargin{1}),'a');
        fprintf(fp,'\n   %s:  ',class(legust));
	fprintf(fp,'\n\tpolitrpokus kitevo:          %g [-]',legust.n);
        fprintf(fp,'\n\tkezdeti légtérfogat:         %g [m^3]',legust.V0);        
        fprintf(fp,'\n\tkezdeti légnyomás:           %g [bar]',legust.p0/1e5);
	fprintf(fp,'\n\tfelület:                     %g [m^2]',legust.A);
	fprintf(fp,'\n\tmagassag:                    %g [m]',legust.H);		
	fprintf(fp,'\n\ttalppont és bemenet közti magasságkülönbség: %g [m]',legust.l);
	fprintf(fp,'\n-----------------------------------------------------------------\n');
	fclose(fp);

    otherwise
        error('Ismeretlen opcio!');
end
    
    