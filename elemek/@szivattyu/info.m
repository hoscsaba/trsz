function info(sziv,flag,varargin)

% flag=0 -> semmi
% flag=1 -> kÃ©pernyÅ‘
% flag=2 -> file, ami a varargin{1}-ben van
% flag=3 -> futÃ¡si kÃ¶zbeni informÃ¡ciÃ³

info(sziv.tranziens_agelem_2csp,flag,varargin)

switch flag
    case 1
        fprintf('\n   %s:  ',class(sziv));        
        
    case 2
        fp=fopen(char(varargin{1}),'a');
        fprintf(fp,'\n   %s:  ',class(sziv));
        fprintf(fp,'\n\tQ [m^3/s] : ');
        for i=1:length(sziv.jgQ), fprintf(fp,' %7.4f ',sziv.jgQ(i)); end
        fprintf(fp,'\n\tH [m]     : ');
        for i=1:length(sziv.jgQ), fprintf(fp,' %7.1f ',sziv.jgH(i)); end
        if sziv.tranziens>0
	  fprintf(fp,'\n\tPbe [kW]  : ');
	  for i=1:length(sziv.jgQ), fprintf(fp,' %7.1f ',sziv.jgP(i)/1000); end
	end
	
%	fprintf('\n tranziens: %g\n',sziv.tranziens);
	
	switch sziv.tranziens
	 case 0
	  fprintf(fp,'\n\tSzimuláció típusa: stacioner');
	  fprintf(fp,'\n-----------------------------------------------------------------\n');
	  fclose(fp);
	 case 1
	  fprintf(fp,'\n\tSzimuláció típusa: kifutÃ¡s');
	  fprintf(fp,'\n\t   n    = %5.3e  [1/min]',sziv.n*60);
	  fprintf(fp,'\n\t   tki  = %5.3e  [s]',sziv.tki);
	  fprintf(fp,'\n\t   teta = %5.3e  [kgm^2]',sziv.teta);
	  fprintf(fp,'\n-----------------------------------------------------------------\n');
	  fclose(fp);
	 case 2
	  fprintf(fp,'\n\tSzimuláció típusa: direkt indÃ­tÃ¡s motor jelleggÃ¶rbÃ©vel');
	  fprintf(fp,'\n\t   nüz  = %5.3e  [1/min]',sziv.n);
	  fprintf(fp,'\n\t   teta = %5.3e  [kgm^2]',sziv.teta);
	  fprintf(fp,'\n\t   Mb   = %5.3e  [Nm]',sziv.Mb);
	  fprintf(fp,'\n\t   nb   = %5.3e  [1/min]',sziv.nb);
	  fprintf(fp,'\n\t   nsz  = %5.3e  [1/min]',sziv.nsz);
	  fprintf(fp,'\n-----------------------------------------------------------------\n');
	  fclose(fp);
	 case 4
	  fprintf(fp,'\n\tSzimuláció típusa: frekvenciaváltós indítás közelítö motor jelleggörbével');
	  fprintf(fp,'\n\t   nuz  = %5.3e  [1/min]',sziv.n);
	  fprintf(fp,'\n\t   teta = %5.3e  [kgm^2]',sziv.teta);
	  fprintf(fp,'\n\t   Mb   = %5.3e  [Nm]',sziv.Mb);
	  fprintf(fp,'\n\t   nb   = %5.3e  [1/min]',sziv.nb);
	  fprintf(fp,'\n\t   Mi   = %5.3e  [Nm]',sziv.Mi);
	  fprintf(fp,'\n\t   nsz  = %5.3e  [1/min]',sziv.nsz);
	  fprintf(fp,'\n\t   tbe  = %5.3e  [s]',sziv.tbe);	  
	  fprintf(fp,'\n\t   tind = %5.3e  [s]',sziv.tind);	  	  
	  fprintf(fp,'\n-----------------------------------------------------------------\n');
	  fclose(fp);
	end	
        
    case 5
        figure
        plot(sziv.jgQ,sziv.jgH,'ro',sziv.jgQ,sziv.jgH);
        xlabel('Q [m^3/s]'),ylabel('H [m]'), grid on
        % axis([0 1.1*max(sziv.jgQ) 0 1.2*max(sziv.jgH)])
        title(['A "',sziv.tranziens_agelem_2csp.nev,'" szivattyu jelleggorbeje'])    
    otherwise
        error('Ismeretlen opcio!');
    end
    
    