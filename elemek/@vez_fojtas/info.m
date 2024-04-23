function info(vfojtas,flag,varargin)

% flag=0 -> semmi
% flag=1 -> k�perny�
% flag=2 -> file, ami a varargin{1}-ben van
% flag=3 -> fut�s k�zbeni inform�ci�

info(vfojtas.tranziens_agelem_2csp,flag,varargin)

switch flag
    case 1
        fprintf('\n   %s...  ',class(vfojtas));
	fprintf('\n\t !! Kepernoyre info-t meg meg kell csinalni...\n\n');
        
    case 2
        fp=fopen(char(varargin{1}),'a');
        fprintf(fp,'\n   %s:',class(vfojtas));
        fprintf(fp,'\n\tjellemz� fel�let:    %g [m^2]',vfojtas.A);
        fprintf(fp,'\n\tfojt�si egy�tthat�:');
        fprintf(fp,'\n\t  e/d [-]: ');
        for i=1:length(vfojtas.epsK), fprintf(fp,'   %7.1f ',vfojtas.epsK(i)); end
        fprintf(fp,'\n\t  K   [-]: ');
        for i=1:length(vfojtas.K), fprintf(fp,' %5.3e ',vfojtas.K(i)); end
        fprintf(fp,'\n\tz�r�si f�ggv�nye:');
        fprintf(fp,'\n\t  t   [s]: ');
        for i=1:length(vfojtas.t_vf), fprintf(fp,' %7.1f ',vfojtas.t_vf(i)); end
        fprintf(fp,'\n\t  e/D [-]: ');
        for i=1:length(vfojtas.epst), fprintf(fp,' %7.1f ',vfojtas.epst(i)); end
        fprintf(fp,'\n-----------------------------------------------------------------\n');
        fclose(fp);
        
 case 5
  
  error('Nincs k�sz: @vez_fojtas/info')
        figure
        subplot(2,1,1), plot(fojtas.K_epD,fojtas.K,'-o')
        xlabel('e/D [-]'),ylabel('K [-]'), grid on
        title(['A "',fojtas.tranziens_agelem_2csp.nev,'" elem fojtasi t�nyez�je']);
        subplot(2,1,2), plot(fojtas.t,fojtas.t_epD,'-o')
        xlabel('t [s]'),ylabel('e/D [-]'), grid on
        title('Z�r�si f�ggv�ny');
        
    otherwise
        error('Ismeretlen opcio!');
    end
    