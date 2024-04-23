function info(akna,flag,varargin)

info(akna.tranziens_agelem_1csp,flag,varargin);

switch flag
    case 1
        str=info_string(akna);
        for i=1:length(str)
            fprintf('\n%s',str{i});
        end
        fprintf('\n\n');
    case 2
        fp=fopen(char(varargin{1}),'a'); 
        str=info_string(akna);
        for i=1:length(str)
            fprintf(fp,'\n%s',str{i});
        end
        fprintf(fp,'\n\n');
        fclose(fp);
        
    otherwise
        error('Ismeretlen opcio!');
end

function out=info_string(akna) 

out{1}=[' ',class(akna)];
out{2}=['    alapterulet:   ',num2str(akna.A),' [m^2]'];
out{3}=['    kezdeti szint: ' ,num2str(akna.y0-akna.hmin),' [m]'];
out{4}=['    fenekszint:    ',num2str(akna.hmin),' [m]'];
out{5}=['    fedlapszint:   ',num2str(akna.hmax),' [m]'];