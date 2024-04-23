function info(buko,flag,varargin)

info(buko.tranziens_agelem_2csp,flag,varargin)

switch flag
    case 1
        str=info_string(buko);
        for i=1:length(str)
            fprintf('\n%s',str{i});
        end
    case 2
        fp=fopen(char(varargin{1}),'a'); 
        str=info_string(buko);
        for i=1:length(str)
            fprintf(fp,'\n%s',str{i});
        end
        fclose(fp);
        
    otherwise
        error('Ismeretlen opcio!');
end

function out=info_string(buko)    

out{1}=[' ',class(buko)];
out{2}=['   bukoszint:   ',buko.h0,' [m]'];
out{3}=['   atfolyasi tenyezo (Cd): ',buko.Cd,' [-]'];
out{4}=['   jellemzo meret     ',buko.B,' [m]'];
out{4}=['   bukoszint kitevo    ',buko.kitevo,' [-]'];
