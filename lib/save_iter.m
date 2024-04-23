function save_iter(fout,flag,varargin)

fp=fopen(fout,'a');

switch flag
    case 1
        fprintf(fp,'\n\n*******************************************************************************************************\n');
        fprintf(fp,'\n Lépés:  %d',varargin{3});
        fprintf(fp,'\n t  = %7.5e [s]',varargin{4}{1}(1));
        fprintf(fp,'\n dt = %7.5e [s]',varargin{4}{1}(2));
        fprintf(fp,'\n Iterációs hibahatárok:');
        fprintf(fp,'\n      rms max ágegy. = %5.3e [bar] ',varargin{1}(1));
        fprintf(fp,'\n      rms max konti. = %5.3e [kg/s]',varargin{1}(2));
        fprintf(fp,'\n      rms max rug.cs.= %5.3e [bar] ',varargin{1}(3));
        fprintf(fp,'\n Maximális lépésszám:  %g\n',varargin{2});
        fprintf(fp,'\n Az iteráció részletesen:\n\n');
        fprintf(fp,  '  i / imax |  rms ágegy. | rms konti.  | rms rug.cs. |');
        fprintf(fp,'\n-----------+-------------+-------------+-------------+');
    case 2
        fprintf(fp,'\n%3d / %3d  | %6.4e | %6.4e | %6.4e | %6.4e |',varargin{1},varargin{2},varargin{3}(1),varargin{3}(2),varargin{3}(3));
    case 3
        fprintf(fp,'\n-----------+-------------+-------------+-------------+');
        fprintf(fp,varargin{1});        
    otherwise
        error('Bajjj van...');
end
fclose(fp);