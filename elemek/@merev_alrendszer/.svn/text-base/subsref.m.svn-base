function out = subsref(elem,index)

switch index(1).type
    case '.'            
        switch index(1).subs
            case 'elemek'
                if length(index) == 1
                    out = elem.elemek;
                else
                    switch index(2).type
                        case '{}'
                            if length(index) == 2
                                out = elem.elemek{index(2).subs{1}};
                            elseif length(index) == 3

                                out = subsref( elem.elemek{index(2).subs{1}}, index(3) );
                            end
                    end
                end
            case 'csp'
                out = elem.csp;
            case 'plot_it'
                out = elem.plot_it;
            case 'nev'
                out = elem.nev;  
            case 'x'
                out = elem.x;
	        case 't'
                out = elem.t;
            case 'dtki'
                out = elem.dtki;
            otherwise
                error('Ismeretlen mezo: ',index.type);
        end
    otherwise
        error('HIBA!!! mezohivatkozas .-al!')
end