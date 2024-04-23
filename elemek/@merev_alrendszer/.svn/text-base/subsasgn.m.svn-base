function mar = subsasgn(mar,index,val)
switch index.type
    case '.'
        switch index.subs
            case 'elemek'
                mar.elemek = val;
            case 'gorbek'
                mar.gorbek = val;
            case 'dtki'
                mar.dtki = val;
            case 'dtkiorig'
                mar.dtkiorig = val;
            otherwise
                error('mar: Invalid field name')
        end
end