function fy = get_K(y,x,type)
% Nedves�tett ker�let meghat�roz�sa

switch type
    case 1
        fy=2*y+x;
    case 2
        r = x/2; theta=acos(1-y/r);
        fy=2*r*theta;
end