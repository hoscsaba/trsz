function out = subsref(nysz,index)

switch index.type
 case '.'
  switch index.subs
      case 'xx', out=nysz.xx;
      case 'x',  out=nysz.xx(4);    
      case 'e',  out=nysz.e;              
   otherwise
    out=subsref(nysz.tranziens_agelem_2csp,index);      
  end      
 otherwise 
  error('HIBA!!! mezohivatkozas .-al (pl. elem.ro)')
end