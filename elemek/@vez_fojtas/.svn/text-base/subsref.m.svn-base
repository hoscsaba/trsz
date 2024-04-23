function out = subsref(vfojt,index)

switch index.type
 case '.'
  switch index.subs
      case 'epsK',  out=vfojt.epsK;
      case 'K',     out=vfojt.K;    
      case 't_vf',  out=vfojt.t_vf;
      case 'epst',  out=vfojt.epst;
   otherwise
    out=subsref(vfojt.tranziens_agelem_2csp,index);      
  end      
 otherwise 
  error('HIBA!!! mezohivatkozas .-al (pl. elem.ro)')
end