function vena = set(vena,varargin)

property_argin = varargin;

while length(property_argin) >= 2,
  prop = property_argin{1};
  val = property_argin{2};
  property_argin = property_argin(3:end);
  
  switch prop
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Minden elemnel lennie kell 'fluid' opcionak.
    % Csak egyszeruen be kell masolni.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
   % Az elem folyadékbeállításai:        
     case 'fluidname'
        vena.fluidname=val;
     case 'ro'
        vena.ro=val;
     case 'nu'
        vena.nu=val;
     case 'mu'
        vena.mu=val;
    case 'B'
        vena.B=val;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	    
    
   case 'ido'
    vena.ido = val;    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   case 'hely'
    vena.hely = val;
    N=vena.N;
    dx=vena.dx;
    switch vena.hely
     
     case 'moc'
      vena.ido='moc';
      
     case 'moc'
      vena.ido='moc';
      
     case 'cds3'
      DD=zeros(N+1,N+1); 
      DD(1,1)=-2; DD(1,2)=2; DD(N+1,N)=-2; DD(N+1,N+1)=2;
      for i=2:N
	DD(i,i-1)=-1; DD(i,i+1)=1;
      end
      vena.DD=1/2/dx*DD;
      
     case 'cds5'
      DD=zeros(N+1,N+1); 
      DD(1,1)=-12;  DD(1,2)=12;  DD(2,1)=-6;    DD(2,3)=6; 
      DD(N,N-1)=-6; DD(N,N+1)=6; DD(N+1,N)=-12; DD(N+1,N+1)=12;
      for i=3:N-1
	DD(i,i-2)=1; DD(i,i-1)=-8; DD(i,i+2)=-1; DD(i,i+1)=8;
      end
      vena.DD=1/12/dx*DD;
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%
      %          FEM            %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%    
      
     case 'fem'
      AA=zeros(N+1,N+1); BB=zeros(N+1,N+1);
      AA(1,1)=2; AA(1,2)=1; AA(N+1,N)=1; AA(N+1,N+1)=2;
      for i=2:N
        AA(i,i-1)=1; AA(i,i)=4; AA(i,i+1)=1;
      end
      
      BB(1,1)=-1; BB(1,2)=1; BB(N+1,N)=-1; BB(N+1,N+1)=1;
      for i=2:N
        BB(i,i-1)=-1; BB(i,i+1)=1;
      end
      vena.DD=3/dx*inv(AA)*BB;
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%
      %    INTERLACING FEM      %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%    
      
     case 'femi'
      for i=1:N+1, x(i)=(i-1)*L/N; end
      dxx=L/(2*N-1);
      
      AA=zeros(N+1,N+1); 
      AA(1,1)=2; AA(1,2)=1; AA(2,1)=1; AA(2,2)=6; AA(2,3)=2; 
      AA(N+1,N)=2; AA(N+1,N+1)=4;
      for i=3:N
        AA(i,i-1)=2; AA(i,i)=8; AA(i,i+1)=2;
      end
      AA=AA/6*dxx;
      
      BB=zeros(N+1,N+1);
      BB(1,1)=-2; BB(1,2)=2; BB(2,1)=-5; BB(2,2)=4; BB(2,3)=1;
      BB(N,N-2)=-1; BB(N,N-1)=-5; BB(N,N)=4; BB(N,N+1)=2;
      BB(N+1,N-1)=-1; BB(N+1,N)=-5; BB(N+1,N+1)=6;
      for i=3:N-1
        BB(i,i-2)=-1; BB(i,i-1)=-5; BB(i,i)=5; BB(i,i+1)=1;
      end
      BB=BB/8;
      
      CC=zeros(N+1,N+1);
      CC(1,1)=4; CC(1,2)=2;
      CC(N,N-1)=2; CC(N,N)=6; CC(N,N+1)=1;
      CC(N+1,N)=1; CC(N+1,N+1)=2;
      for i=2:N-1
        CC(i,i-1)=2; CC(i,i)=8; CC(i,i+1)=2;
      end
      CC=CC/6*dxx;
      
      DDD=zeros(N+1,N+1);
      DDD(1,1)=-6; DDD(1,2)=5; DDD(1,3)=1;
      DDD(2,1)=-2; DDD(2,2)=-4; DDD(2,3)=5; DDD(2,4)=1;
      DDD(N,N-1)=-1; DDD(N,N)=-4; DDD(N,N+1)=5;
      DDD(N+1,N)=-2; DDD(N+1,N+1)=2;
      for i=3:N-1
        DDD(i,i-1)=-1; DDD(i,i)=-5; DDD(i,i+1)=5; DDD(i,i+2)=1;
      end
      DDD=DDD/8;

      EE=zeros(N+1,N+1);
      EE(1,1)=-1; EE(1,2)=1;
      EE(N+1,N)=-1; EE(N+1,N+1)=1;
      for i=2:N
        EE(i,i-1)=-1; EE(i,i+1)=1;
      end
      EE=EE*2;
      
      DD1=inv(AA)*BB; DD2=inv(CC)*DDD; DD3=inv(CC)*EE;
      vena.DD=[zeros(N+1,N+1) DD1; DD2 zeros(N+1,N+1)];
      
     otherwise
      fprintf('\n\n ERROR!! -> Unknown spatial disretization in %s',vena.name);
      fprintf('\n\t moc, cds3, cds5, fem, femi\n\n');
    end
   case 'p'
    for i=1:vena.N+1
      vena.p(i) = val;
    end
   case 'v'
    for i=1:vena.N+1
      vena.v(i) = val;
    end
   case 'D'
    vena.D = val;
   case 'nu'
    vena.nu = val;
   case 'adopt'
    vena.adopt = val;
   case 'reltol'
    vena.reltol = val;
   case 'friction'
    vena.friction = val;
    switch vena.friction
     case 'unsteady'                    
      switch vena.hely
       case 'moc'
	for i=1:vena.N+1
	  vena.vregi(i)=vena.v(i);
	  for j=1:vena.zielke_term
	    vena.I(i,j)=0;
	  end      
	end
      end
    end        
    
   case 'zielke_term'
    vena.zielke_term = val;
    for i=1:vena.N+1                
      for j=1:vena.zielke_term
	vena.I(i,j)=0;
      end      
    end
    
   otherwise
    error('Invalid pipe property');
    fprintf('\n\t adopt    - on,off');
    fprintf('\n\t ido      - beu,bdf2,esdirk3,esdirk4');
    fprintf('\n\t hely     - moc,cds3,cds5,fem,femi');
    fprintf('\n\t reltol   - <value>')
    fprintf('\n\t friction - none, steady, unsteady')            
  end
end
