function [ Kp] = slanted_shell_stiffness_cisaillement_borne_1(n,E1,E3,v13,v31,rm,t,L,phi,r2)


syms  B D ;

coef=5/6;  %coeff de correction cisaillement
r=rm+r2*n;
%r=rm+L/2*n

    % Définition des composantes de la matrice B
    B(1,1)=-1/L;
    B(1,2)=0;
    B(1,3)=0;
    B(1,4)=1/L;
    B(1,5)=0;
    B(1,6)=0;
    
    
    B(2,1)=((1-n)*sin(phi))/(2*r);
    B(2,2)=((1-n)*cos(phi))/(2*r);
    B(2,3)=0;
    B(2,4)=((1+n)*sin(phi))/(2*r);
    B(2,5)=((1+n)*cos(phi))/(2*r);
    B(2,6)=0;
    
    B(3,1)=0;
    B(3,2)=0;
    B(3,3)=1/L;
    B(3,4)=0;
    B(3,5)=0;
    B(3,6)=-1/L;
    
    
    B(4,1)=0;
    B(4,2)=0;
    B(4,3)=-((1-n)*sin(phi))/(2*r);
    B(4,4)=0;
    B(4,5)=0;
    B(4,6)=-((1+n)*sin(phi))/(2*r);

    
    B(5,1)=0;
    B(5,2)=-1/L;
    B(5,3)=-(1-n)/2;
    B(5,4)=0;
    B(5,5)=1/L;
    B(5,6)=-(1+n)/2;
  
 %passage au repere global 
   %P1=[cos(phi),-sin(phi),0;sin(phi),cos(phi),0;0,0,-1];
   P1=[cos(phi),-sin(phi),0;sin(phi),cos(phi),0;0,0,-1];
   P2=zeros(3,3);
   Passage=[P1,P2;P2,P1];
   P=inv(Passage);

    %B1=Passage.*B.*inv(Passage);
    % Définition des composantes de la matrice D
    D(1,1)=E1*t/(1-v13*v31);
    D(1,2)=v13*E3*t/(1-v13*v31);
    D(1,3)=0;
    D(1,4)=0;
    D(1,5)=0;
    
    D(2,1)=v31*E1*t/(1-v13*v31);
    D(2,2)=E3*t/(1-v13*v31);
    D(2,3)=0;
    D(2,4)=0;
    D(2,5)=0;
    
    D(3,1)=0;
    D(3,2)=0;
    D(3,3)=E1*t^3/(12*(1-v13*v31));
    D(3,4)=v13*E3*t^3/(12*(1-v13*v31));
    D(3,5)=0;
    
    
    D(4,1)=0;
    D(4,2)=0;
    D(4,3)=v31*E1*t^3/(12*(1-v13*v31));
    D(4,4)=E3*t^3/(12*(1-v13*v31));
    D(4,5)=0;
    
    D(5,1)=0;
    D(5,2)=0;
    D(5,3)=0;
    D(5,4)=0;
    D(5,5)=(E1*coef*t)/(2*(1+v13));   %E1 et v13 a corig si composite

    % Matrice de raideur d'un élement coque 
    K=B'*D*B*r*L*pi();
   Kp=Passage*K*P;

   

end

