
clear all;
close all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% F. LACHAUD  
% Carac Mecanique
E=45000;
nu=0.3;
%Dimension
Lp=250;
b=20;
h=50;
S=b*h;
I=b*h*h*h/12;
EI=E*I;
% Rupture , définition de la fonction d'endommagement
EpsR=0.001;
Y0=0.5*EpsR*E*EpsR;
YC=25*Y0;

% Maillage
%
Nbe=3;
L=Lp/Nbe;

% Tableau de coordonnées des noeuds
noeuds=zeros(Nbe+1,1);
noeuds=linspace(0,Lp,Nbe+1)';



%% MATRICE de raideur de l'élément de poutre
%% Faire appel à la fonction en bas du prog
Ke=zeros(4,4,Nbe);
Ke=[12*EI/L^3 6*EI/L^2 -12*EI/L^3 6*EI/L^2 ;  6*EI/L^2 4*EI/L -6*EI/L^2 2*EI/L ; ...
    -12*EI/L^3 -6*EI/L^2 12*EI/L^3 -6*EI/L^2 ; 6*EI/L^2 2*EI/L -6*EI/L^2 4*EI/L ]
Kg=zeros(2*Nbe+2,2*Nbe+2);
Kr=zeros(2*Nbe,2*Nbe);

% Assemblage
for i=1:Nbe
    Kg(2*i-1:2*i+2,2*i-1:2*i+2) = Kg(2*i-1:2*i+2,2*i-1:2*i+2) + Ke(:,:);
end
% Matrice réduite conditions encastrée : suppression 2 ligne et deux
    Kr = Kg(3:end,3:end);

U=zeros(2*Nbe+2,1);
Ur=zeros(2*Nbe,1);
V=zeros(2*Nbe+2,1);
Vr=zeros(2*Nbe,1);
F=zeros(2*Nbe+2,1);
Fr=zeros(2*Nbe,1);
Res=zeros(2*Nbe,1);
Resg=zeros(2*Nbe+2,1);
EpsX=zeros(Nbe,1);
SigX=zeros(Nbe,1);
Energie=zeros(Nbe,1);
Energien=zeros(Nbe,1);
Energien(:,1)=(Y0);
d=zeros(Nbe,1);
dn=zeros(Nbe,1);
%Fn=zeros(4,4,Nbe)
%% Application du chargement extérieur
maxiter=15;
iter = 0;
te=0;
test=1.e-4; 
vmax=2.;            % si dep imposé dep max
Fmax=6000.;         % si effort imposé effort max

% Fonction de chargement
cycle=1;            % Nombre de cycle de charge-décharge
pas=20;             % Nombre de Pas par demi cycle (montée par exemple)
[pmax,lambda]=charge(cycle,pas);

% Boucle sur les pas de temps d'incrément de charge modele Effort imposé
for i=1:pmax+1
    iter = 0;
    
    % Effort imposé
    Fi=lambda(1,i)*Fmax;
    
    % Vecteur effort de F=K.q
    F(2*Nbe+1,1)=Fi;
    
    % Residu des efforts
    Resg=Kg*U-F;
    numpas=i;
    fprintf('pas numero %d \n',numpas);
    
    while (norm(Resg)>test && iter < maxiter) 
        iter = iter + 1;

         %Résidu reduit
         Res=Resg(3:end,1);
       
        
         %Matrice de raideur réduite
         Kr = Kg(3:end,3:end);
        
         % incrément de déplacement pour équilibre       
         deltaUr=-Kr\Res;
  
         % mise à jour des déplacements structure sur systeme reduit
         Ur=Ur+deltaUr;
  
         % Conditions aux limites encastrée-libres
         U(1,1)=0;
         U(2,1)=0;
         
         U(3:end,1)=Ur(1:end,1);        
         
        %==================================================================
        % Mise à jour des endommagements 
        % Calcul des dommages par éléments avec l'énergie
        %==================================================================
        % Epsilon=B*q
        % Energie=0.5*Eps*E*Eps
        % d=(Y-Y0)/(Yc-Y0)
        % E=E0*(1-d)
        % A FAIRE!
        % Fonction calcul des déformations par élément
        % Fonction pour calculer l'énergie
        % Fonction pour calculer l'endommagement
         Kg=zeros(2*Nbe+2,2*Nbe+2);
         
          for ne=1:Nbe
              
%              % Choix de la déformation moyenne par élément (Possible de
%              % faire le max... etc...              
              X1=noeuds(ne,1);
              X2=noeuds(ne+1,1);
              X=(X2-X1)/2;
%              
%              X=X1;
%              EpsX1=-0.5*h*(6/L^2*(-1+2*X/L)*U(ne,1)+1/L*(-4+6*X/L)*U(ne+1,1)...
%                      +6/L^2*(1-2*X/L)*U(ne+2,1)+1/L*(-2+6*X/L)*U(ne+4,1));
%              X=X2;
%              EpsX2=-0.5*h*(6/L^2*(-1+2*X/L)*U(ne,1)+1/L*(-4+6*X/L)*U(ne+1,1)...
%                      +6/L^2*(1-2*X/L)*U(ne+2,1)+1/L*(-2+6*X/L)*U(ne+4,1));
%              

%              EpsX(ne,1)=(EpsX1+EpsX2)/2;
                                    
%               X=L/2;
%               EpsX(ne,1)=0.5*h*(6/L^2*(-1+2*X/L)*U(ne,1)+1/L*(-4+6*X/L)*U(ne+1,1)...
%                      +6/L^2*(1-2*X/L)*U(ne+2,1)+1/L*(-2+6*X/L)*U(ne+4,1));
%              
%                Fi
%                  
%                EpsX
%                  
%              Energie(ne,1)=0.5*EpsX(ne,1)*E*EpsX(ne,1);
           
               Energie(ne,1)=0.5*EpsX(ne,1)*E*EpsX(ne,1);
               %Fi=Fr(2*Nbe-1,1);
               SigX(ne,1)=Fi*(Lp-X)/I*h/2;
               Energie(ne,1)=0.5*SigX(ne,1)*SigX(ne,1)/E;
               EpsX(ne,1)=SigX(ne,1)/E;
               Energie(ne,1)=0.5*EpsX(ne,1)*E*EpsX(ne,1);
             
             if (sqrt(Energie(ne,1)) > sqrt(Energien(ne,1)))
                 
                  d(ne,1)=(sqrt(Energie(ne,1))-sqrt(Y0))/(sqrt(YC)-sqrt(Y0));
                  
                 if (d(ne,1)>1) 
                      d(ne,1)=0.99;
                 end 
             else
                 
                  d(ne,1)=dn(ne,1);
             end
             
           
                          
             EI=E*I*(1-d(ne,1));
             Ke=[12*EI/L^3 6*EI/L^2 -12*EI/L^3 6*EI/L^2 ;  6*EI/L^2 4*EI/L -6*EI/L^2 2*EI/L ; ...
            -12*EI/L^3 -6*EI/L^2 12*EI/L^3 -6*EI/L^2 ; 6*EI/L^2 2*EI/L -6*EI/L^2 4*EI/L ];

             Kg(2*ne-1:2*ne+2,2*ne-1:2*ne+2) = Kg(2*ne-1:2*ne+2,2*ne-1:2*ne+2) + Ke(:,:);
             
          end

       
         % Mise à jour du Residu global en effort
          Resg=Kg*U-F;

          if (iter > maxiter)
              printf('attention dépassement des itérations : vérifier les résultats');
          end

    end
        dn=d;
        Energien=Energie; 


    
     
     % sauvegarde
    
    inc(i,1)=i;
    dep(i,1)=U(2*Nbe+1,1);
    theta(i,1)=U(2*Nbe+2,1)*180/pi();
    Effort(i,1)=max(F(2*Nbe+1,1));
    Epsilon1(i,1)=EpsX(1,1);
    Epsilonn(i,1)=EpsX(Nbe,1);
    Ener1(i,1)=(Energie(1,1));
    Enern(i,1)=(Energie(Nbe,1));
    module(i,1)=E*(1-d(1,1));
    endom(i,1)=d(1,1);
end 

    figure(1);
    plot(inc,lambda);
    figure(2);
    plot(inc,Effort);

    figure(3);
    plot(dep,Effort,'--go',theta,Effort,'-c*');
    abs=xlabel('Angle (°),   Deplacement (mm)');
    ord=ylabel('Effort');
    
    figure(4);
    plot(inc,Ener1,inc,Enern);
    abs=xlabel('Inc');
    ord=ylabel('Energie');
    
    figure(5);
    plot(inc,Epsilon1,inc,Epsilonn);
    abs=xlabel('increment');
    ord=ylabel('Epsilon');
    
    figure(6);
    plot(inc,endom);
    abs=xlabel('increment');
    ord=ylabel('Endommagement');
    
    figure(7);
    plot(Ener1,Epsilon1,Enern,Epsilonn);
    abs=xlabel('Energie');
    ord=ylabel('Deformation');
    
    
% Valeur theorique    
   
    Fmm=max(Effort);
    Fmm
    vmm=Fmax*Lp^3/(3*EI)

% Valeur calculée
    vc=max(dep)
    
    
function [pmax,lambda]=charge(cycle,pas);
         T=1/cycle;  
         Te=T/pas;
         t1=0:Te:(T-Te);t2=Te:Te:(T-Te);
         triang=[t1 T fliplr(t2)]'; % Contient la fonction triangle
         croiss=[1:cycle];
         lambda(1,:)=reshape(triang*croiss,1,[]); % Contient le cycle
         pmax=length(lambda)+1; % Nombre de pas de temps total
         lambda(1,pmax+1)=0.;
end


% Fonction de calcul de l'énergie max dans chaque élément
function [Energie,EpsX] = energie_beam(Le,U)

end

% Exemple de fonction pour le calcul de Ke et K a faire
function [Ke,Kg] = raideur_beam(Le,EI)
if ~(length(kb)==length(l))
    error('Base stiffness and length arrays must have equal number of elements.');
end

num_elem = length(kb);
num_node = num_elem + 1;

ke = zeros(4,4,length(kb));
for cnt=1:num_elem
    ke(:,:,cnt) = kb(cnt) * [12 6*l(cnt) -12 6*l(cnt);
        6*l(cnt) 4*l(cnt)^2 -6*l(cnt) 2*l(cnt)^2;
        -12 -6*l(cnt) 12 -6*l(cnt);
        6*l(cnt) 2*l(cnt)^2 -6*l(cnt) 4*l(cnt)^2];
end

K = zeros(2*num_node,2*num_node);
for cnt=1:num_elem
    K(2*cnt-1:2*cnt+2,2*cnt-1:2*cnt+2) = K(2*cnt-1:2*cnt+2,2*cnt-1:2*cnt+2) + ...
        ke(:,:,cnt);
end

 Kg=K;
 Ke=ke;
end
