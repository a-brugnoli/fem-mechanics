clear all;
close all;
clc;
%% Construction de la geometrie et maillage

%% %%%%%%% %%
%% Données %%
%% %%%%%%% %%
Fmax=480664;
E1=210000;
E2=210000;
v12=0.33;
v21=v12*E2/E1;
r_int=1020;
r_ext=1064;
rr=4; %rayon de raccordement
L=40;
h=296; %300;
t=5;
N=20;
amp=100;
angle_virole=90*pi/180;

taille_elt_pied_bride=2;
taille_elt_virole=2;
nbr_element_raccr=6;
angle_sect_maillage_raccr=(90/nbr_element_raccr)*pi/180;
nbr_elem_pied_bride=round(L/taille_elt_pied_bride);
nbr_elem_virole=round(h/taille_elt_virole);




%% %%%%%%%% %%
%% Maillage %%
%% %%%%%%%% %%

noeud=[];
w=0;
elem=[];
noeud(1).r=r_int;
noeud(1).z=300;

%Maillage virole
for k=2:1:nbr_elem_virole+1;
    noeud(k).r=r_int;
    noeud(k).z=noeud(1).z-sin(angle_virole)*(h/nbr_elem_virole*(k-1)) ; 
    elem(k-1).noeud1=k-1;
    elem(k-1).noeud2=k;
    elem(k-1).E1=E1;
    elem(k-1).E2=E2;
    elem(k-1).nu12=v12;
    elem(k-1).type=1;% coque cylindrique
    elem(k-1).epaiss=t; % epaisseur reel
    elem(k-1).longueur=h/nbr_elem_virole;
    elem(k-1).rd=noeud(k).r; % pour simplifier le repérage de rtube dans le cas de virole à 90°
    elem(k-1).rg=noeud(k-1).r;
    elem(k-1).angle=0;
end
d=1;
angle=180-[168.749161 ;146.250778;123.750511;101.249802];
%maillage raccord
for j=(k+1):1:(nbr_element_raccr+nbr_elem_virole+1)
    noeud(j).r=noeud(k).r+(rr-cos(angle_sect_maillage_raccr*(j-k))*rr) ;   %r_ext-L-(j-i)*(dc/nbr_elem_rayon);
    noeud(j).z=noeud(k).z-(sin(angle_sect_maillage_raccr*(j-k))*rr) ;        %noeud(j-1).z+dc/nbr_elem_rayon;
    elem(j-1).noeud1=j-1;
    elem(j-1).noeud2=j;
    elem(j-1).E1=E1;
    elem(j-1).E2=E2;
    elem(j-1).nu12=v12;
    elem(j-1).type=2;
    elem(j-1).epaiss=t;
    elem(j-1).angle=(atan((noeud(j).r-noeud(j-1).r)/abs((noeud(j).z-noeud(j-1).z))));
    elem(j-1).longueur=sqrt((noeud(j).r-noeud(j-1).r)^2+(noeud(j).z-noeud(j-1).z)^2);
    elem(j-1).rd=noeud(j).r;
    elem(j-1).rg=noeud(j-1).r;%Geom2 est la position radiale du deuxieme noeud de l element (inférieure)
end

%Maillage pied de bride
for i=(j+1):1:nbr_elem_pied_bride+nbr_element_raccr+nbr_elem_virole+1;
    noeud(i).r=noeud(j).r+(i-j)*(L/nbr_elem_pied_bride);
    r=noeud(i).r;
    noeud(i).z=0;
    elem(i-1).noeud1=i-1;
    elem(i-1).noeud2=i;
    elem(i-1).E1=E1;
    elem(i-1).E2=E2;
    elem(i-1).nu12=v12;
    elem(i-1).type=3;
    elem(i-1).epaiss=t;
    elem(i-1).longueur=L/nbr_elem_pied_bride;
    elem(i-1).angle=90;
    elem(i-1).rg=noeud(i-1).r;%Geom2 est la position radiale du deuxieme noeud de l element (inférieure)
    elem(i-1).rd=noeud(i).r
end


xsm1=[noeud(:).r]; % Cordoonnes en x des noeuds de maillage de la bride supérieure
ysm1=[noeud(:).z]; % Cordoonnes en y des noeuds de maillage de la bride supérieure

figure
scatter(xsm1,ysm1,'r')  %scatter
axis([r_int-100 r_ext+100 -100 h+100])
axis equal

figure
plot(xsm1,ysm1,'r-o')
axis([r_int-100 r_ext+100 -100 h+100])
axis equal
  
%%%%%%%%%%%%%%%%%%%%%%%%
% Matrices de raideur %%
%%%%%%%%%%%%%%%%%%%%%%%%

% [Kcoque_sym] = slanted_shell_stiffness_cisaillement_borne_1;
% [Kcoque_sym_v]=slanted_shell_stiffness_cisaillement_borne_1_virole;



Nombre_noeuds=length(noeud);
ddl=3;
taille=Nombre_noeuds*ddl;
Global_stiffness=sparse(taille,taille);
%   Matrices de raideurs des elements dans leurs repère global

for i =1:length(elem)
    
    switch (elem(i).type)
        case 3  %pied de bride
           E1=elem(i).E1;
            E3=elem(i).E2;
            v13=elem(i).nu12;
            v31=elem(i).nu12*E3/E1;
            t=elem(i).epaiss;
            L=taille_elt_pied_bride;
            phi=elem(i).angle*pi/180;
            rm=elem(i).rg+(elem(i).rd-elem(i).rg)/2;
            r2=(elem(i).rd-elem(i).rg)/2;
            
            
           [k_elem] = integral_gauss(E1,E3,v13,v31,rm,t,L,phi,r2);           
        
        case 2  %raccord
             E1=elem(i).E1;
            E3=elem(i).E2;  
            v13=elem(i).nu12;
            v31=elem(i).nu12*E3/E1;
            rm=elem(i).rg+(elem(i).rd-elem(i).rg)/2;
            r2=(elem(i).rd-elem(i).rg)/2;
            t=elem(i).epaiss;
            L=elem(i).longueur;
            phi=elem(i).angle;     
            [k_elem] = integral_gauss(E1,E3,v13,v31,rm,t,L,phi,r2);
            
            
        case 1  %virole
            E1=elem(i).E1;
            E3=elem(i).E2;
            v13=elem(i).nu12;
            v31=elem(i).nu12*E3/E1;
            t=elem(i).epaiss;
            L=elem(i).longueur;
            phi=elem(i).angle*pi/180;
            rm=elem(i).rg;
           
            
              [k_elem] = integral_gauss_v(E1,E3,v13,v31,rm,t,L,phi);
  
            
            
       
    end    
        for l=1:6;
            if(l<=3);
                a=3*(elem(i).noeud1-1)+l;
            else
                a=3*(elem(i).noeud2-1)+(l-3);
            end
            
            for m=1:6
                if(m<=3)
                    b=3*(elem(i).noeud1-1)+m;
                else
                    b=3*(elem(i).noeud2-1)+(m-3);
                end
                Global_stiffness(a,b)=Global_stiffness(a,b)+k_elem(l,m);
            end
        end
end
K=Global_stiffness;

%Création du vecteur effort
Effort=zeros(length(Global_stiffness),1);
Effort(length(Effort),:)=[];
Effort(length(Effort),:)=[];
Effort(length(Effort),:)=[];
K(size(K,1),:)=[];
K(size(K,1),:)=[];
K(size(K,1),:)=[];
K(:,size(K,2))=[];
K(:,size(K,2))=[];
K(:,size(K,2))=[];
%ici on retire le premier noeud de l'assemblage (CL:blocage noeud 1) car sinon la matrice de rigidité sera égale à 0
%% Application du chargement extérieur
fileID=fopen('Resultats_BrideSup_Metallique_raccord_4.dat','w');
fprintf(fileID,'%s\t%s\t%s\t\r\n','Fe(N)','U2_pt_effort(mm)','U2_pt_decol(mm)');
fclose(fileID);

% close figure(1)
k=1;
for j=0:Fmax/N:Fmax
    Effort(1)=-j;  %effort sur le dernier noeud dans la direction v
    % Effort(length(Effort)-2)=16000;
    F=Effort;
    U=K\F;
    U=[U(1:length(U));0;0;0]; %ici on remet les dep du premiers noeud à 0 pcq on l'a sup ligne 207 à 215 
    fileID=fopen('Resultats_BrideSup_Metallique_raccord_4.dat','a+');
    fprintf(fileID,'%d\t%.4f\t%.4f\t\r\n',j,U(length(U)-1),U(3*i-1));
    fclose(fileID);
    U1(k,1)=U(length(U)-2);
    U2(k,1)=U(length(U)-1);
    U3(k,1)=U(length(U));
   % xsm2=xsm1'+U(1:3:length(U))*amp; % Cordoonnes en x des noeuds de maillage de la bride supérieure
    %ysm2=ysm1'+U(2:3:length(U))*amp; % Cordoonnes en y des noeuds de maillage de la bride supérieure
   xsm2=xsm1'+U(2:3:length(U))*amp; % Cordoonnes en x des noeuds de maillage de la bride supérieure
   ysm2=ysm1'-U(1:3:length(U))*amp; % Cordoonnes en y des noeuds de maillage de la bride supérieure
    U1(k,1)=U(1);
    U2(k,1)=U(2);
    U3(k,1)=U(3);
    k=k+1;
end
U(length(U)-1);
figure
    hold on
    plot(xsm1,ysm1,'-','Linewidth',2)
    hold on
    plot(xsm2,ysm2,':','Linewidth',2)
    axis equal
    axis([r_int-100 r_ext+100 -300 h+100])
    title('Appercu déformation')
    hold off



%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EF Vs Semi-Analytique %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

FN1=sprintf('EF_1D_Axi_Bride_Metallique_Raccord_4_U1.dat');
FN2=sprintf('EF_1D_Axi_Bride_Metallique_Raccord_4.dat');
FN3=sprintf('EF_1D_Axi_Bride_Metallique_Raccord_4_U3.dat');
FN4=sprintf('Resultats_BrideSup_Metallique_raccord_4.dat');

R_EF_1D_U1=importdata(FN1);
R_EF_1D=importdata(FN2);
R_EF_1D_U3=importdata(FN3);
R_Sa=importdata(FN4);
R_SA=R_Sa.data;

% 
%%%%calcul erreur%%%%
%erreur_1D_axi=abs(R_SA(:,2)-R_EF_1D(:,2))./R_EF_1D(:,2);
% erreur_2D_axi=abs(R_SA(:,2)-R_EF_2D_axi(:,2))./R_EF_2D_axi(:,2);
% erreur_3D=abs(R_SA(:,2)-R_EF_3D(:,2))./R_EF_3D(:,2);
%R_SA(:,2) Deplacement 
%R_SA(:,1) Effort

figure %U2
plot(R_EF_1D(:,1)*1e-3*Fmax/N,R_EF_1D(:,2),'-.m','Linewidth',1)
hold on
plot(R_SA(:,1)*1e-3,-U1,'-.r','Linewidth',1)
xlabel('Effort extérieur (kN)')
ylabel('Déplacement axial (mm)')
% yyaxis right
% hold on
% plot(R_SA(:,1)*1e-3,100*erreur_1D(:,1))
% ylabel('Erreur (%)')
% axis ([0 5000 -5 100])
legend('EF-1D-axi','SA','Erreur par/ à EF-1D Axi','Location','northwest')
title('deplacement axial')
hold off
% 
figure %U1
plot(R_EF_1D(:,1)*1e-3*Fmax/N,R_EF_1D_U1(:,2),'-.m','Linewidth',1)
hold on
plot(R_SA(:,1)*1e-3,U2,'-.r','Linewidth',1)
xlabel('Effort extérieur (kN)')
ylabel('Déplacement (mm)')
% yyaxis right
% hold on
% plot(R_SA(:,1)*1e-3,100*erreur_1D_U1(:,1))
% ylabel('Erreur (%)')
% axis ([0 5000 -5 100])
legend('EF-1D-axi','SA','Erreur par/ à EF-1D Axi','Location','northwest')
title('deplacement radial')
hold off

figure %U3
plot(R_EF_1D(:,1)*1e-3*Fmax/N,R_EF_1D_U3(:,2),'-.m','Linewidth',1)
hold on
plot(R_SA(:,1)*1e-3,-U3,'-.r','Linewidth',1)
xlabel('Effort extérieur (kN)')
ylabel('Déplacement (mm)')
% yyaxis right
% hold on
% plot(R_SA(:,1)*1e-3,100*erreur_1D_U1(:,1))
% ylabel('Erreur (%)')
% axis ([0 5000 -5 100])
legend('EF-1D-axi','SA','Location','northwest')
title('rotation')
hold off
t=toc
