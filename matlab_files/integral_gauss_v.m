function [k_elem] = integral_gauss_v(E1,E3,v13,v31,rm,t,L,phi)

 
%%%parametre de gauss 2 points

w1=1;
w2=1;
x1=-0.5773502691896257;
x2=0.5773502691896257;
[kpv1]=slanted_shell_stiffness_cisaillement_borne_1_virole(x1,E1,E3,v13,v31,rm,t,L,phi);
[kpv2]=slanted_shell_stiffness_cisaillement_borne_1_virole(x2,E1,E3,v13,v31,rm,t,L,phi);

k_elem=w1*kpv1+w2*kpv2;

