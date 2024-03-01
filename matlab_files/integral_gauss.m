function [k_elem] = integral_gauss(E1,E3,v13,v31,rm,t,L,phi,r2)

 
%%%parametre de gauss 2 points

w1=1;
w2=1;
x1=-0.5773502691896257;
x2=0.5773502691896257;
[ kp1] = slanted_shell_stiffness_cisaillement_borne_1(x1,E1,E3,v13,v31,rm,t,L,phi,r2);
[ kp2] = slanted_shell_stiffness_cisaillement_borne_1(x2,E1,E3,v13,v31,rm,t,L,phi,r2);
k_elem=w1*kp1+w2*kp2;


end

