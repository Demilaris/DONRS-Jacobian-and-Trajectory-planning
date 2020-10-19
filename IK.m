P = [1 0 1]
B = IK(P)


function q = IK(p)

l1 = 1;
l2 = 1;
l3 = 1;
p_x=p(1);p_y=p(2);p_z=p(3);

theta_1=atan2(p_y,p_x);
s=sqrt(p_x^2+p_y^2);

theta_3 = acos( ( s^2 + (l1-p_z)^2 - l2^2 - l3^2)/ (2 * l2 * l3) );
theta_2 = - atan2( l3 * sin(theta_3), (l2 + l3 * cos(theta_3))) + atan2((l1-p_z),s);
q = [theta_1 theta_2 theta_3];

end
