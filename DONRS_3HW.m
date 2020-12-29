clear all
syms theta_1 theta_2 d1 a2 d3 real
%tests
d1_t = 5;
a2_t = 15;
d3_t = 10;

q1_t = -10*pi/180;
%check singulatity
q2_t = -pi/2 +pi;


%FK
H = Rz(theta_1)*Tz(d1)*Ry(-theta_2)*Tx(a2)*Tx(d3);
H = simplify(H);

%for test
Test_1 = Rz(q1_t)*Tz(d1_t)*Ry(-q2_t)*Tx(a2_t)*Tx(d3_t);

px_t = Test_1(1,4);
py_t = Test_1(2,4);
pz_t = Test_1(3,4);

px = H(1,4);
py = H(2,4);
pz = H(3,4);

%IK
r = sqrt(px_t^2+py_t^2);
s = pz_t - d1_t;

d3_k = sqrt(r^2+s^2) - a2_t;

theta_1_t = atan2(py_t, px_t);
theta_2_t = atan2(s,r);

%Test IK
T2 = Rz(theta_1_t)*Tz(d1_t)*Ry(-theta_2_t)*Tx(a2_t)*Tx(d3_k);

%Jacob_class

D = H(1:3,4);

Jacobian_diff = [diff(px, theta_1), -diff(px, theta_2), diff(px, d3);
                diff(py, theta_1), -diff(py, theta_2), diff(py, d3);
                diff(pz, theta_1), -diff(pz, theta_2), diff(pz, d3);];  

J_rot = [0, -sin(theta_1), 0;
        0, cos(theta_1), 0;
        1, 0, 0];
Jacobian_class = simplify([Jacobian_diff; J_rot])

%Jacob_skew_theory
T00 = eye(3);
T01 = Rz(theta_1)*Tz(d1);
T02 = Rz(theta_1)*Tz(d1)*Ry(-theta_2)*Tx(a2);
T03 = Rz(theta_1)*Tz(d1)*Ry(-theta_2)*Tx(a2)*Tx(d3);

o0 = T00(1:3,3);
o1 = T01(1:3,4);
o2 = T02(1:3,4);
o3 = T03(1:3,4);

z0 = T00(1:3,3);
z1 = T01(1:3,2);
z2 = T02(1:3,1);

j1 = [cross(z0,(o3-o0));z0];
j2 = [cross(z1,(o3-o1));z1];
j3 =[z2;0;0];

Jacobian_geom = [simplify(j1),simplify(j2),simplify(j3)]

%Numerical
R = simplify(H(1:3,1:3));
Td = Rzd(theta_1)*Tz(d1)*Ry(-theta_2)*Tx(a2)*Tx(d3)*[R^-1 zeros(3,1);0 0 0 1];
Jn1 = [Td(1,4), Td(2,4),Td(3,4),Td(3,2),Td(1,3),Td(2,1)]';

Td1 = Rz(theta_1)*Tz(d1)*Ryd(-theta_2)*Tx(a2)*Tx(d3)*[R^-1 zeros(3,1);0 0 0 1];
Jn2 = [Td1(1,4), Td1(2,4),Td1(3,4),Td1(3,2),Td1(1,3),Td1(2,1)]';

Td2 = Rz(theta_1)*Tz(d1)*Ry(-theta_2)*Tx(a2)*Txd(d3)*[R^-1 zeros(3,1);0 0 0 1];
Jn3 = [Td2(1,4), Td2(2,4),Td2(3,4),Td2(3,2),Td2(1,3),Td2(2,1)]';

Jacobian_numerical = [simplify(Jn1),simplify(Jn2),simplify(Jn3)]


%Singularities
Singularity = Jacobian_geom((1:3),:);
Determ_singularity = simplify(det(Singularity))

%Compute the velocity of the tool frame 
syms t theta_1_t theta_2_t d3_t real

theta_1_t = sin(t);
theta_2_t = cos(2*t);
d3_t = sin(3*t);

J_in_time = simplify(subs(Jacobian_geom, {d1, a2, theta_1, theta_2, d3}, {1, 1, theta_1_t, theta_2_t, d3_t}));

q = [theta_1_t theta_2_t d3_t]';
q_diff = diff(q);

xi = simplify(J_in_time * q_diff);

time = 0:0.1:10;
xi_in_time = subs(xi, {t}, {time});
%Linear velocity plot

figure;
plot(time, xi_in_time(1:3, :))
xlabel('Time')
ylabel('Velocity value')
title('Linear velocity plot')
grid on;
%Angular velocity plot
figure;
plot(time, xi_in_time(4:6, :))
xlabel('Time')
ylabel('Velocity value')
title('Angular velocity plot')
grid on;






