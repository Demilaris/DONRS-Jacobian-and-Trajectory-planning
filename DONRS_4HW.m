clear all
syms theta_1 theta_2 theta_3 real
l1 = 1;
l2 = 1;
l3 = 1;

%FK
H = Rz(theta_1)*Tz(l1)*Ry(theta_2)*Tx(l2)*Ry(theta_3)*Tx(l3);
H = simplify(H)


x1 = cos(theta_1) * (l2 * cos(theta_2) + l3 * cos(theta_2 * theta_3));
y1 = sin(theta_1)*(l2 * cos(theta_2) + l3 * cos(theta_2 * theta_3));
z1 = l1 - l2*sin(theta_2) - l3*sin(theta_2*theta_3)
%numerical
R = simplify(H(1:3,1:3));
Td = Rzd(theta_1)*Tz(l1)*Ry(theta_2)*Tx(l2)*Ry(theta_3)*Tx(l3)*[R^-1 zeros(3,1);0 0 0 1];
Jn1 = [Td(1,4), Td(2,4),Td(3,4),Td(3,2),Td(1,3),Td(2,1)]';

Td1 = Rz(theta_1)*Tz(l1)*Ryd(theta_2)*Tx(l2)*Ry(theta_3)*Tx(l3)*[R^-1 zeros(3,1);0 0 0 1];
Jn2 = [Td1(1,4), Td1(2,4),Td1(3,4),Td1(3,2),Td1(1,3),Td1(2,1)]';

Td2 = Rz(theta_1)*Tz(l1)*Ry(theta_2)*Tx(l2)*Ryd(theta_3)*Tx(l3)*[R^-1 zeros(3,1);0 0 0 1];
Jn3 = [Td2(1,4), Td2(2,4),Td2(3,4),Td2(3,2),Td2(1,3),Td2(2,1)]';

Jacobian_numerical = [simplify(Jn1),simplify(Jn2),simplify(Jn3)];

%Joint trajectory q(t) from q(0) = (0, 0, 0) to q(2) = (2, 3, 4) 
clear all
close all
t0 = 0; tf = 2;
q0 = 0; qf = 2; 
v0 = 0; vf = 0;
acc0 = 0; accf = 0;
A = MAT6_a(t0,tf);
c = [q0;v0;acc0;qf;vf;accf];
b = A\c
c1 = [q0;v0;acc0;3;vf;accf];
b1 = A\c1
c2 = [q0;v0;acc0;4;vf;accf];
b2 = A\c2
 
a0 = b(1); a1 = b(2); a2 = b(3); a3 = b(4); a4 = b(5); a5 = b(6);
a10 = b1(1); a11 = b1(2); a12 = b1(3); a13 = b1(4); a14 = b1(5); a15 = b1(6);
a20 = b2(1); a21 = b2(2); a22 = b2(3); a23 = b2(4); a24 = b2(5); a25 = b2(6);
t0 = 0 ; tf = 2;
t = t0:0.1:tf;
q = a0+a1.*t+a2.*t.^2+a3.*t.^3+a4.*t.^4+a5.*t.^5;
v = a1+2*a2.*t+3*a3.*t.^2+4*a4.*t.^3+5*a5.*t.^4;
acc = 2*a2+6*a3.*t+12*a4.*t.^2+20*a5.*t.^3;
 
q1 = a10+a11.*t+a12.*t.^2+a13.*t.^3+a14.*t.^4+a15.*t.^5;
v1 = a11+2*a12.*t+3*a13.*t.^2+4*a14.*t.^3+5*a15.*t.^4;
acc1 = 2*a12+6*a13.*t+12*a14.*t.^2+20*a15.*t.^3;
 
q2 = a20+a21.*t+a22.*t.^2+a23.*t.^3+a24.*t.^4+a25.*t.^5;
v2 = a21+2*a22.*t+3*a23.*t.^2+4*a24.*t.^3+5*a25.*t.^4;
acc2 = 2*a22+6*a23.*t+12*a24.*t.^2+20*a25.*t.^3;

 figure
 plot(t,q,'r-')
 hold on
 plot(t,q1,'g-')
  hold on
 plot(t,q2,'b-')
 title('position vs time')
 legend('joint_1','joint_2','joint_3')
 grid on
 
 figure
 plot(t,v,'r-')
 hold on
 plot(t,v1,'g-')
 hold on
 plot(t,v2,'b-')
 title('velocity vs time')
 legend('joint_1','joint_2','joint_3')
 grid on
 
 figure
 plot(t,acc,'r-')
 hold on
 plot(t,acc1,'g-')
 plot(t,acc2,'b-')
 hold on
 title('acceleration vs time')
 legend('joint_1','joint_2','joint_3')
 grid on
 %3 task
 q10 = 0; q1f = 2; v1 = 1; a1 = 10;
q20 = 0; q2f = 3; v2 = 1; a2 = 10;
q30 = 0; q3f = 4; v3 = 1; a3 = 10;
v10 = 0; v20 = 0; v30 = 0;
delta_t = 0.01;
n = 0;

while (floor(delta_t*10^n)~=delta_t*10^n)
    n=n+1;
end
E = 1*10^-n;

t1a = v1/a1;
if rem(t1a,delta_t)~=0
    t1a_new = round(t1a,n)+E;
else
    t1a_new = round(t1a,n);
end
t1f = (q1f-q10)/v1 + t1a_new;
if rem(t1f,delta_t)~=0
    t1f_new = round(t1f,n)+E;
else
    t1f_new = round(t1f,n);
end

t2a = v2/a2;
if rem(t2a,delta_t)~=0
    t2a_new = round(t2a,n)+E;
else
    t2a_new = round(t2a,n);
end
t2f = (q2f-q20)/v2 + t2a_new;
if rem(t2f,delta_t)~=0
    t2f_new = round(t2f,n)+E;
else
    t2f_new = round(t2f,n);
end

t3a = v3/a3;
if rem(t3a,delta_t)~=0
    t3a_new = round(t3a,n)+E;
else
    t3a_new = round(t3a,n);
end
t3f = (q3f-q30)/v3 + t3a_new;
if rem(t3f,delta_t)~=0
    t3f_new = round(t3f,n)+E;
else
    t3f_new = round(t3f,n);
end

Maximtf = [t1f_new t2f_new t3f_new]
tf_new = max(Maximtf);
ta_new = t3a_new

v1_new = ((q1f-q10)/(tf_new-ta_new));
a1_new = v1_new/ta_new;

v2_new = ((q2f-q20)/(tf_new-ta_new));
a2_new = v2_new/ta_new;

v3_new = ((q3f-q30)/(tf_new-ta_new));
a3_new = v3_new/ta_new;

% joint 1 - coefficients:
% t0 --> ta:
a10 = q10;
a11 = v10;
a12 = 0.5*a1_new;

% tf-ta --> tf:
a20 = q10 + 0.5*a1_new*ta_new^2 - v1_new*ta_new;
a21 = v1_new;

% tf-ta --> tf:
a30 = q1f - 0.5*a1_new*tf_new^2;
a31 = a1_new*tf_new;
a32 = -0.5*a1_new;


% joint 2 - coefficients:
% t0 --> ta:
b10 = q20;
b11 = v20;
b12 = 0.5*a2_new;


% ta --> tf-ta:
b20 = q20 + 0.5*a2_new*ta_new^2 - v2_new*ta_new;
b21 = v2_new;

% tf-ta --> tf:
b30 = q2f - 0.5*a2_new*tf_new^2;
b31 = a2_new*tf_new;
b32 = -0.5*a2_new;


% joint 3 - coefficients:
% t0 --> ta:
c10 = q30;
c11 = v30;
c12 = 0.5*a3_new;

% ta --> tf-ta:
c20 = q30 + 0.5*a3_new*ta_new^2 - v3_new*ta_new;
c21 = v3_new;

% tf-ta --> tf:
c30 = q3f - 0.5*a3_new*tf_new^2;
c31 = a3_new*tf_new;
c32 = -0.5*a3_new;



t = 0:delta_t:tf_new;
q1 = (a10+a11.*t+a12.*t.^2).*(t<=ta_new)...
    +(a20+a21.*t).*(t>ta_new).*(t<=(tf_new-ta_new))...
    +(a30+a31.*t+a32.*t.^2).*(t>(tf_new-ta_new)).*(t<=tf_new);
v1 = (a11+2*a12.*t).*(t<=ta_new)...
    +(a21).*(t>ta_new).*(t<=(tf_new-ta_new))...
    +(a31+2*a32.*t).*(t>(tf_new-ta_new)).*(t<=tf_new);
acc1 = (2*a12).*(t<=ta_new)...
    +(0).*(t>ta_new).*(t<=(tf_new-ta_new))...
    +(2*a32).*(t>(tf_new-ta_new)).*(t<=tf_new);


q2 = (b10+b11.*t+b12.*t.^2).*(t<=ta_new)...
    +(b20+b21.*t).*(t>ta_new).*(t<=(tf_new-ta_new))...
    +(b30+b31.*t+b32.*t.^2).*(t>(tf_new-ta_new)).*(t<=tf_new);
v2 = (b11+2*b12.*t).*(t<=ta_new)...
    +(b21).*(t>ta_new).*(t<=(tf_new-ta_new))...
    +(b31+2*b32.*t).*(t>(tf_new-ta_new)).*(t<=tf_new);
acc2 = (2*b12).*(t<=ta_new)...
    +(0).*(t>ta_new).*(t<=(tf_new-ta_new))....
    +(2*b32).*(t>(tf_new-ta_new)).*(t<=tf_new);


q3 = (c10+b11.*t+c12.*t.^2).*(t<=ta_new)...
    +(c20+c21.*t).*(t>ta_new).*(t<=(tf_new-ta_new))...
    +(c30+c31.*t+c32.*t.^2).*(t>(tf_new-ta_new)).*(t<=tf_new);
v3 = (c11+2*c12.*t).*(t<=ta_new)...
    +(c21).*(t>ta_new).*(t<=(tf_new-ta_new))...
    +(c31+2*c32.*t).*(t>(tf_new-ta_new)).*(t<=tf_new);
acc3 = (2*c12).*(t<=ta_new)...
    +(0).*(t>ta_new).*(t<=(tf_new-ta_new))....
    +(2*c32).*(t>(tf_new-ta_new)).*(t<=tf_new);

figure
plot(t,q1,'r')
hold on
plot(t,q2,'b')
grid on
plot(t,q3,'g')
grid on
title('position vs time')
legend('joint_1', 'joint_2','joint_3')
axis([0 tf_new -inf inf])

figure
plot(t,v1,'r')
hold on
plot(t,v2,'b')
hold on
plot(t,v3,'g')
title('velocity vs time')
legend('joint_1' , 'joint_2','joint_3')
grid on
axis([0 tf_new -inf inf])

figure
plot(t,acc1,'r')
hold on
plot(t,acc2,'b')
hold on
plot(t,acc3,'g')
title('acceleration vs time')
legend('joint_1','joint_2','joint_3')
grid on
axis([0 tf_new -inf inf])

