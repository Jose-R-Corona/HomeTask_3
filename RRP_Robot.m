clear all
close all

%symbolical variables
syms q1 q2 q3  d1 a2  t real

% set angles 
q = [q1 q2 q3];
q_test = [0 pi/4 1];

%Link lengths
L=[d1 a2];
L_test=[10, 10];


%%
%FK
FK=Rz(q(1))*Tz(L(1))*Ry(-q(2))*Tx(L(2))*Tx(q(3));
simplify(FK);

T=double(subs(FK, [q L], [q_test L_test]))

draw_robot(q_test,L_test)


%%
%IK
Pos=T(1:3,4) % extract position of EE
%Pos=[Px Py Pz]
%L_test=[d1 a2];
l_p=sqrt(Pos(1)^2+Pos(2)^2);

%%%%%%%%%%%%%%%%%% 1st solution
q_11=atan2(Pos(2),Pos(1));
q_3=sqrt((Pos(3)-L_test(1))^2+l_p^2)-L_test(2);
q_21=atan2((Pos(3)-L_test(1)),l_p);

q_IK_1=[q_11 q_21 q_3];

%%%%%%%%%%%%%%%% 2nd solution
q_12=pi+atan2(Pos(2),Pos(1));
q_22=pi-atan2((Pos(3)-L_test(1)),l_p);

q_IK_2=[q_12 q_22 q_3];


T_ik_21=double(subs(FK, [q L], [q_IK_1 L_test] ))
T_ik_22=double(subs(FK, [q L], [q_IK_2 L_test] ))

%%
%Jacobian Numerial method

H = Rz(q(1))*Tz(L(1))*Ry(-q(2))*Tx(L(2))*Tx(q(3));% forward kinematics
H=simplify(H);
R = simplify(H(1:3,1:3));  % extract rotation matrix
% diff by q1
Td=Rzd(q(1))*Tz(L(1))*Ry(-q(2))*Tx(L(2))*Tx(q(3))*...
    [R^-1 zeros(3,1);0 0 0 1];
J_1 = [Td(1,4), Td(2,4), Td(3,4), Td(3,2), Td(1,3), Td(2,1)]' ; % extract 6 components from 4x4 Td matrix to Jacobian 1st column
% diff by q2
Td=Rz(q(1))*Tz(L(1))*Ryd(-q(2))*Tx(L(2))*Tx(q(3))*...
    [R^-1 zeros(3,1);0 0 0 1];
J_2 = [Td(1,4), Td(2,4), Td(3,4), Td(3,2), Td(1,3), Td(2,1)]' ; % extract 6 components from 4x4 Td matrix to Jacobian 1st column
% diff by q3
Td=Rz(q(1))*Tz(L(1))*Ry(-q(2))*Tx(L(2))*Txd(q(3))*...
    [R^-1 zeros(3,1);0 0 0 1];
J_3 = [Td(1,4), Td(2,4), Td(3,4), Td(3,2), Td(1,3), Td(2,1)]' ; % extract 6 components from 4x4 Td matrix to Jacobian 1st column

% Full Jacobian 6x3
Jq1 = [simplify(J_1), simplify(J_2), simplify(J_3)]

%%
%Jacobian: Screw Theory
%Finds FKs for each joint
T00= eye(4);  
T01= Rz(q(1))*Tz(L(1));
T02= Rz(q(1))*Tz(L(1))*Ry(-q(2))*Tx(L(2));
T03= Rz(q(1))*Tz(L(1))*Ry(-q(2))*Tx(L(2))*Tx(q(3));

%Find Origins
O0 = T00(1:3,4);
O1 = T01(1:3,4);
O2 = T02(1:3,4);
O3 = T03(1:3,4);

%rotation (translation in case of prismatic joint) axis Z from transformation
Z0 = T00(1:3,3); % 3rd coloumn corresponds to Rz
Z1 = T01(1:3,2); % 2nd coloumn corresponds to Ry
Z2 = T02(1:3,1); % 1rd coloumn corresponds to Tx

J_1 = [cross(Z0,(O3-O0));Z0];
J_2 = [cross(Z1,(O3-O1));Z1];
J_3 = [Z2;0; 0; 0];

% Full Jacobian 6x3
Jq2 = [simplify(J_1), simplify(J_2), simplify(J_3)]

simplify(Jq1-Jq2)

%%
%Jacobian: clasical derivate method
H = Rz(q(1))*Tz(L(1))*Ry(-q(2))*Tx(L(2))*Tx(q(3));% forward kinematics

Px=H(1,4);
Py=H(2,4);
Pz=H(3,4);

%w=[-q(2)*sin(q(1)),q(2)*cos(q(1)),q(1)]';

%j=[Jlinear velocity ; J angular velocity]
J_1 = [diff(Px,q(1)),diff(Py,q(1)),diff(Pz,q(1)),0,0,1]';

%since -q2
J_2=[-diff(Px,q(2)),-diff(Py,q(2)),-diff(Pz,q(2)),-sin(q(1)),cos(q(1)),0]';

J_3=[diff(Px,q(3)),diff(Py,q(3)),diff(Pz,q(3)),0,0,0]';

% Full Jacobian 6x3
Jq3 = [simplify(J_1), simplify(J_2), simplify(J_3)]

simplify(Jq2-Jq3)

%%
%Singularities
%[R,p] = rref(Jq2);
R = rref(Jq2);

S_1=double(subs(Jq2, [q L], [[pi/2 pi/2 1] L_test] ))
rank(S_1)   %is always 3
%det(S_1)

%%
%velocities
q_position1 = sin(t); q_position2 = cos(2*t); q_position3 = sin(3*t);
q_position = [q_position1,q_position2,q_position3];

v_q = diff(q_position,t)'
v_tool=Jq1*(v_q)
T = 0:pi/100:2*pi;
L_test=[10, 10];

V_tool= subs(v_tool, [q L], [q_position L_test] );
V_tool=double(subs(V_tool, t, T ));

figure(2)
plot(T,V_tool(1,:),T,V_tool(2,:),T,V_tool(3,:))
title('Cartesian Velocities of tool frame')
xlabel('0 < t < 2pi') 
ylabel('Velocity') 
legend({'Vx tool','Vy tool','Vz tool'},'Location','southwest')

figure(3)
plot(T,V_tool(4,:),T,V_tool(5,:),T,V_tool(6,:))
title('Angular Velocities of tool frame')
xlabel('0 < t < 2pi') 
ylabel('Velocity') 
legend({'Wx tool','Wy tool','Wz tool'},'Location','southwest')


Q_position=subs(q_position',t,T);
figure(4)
plot(T,Q_position(1,:),T,Q_position(2,:),T,Q_position(3,:))
title('Joints Possition')
xlabel('0 < t < 2pi') 
ylabel('Possition') 
legend({'q1','q2','q3'},'Location','southwest')