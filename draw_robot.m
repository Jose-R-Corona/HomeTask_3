function draw_robot(q,L)

figure;
hold on
view(200,25)
grid on
axis equal
T1=Rz(q(1))*Tz(L(1));
plot3(0,0,0,'ro','MarkerSize',10,'LineWidth', 10);
plot3([0 T1(1,4)],[0 T1(2,4)],[0 T1(3,4)],'-b','LineWidth', 5);

T2=Rz(q(1))*Tz(L(1))*Ry(-q(2))*Tx(L(2));
plot3(T1(1,4),T1(2,4),T1(3,4),'ro','MarkerSize',10,'LineWidth', 10);
plot3([T1(1,4) T2(1,4)],[T1(2,4) T2(2,4)],[T1(3,4) T2(3,4)],'-b','LineWidth', 10);


T3=Rz(q(1))*Tz(L(1))*Ry(-q(2))*Tx(L(2))*Tx(q(3));
plot3(T2(1,4),T2(2,4),T2(3,4),'ro','MarkerSize',10,'LineWidth', 10);
plot3([T2(1,4) T3(1,4)],[T2(2,4) T3(2,4)],[T2(3,4) T3(3,4)],'-r','LineWidth', 5);
plot3(T3(1,4),T3(2,4),T3(3,4),'ro','MarkerSize',10,'LineWidth', 10);


