function T=FK_stanford(q,L)
T=(Rz(q(1))*Tz(L(1))*Tx(-L(2))*Rx(q(2))*Tz(L(3))*Tz(q(3))*Rz(q(4))*Rx(q(5))*Rz(q(6)));
end