function q=IK_stanford(T,L)

Pos=T(1:3,4);

dPos=sqrt(Pos(1)^2+Pos(2)^2);
dL3_p=sqrt(dPos^2-L(2)^2);
q11=pi+atan2(Pos(2),Pos(1))-atan2(dL3_p,L(2));
q12=pi+atan2(Pos(2),Pos(1))-atan2(-dL3_p,L(2));

dZ=Pos(3)-L(1);
q31=sqrt(dL3_p^2+dZ^2)-L(3);

q21=atan2(dL3_p,dZ);
q22=atan2(-dL3_p,dZ);

T123=Rz(q11)*Tz(L(1))*Tx(-L(2))*Rx(q21)*Tz(L(3))*Tz(q31);
T456bis=inv(T123)*T;
if abs(T456bis(3,3))~=1
    q41=atan2(T456bis(1,3),-T456bis(2,3));
    q61=atan2(T456bis(3,1),T456bis(3,2));
    q51 = atan2(sqrt(T456bis(1,3)^2+T456bis(2,3)^2),T456bis(3,3));
else
    disp('singular')
end

T123=Rz(q12)*Tz(L(1))*Tx(-L(2))*Rx(q22)*Tz(L(3))*Tz(q31);
T456bis=inv(T123)*T;
if abs(T456bis(3,3))~=1;
    q42=atan2(T456bis(1,3),-T456bis(2,3));
    q62=atan2(T456bis(3,1),T456bis(3,2));
    q52 = atan2(sqrt(T456bis(1,3)^2+T456bis(2,3)^2),T456bis(3,3));
else
    disp('singular')
end


T123=Rz(q11)*Tz(L(1))*Tx(-L(2))*Rx(q21)*Tz(L(3))*Tz(q31);
T456bis=inv(T123)*T;
if abs(T456bis(3,3))~=1
    q43=atan2(-T456bis(1,3),T456bis(2,3));
    q63=atan2(-T456bis(3,1),-T456bis(3,2));
    q53 = atan2(-sqrt(T456bis(1,3)^2+T456bis(2,3)^2),T456bis(3,3));
else
    disp('singular')
end

T123=Rz(q12)*Tz(L(1))*Tx(-L(2))*Rx(q22)*Tz(L(3))*Tz(q31);
T456bis=inv(T123)*T;
if abs(T456bis(3,3))~=1
    q44=atan2(-T456bis(1,3),T456bis(2,3));
    q64=atan2(-T456bis(3,1),-T456bis(3,2));
    q54 = atan2(-sqrt(T456bis(1,3)^2+T456bis(2,3)^2),T456bis(3,3));
else
    disp('singular')
end

q=[q11 q21 q31 q41 q51 q61;
    q12 q22 q31 q42 q52 q62;
    q11 q21 q31 q43 q53 q63;
    q12 q22 q31 q44 q54 q64];

end