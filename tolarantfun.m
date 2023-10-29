function dx=tolarantfun(t,x)
global S  m d   k1 I b   k2  yr yd1 yd2 yd3 e gama1 liao r u1 gama2 p
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1=x(1);
x2=x(2);
x3=x(3);
x4=x(4);
x5=x(5);

x(16);x(17);
x6=[x(16);x(17)];
theta21=x(6:10)';
theta22=x(11:15)';
theta2=[theta21;theta22];
% p1=20*exp(-0.1*t)+1.53;
% p=0.05*[1,1,1]';%%p是已知的矩阵3*1
X1=[x1,x2,x3]';
X2=[x4,x5]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yd1=2*sin(0.5*t);
dyd1=1*cos(0.5*t);
ddyd1=-0.5*sin(0.5*t);
yd2=2*cos(0.5*t);
dyd2=-1*sin(0.5*t);
ddyd2=-0.5*cos(0.5*t);
yd3=-0.5*t;
dyd3=-0.5;
ddyd3=0;
yr=[yd1,yd2,yd3]';
e=X1-yr;
dyr=[dyd1,dyd2,dyd3]';
ddyr=[ddyd1,ddyd2,ddyd3]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S=[cos(x3),0;sin(x3),0;0,1];
dS=[-sin(x3)*x5,0;cos(x3)*x5,0;0,0];
invHE=[1/(2*gama1) 1/(2*b*gama1);1/(2*gama1) -1/(2*b*gama1)];
HM=[m 0;0 d^2*m+I];
HV=[0 -d*m*x5;d*m*x5 0];
invHM=[1/m 0;0 1/(d^2*m+I)];
% rou1=20*exp(-0.1*t)+1.53;
% drou1=-2*exp(-0.1*t)+1.53;
% rou=rou1*[1 1 1]';
z1=e;
% p=dyr+(drou1/rou1)*e;
alpha1=(-S')*(k1*z1-dyr);
dalpha1=-(dS)'*(k1*z1-dyr)-S'*(k1*(S*X2-dyr)-ddyr);
z2=X2-alpha1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:5
    c2{i}=[-3+i;-3+i]';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:5
    fai2(i)=exp(-(X2-c2{i}')'*(X2-c2{i}')*0.2);    

end
% alpha2=-k2*z2-(theta2)*(fai2)'-0.5*z2-(S')*z1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u=[x(16),x(17)].*(invHE)*HM*(invHM*HV*X2+dalpha1-(theta2)*(fai2)'-S'*z1-k2*z2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(t)
   if (t(i)>=10|t(i)<=20) & (t(i)>=30|t(i)<=40) & (t(i)>=50|t(i)<=60) & (t(i)>=70|t(i)<=80)
    uc1(i)=0.1*u(i);%1.4
else
    uc1(i)=1*u(i);%0.5
end 

if (t(i)>=10|t(i)>=20) & (t(i)>=30|t(i)<=40) & (t(i)>=50|t(i)<=60) & (t(i)>=70|t(i)<=80)
    uc2(i)=0.1*u(i);%0.95
else
    uc2(i)=1*u(i);
end  
end

C=[1 0];C1=[1 0 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dX21=d*x5*x5+(1/m)*(x(16).*uc1+x(17).*uc2);
dX22=-d*m/(d^2*m+I)*x4+b/(d^2*m+I)*(u(1)-u(2));
liao=[liao, dX21];
dtheta2=gama2*z2*fai2-k2*theta2; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx=[x4*cos(x3),x4*sin(x3),x5,dX21,dX22,dtheta2(1,:),dtheta2(2,:),-0.001*x(16),-0.001*x(17)]';
end