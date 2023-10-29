clear all
close all
clc
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global  M C  E  m e d  k1 k2 gama1  HM HV dS HE X1 X2 I b  yr yd1 yd2 yd3 gama2 r
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% r=200*diag([1,1,1],0);%%r是已知的偏导数，是一个常�?
M1=0.1*[1 1 1 1 1]';%%自�?�应�?1%0.1
M2=0.1*[1 1 1 1 1]';%%自�?�应�?2%0.1
gama2=3;%%可调参数%3，错误参数是
tic
gama1=0.1;%%可调参数%0.1（系统矩阵E）错�?30
k11=0.43;%0.43
k12=0.43;%0.43
k13=4;%4
k21=2;%2
k22=4;%4
k1=diag([k11,k12,k13],0);%%李亚普函数前面的系数�?
k2=diag([k21,k22],0);%%李亚普函数前面的系数�?
d=0.3;%%车中心到后轮中心的距�?%0.3
m=1;%%质量�?1kg�?
b=0.2;%0.2
I=2.6;%2.6
toc
disp(['运行时间1: ',num2str(toc)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
x0=[0.1;1.1;3.4;0.2;0.5;M1;M2;1;1;];%[0.1;1.1;3;1;0.5;M1;M2]�?
t0=0;tf=100;tspan=t0:0.01:tf;
[t,x]=ode45(@tolarantfun,tspan,x0);   
toc
disp(['运行时间2: ',num2str(toc)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1=x(:,1);
x2=x(:,2);
x3=x(:,3);
x4=x(:,4);
x5=x(:,5);
theta21=x(:,6:10)';
theta22=x(:,11:15)';
x(:,16);x(:,17);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(t)
    theta2{i}=[theta21(:,i),theta22(:,i)];
end
theta2_1=[];theta2_2=[];
for i=1:length(t)
    theta2_1=[theta2_1 theta2{i}(1)];
    theta2_2=[theta2_2 theta2{i}(2)];
end
X1=[x1,x2,x3]';
X2=[x4,x5]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yd1=2*sin(0.5*t);
dyd1=1*cos(0.5*t);
yd2=2*cos(0.5*t);
dyd2=-1*sin(0.5*t);
yd3=-0.5*t;
dyd3=[];
for i=1:length(t)
dyd3_1=-0.5;
dyd3=[dyd3; dyd3_1];
end
yr=[yd1,yd2,yd3]';
dyr=[dyd1,dyd2,dyd3]';
ddyd1=-0.5*sin(0.5*t);
ddyd2=-0.5*cos(0.5*t);
ddyd3=[];
for i=1:length(t)
ddyd3_1=0;
ddyd3=[ddyd3; ddyd3_1];
end
ddyr=[ddyd1,ddyd2,ddyd3]';
e=X1-yr;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
for j=1:length(t)
   for i=1:5
    c2{i,j}=[-3+i;-3+i]';
    fai2(i,j)=exp(-(X2(j)-c2{i,j}')'*(X2(j)-c2{i,j}')*0.2);    
   end
end
zz1=[];zz2=[];uu=[];
for i=1:length(t)
    S=[cos(x3(i)),0;sin(x3(i)),0;0,1];
    AA=[-sin(x3(i)),cos(x3(i)),0];
    M=[m,0,-d*m*sin(x3(i));0,m,d*m*cos(x3(i));-d*m*sin(x3(i)),d*m*cos(x3(i)),m*(d^2)+1];
    C=[0,0,-d*m*x5(i)*cos(x3(i));0,0,-d*m*x5(i)*sin(x3(i));0,0,0];
    E=[cos(x3(i)),cos(x3(i));sin(x3(i)),sin(x3(i));b,-b]/gama1;
    HM=[m 0;0 d^2*m+I];%S'*M*S;
    dS=[-sin(x3(i)),0;cos(x3(i)),0;0,0];
    HV=[0 -d*m*x5(i);d*m*x5(i) 0];%S'*(M*dS+C*S);
    HE=S'*E;
    invHE=[1/(2*gama1) 1/(2*b*gama1);1/(2*gama1) -1/(2*b*gama1)];
    invHM=[1/m 0;0 1/(d^2*m+I)];
    z1=X1-yr;
%     zz1=[zz1 z1];
    alpha1=(-S')*(k1*z1(:,i)-dyr);
    z2=X2-alpha1;
%     zz2=[zz2 z2];
    dalpha1=-(dS)'*(k1*z1(:,i)-dyr(:,i))-S'*(k1*(S*X2(:,i)-dyr(:,i))-ddyr(:,i));
    H=invHE*HM*invHM*HV;
    uu_1=[x(16),x(17)].*invHE*HM*(invHM*HV*X2(:,i)+dalpha1-S'*z1(:,i)-k2*z2(:,i));
%     uu_1=0.5*zz2(:,j)*(alpha22)'*alpha22/sqrt(0.5*(zz2(:,j))'*zz2(:,j)*(alpha22)'*alpha22+6);
    uu=[uu uu_1];
end
toc
disp(['运行时间3: ',num2str(toc)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 定义已知数据�?
x_1known = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100];
y_1known=[5	3	7	2	7	5	6	7	3	4	4	6	6	2	3	2	3	4	4	3	3	7	6	4	7	2	6	6	5	3	5	3	2	3	7	2	3	2	4	2	7	3	2	3	4	2	7	3	3	2	3	2	5	6	5	2	2	6	7	5	2	6	4	3	6	2	2	6	5	5	6	6	6	3	6	5	4	2	6	4	5	6	2	2	5	4	7	6	6	2	2	2	6	7	6	2	6	2	2	5	3];
% 进行样条插�?�计�?
spline_interp1 = csape(x_1known, y_1known, 'complete', 'variational');

% 定义插�?�的目标�?
x_1interp = 0:0.01:100;

% 计算插�?�结�?
y_1interp = ppval(spline_interp1, x_1interp);

% 定义已知数据�?
x_2known = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100];
y_2known=[-2	-1	-2	-1	-5	-1	1	0	-7	-4	-4	-1	-2	0	-4	0	-1	-2	-6	-4	-3	-5	-2	-3	-6	0	-7	-5	-5	-3	-7	-4	-7	-6	0	-5	-2	1	-4	-1	-1	-4	-2	-7	1	-6	-5	0	-3	-1	-4	-5	-7	-1	-4	-3	-2	-7	-5	-1	-1	-6	-6	-7	-7	-4	-2	-1	-3	-7	-2	-6	-6	-7	-6	-6	-6	-5	-5	-6	-5	1	-1	-2	-6	-6	-7	1	-1	-2	-5	-6	-2	1	-6	-5	-4	-7	-1	-4	1];
% 进行样条插�?�计�?
spline_interp2 = csape(x_2known, y_2known, 'complete', 'variational');

% 定义插�?�的目标�?
x_2interp = 0:0.01:100;

% 计算插�?�结�?
y_2interp = ppval(spline_interp2, x_2interp);

% 定义已知数据�?
x_3known = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100];
y_3known=[-40 -73 -58 -36 -3 -2 -68 -2 -3 -41 -16 -69 -46 -6 -5 -3 -7 -4 -8 -5 -11 -13   -12   -15   -13   -15   -10   -10   -14   -11   -14 -17	-15	-17	-19	-19	-18	-18	-19	-17	-16 -23	-23	-26	-26	-24	-24	-21	-19	-19	-23 -30	-25	-25	-26	-26	-31	-26	-28	-30	-31 -30	-35	-35	-31	-29	-32	-31	-35	-29	-32 -36	-41	-35	-36	-41	-37	-39	-37	-38	-38 -45	-44	-46	-39	-41	-39	-45	-39	-40	-42 -47	-48	-45	-48	-49	-45	-45	-45	-46	-47];
% 进行样条插�?�计�?
spline_interp3 = csape(x_3known, y_3known, 'complete', 'variational');

% 定义插�?�的目标�?
x_3interp = 0:0.01:100;

% 计算插�?�结�?
y_3interp = ppval(spline_interp3, x_3interp);
% % 绘制连续不规则函�?
% figure;
% plot(x_known, y_known, 'o', x_1interp, y_interp);
% legend('已知数据�?', '插�?�结�?');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tfai2=(fai2)';
tic
u_1=[];u_2=[];u_11=[];u_22=[];x_1=[];x_2=[];x_11=[];x_22=[];x_33=[];y_11=[];y_22=[];y_33=[];
for i=1:length(t)
    u_1=[u_1 uu(1,i)];
    u_2=[u_2 uu(2,i)];
    u_11=[u_11 uu(1,i)];
    u_22=[u_22 uu(2,i)];
    x_1=[x_1 x1(i)];
    x_2=[x_2 x2(i)];
    x_11=[x_11 x1(i)];
    x_22=[x_22 x2(i)];
    x_33=[x_33 x3(i)];
    y_11=[y_11 x1(i)];
    y_22=[y_22 x2(i)];
    y_33=[y_33 x3(i)];
end
for i=1:length(t)
   if (t(i)>=20|t(i)<=40) & (t(i)>=30|t(i)<=40) & (t(i)>=50|t(i)<=60) & (t(i)>=70|t(i)<=80)
    u_11(i)=0.1*u_1(i);%0.5
else
    u_11(i)=u_1(i);
end 

if (t(i)>=20|t(i)>=40) & (t(i)>=30|t(i)<=40) & (t(i)>=50|t(i)<=60) & (t(i)>=70|t(i)<=80)
    u_22(i)=0.1*u_2(i);
else
    u_22(i)=u_2(i);
end  
end

% if t<=10
%     u1=u_1;
% else
%     u1=50;
% end 
% 
% if   t<=10
%     u2=u_1;
% else
%     u2=15;
% end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i=6000:6300
%     u_1(i)=10*u_1(i);
% end
% uu_1=[];uu_2=[];mm=0;
% for i=1:length(t)
%    mm=mm+1;
% 
% if mm<=10/0.01
%     uu_1=[uu_1 u_1(i)];
% else
%     uu_1=[uu_1 5];
% end
% end
% 
% nn=0;
% for i=1:length(t)
%    nn=nn+1;
% if nn<=15/0.01
%     uu_2=[uu_2 u_2(i)];
% else
%     uu_2=[uu_2 0.5*u_2(i)];
% end
% end
% er1_1=0.2*ones(20001,1);
% er1_2=-0.2*ones(20001,1);
% er2_1=0.3*ones(20001,1);
% er2_2=-0.3*ones(20001,1);
% er3_1=0.2*ones(20001,1);
% er3_2=-0.2*ones(20001,1);
err_1=x1-yd1;
err_2=x2-yd2;
err_3=x3-yd3;

 err_11=[];err_22=[];err_33=[];
 for i=1:length(t)
     err_11=[err_11 err_1(i)];
     err_22=[err_22 err_2(i)];
 end
for i=length(t)
    x_11(i)=x1(i);
    x_22(i)=x2(i);
    x_33(i)=x3(i);
end
for i=1:1500
    y_11(i)=x1(i);
    y_22(i)=x2(i);
    y_33(i)=x3(i);
end
for i=1501:10001
     y_11(i)=y_1interp(i);
    y_22(i)=y_2interp(i);
    y_33(i)=y_3interp(i);
end
err_11=y_11'-yd1;%+0.25*cos(0.2*t)
err_22=y_22'-yd2;%-0.2*cos(0.1*t)
err_33=y_33'-yd3;
    for  i=1501:10001
    x_11(i)=NaN;
    x_22(i)=NaN;
    x_33(i)=NaN;
    end
    for i=1500
    x_11(i)=8;
    x_22(i)=8;
    x_33(i)=20;
    end

% for i=2000:2050
%     x_1(i)=cos(0.007*i)-0.8;
%     x_2(i)=cos(0.007*i)-1;
%     err_11=20*sin(i);
%     err_22=20*cos(i);
% end
% for i=6000:6050
%    x_1(i)=cos(0.007*i)-1.4;
%    x_2(i)=cos(0.007*i)+0.22;
%    err_11=20*sin(i);
%    err_22=20*cos(i);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(t)
xx(i)=x4(i)*sin(x3(i));
yy(i)=x4(i)*cos(x3(i));
R(i)=(xx(i))^2+(yy(i))^2;
xxd_1(i)=yd1(i);
yyd_2(i)=yd2(i);
Rr(i)=(xxd_1(i))^2+(yyd_2(i))^2;
end
rou1=4*exp(-0.21*t)+0.34;
rou2=-4*exp(-0.21*t)-0.34;
ze=((err_1).^2+(err_2).^2).^0.5;
ze1=((err_11).^2+(err_22).^2).^0.5;
toc
disp(['运行时间4: ',num2str(toc)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ����һЩʾ������
actual = x1;
predicted = yd1;

% ���ú�������MAE
mae_value = calculateMAE(actual, predicted);
mse_value = calculateMSE(actual, predicted);

% ��ʾMAE���
disp(['MAE: ', num2str(mae_value)]);

% ��ʾMAE���
disp(['MSE: ', num2str(mse_value)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
subplot(2,2,1)
plot(t,yd1,'r',t,x1,'--g','LineWidth',1);hold on
h=legend('$\varsigma_{d}$','$\varsigma$','Location','Northeast');
 axis([0  100  -8 8]);
set(h,'interpreter','latex')
xlabel('$Time[s]$','interpreter','latex');
ylabel('Displacement','interpreter','latex');
% grid on;
subplot(2,2,2)
plot(t,yd2,'r',t,x2,'--g','LineWidth',1);hold on
h=legend('$\mu_{d}$','$\mu$','Location','Northeast');
axis([0  100  -8 8]);
set(h,'interpreter','latex')
xlabel('$Time[s]$','interpreter','latex');
ylabel('Displacement','interpreter','latex');
% grid on;
subplot(2,2,3)
plot(t,yd3,'r',t,x3,'--g','LineWidth',1);hold on
h=legend('$\psi_{d}$','$\psi$','Location','Northeast');
axis([0 100  -40 10]);
set(h,'interpreter','latex')
xlabel('$Time[s]$','interpreter','latex');
ylabel('angle','interpreter','latex');
% grid on;
subplot(2,2,4)
plot(x1,x2,'--g',yd1,yd2,'r','LineWidth',1);hold on
axis([-4  4  -3 3]);
h=legend('$actual$','$reference$','Location','Northeast');
set(h,'interpreter','latex')
xlabel('Displacement','interpreter','latex');
ylabel('Displacement','interpreter','latex');
% grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(2)
% plot(t,yd2,'r',t,x2,'--g','LineWidth',2);hold on
% h=legend('$\mu_{d}$','$\mu$','Location','Northeast');
% axis([0  100  -8 8]);
% set(h,'interpreter','latex')
% xlabel('$Time[s]$','interpreter','latex');
% ylabel('Displacement','interpreter','latex');
% % grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(3)
% plot(t,yd3,'r',t,x3,'--g','LineWidth',2);hold on
% h=legend('$\psi_{d}$','$\psi$','Location','Northeast');
% axis([0 100  -40 10]);
% set(h,'interpreter','latex')
% xlabel('$Time[s]$','interpreter','latex');
% ylabel('angle','interpreter','latex');
% % grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
% subplot(3,1,1)
plot(t,ze,'g',t,err_3,'m',t,rou1,'b',t,rou2,'b','LineWidth',1);
h=legend('${z_{e}}$','${\psi_{e}}$','��','Location','Northeast');
set(h,'interpreter','latex')
xlabel('$Time[s]$','interpreter','latex');
ylabel('Displacement','interpreter','latex');
axis([0 100 -5 5]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(7)
% plot(x1,x2,'--g',yd1,yd2,'r','LineWidth',2);hold on
% axis([-4  4  -3 3]);
% h=legend('$actual$','$reference$','Location','Northeast');
% set(h,'interpreter','latex')
% xlabel('Displacement','interpreter','latex');
% ylabel('Displacement','interpreter','latex');
% % grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
subplot(2,1,1)
plot(t,u_1,'m','LineWidth',1);%t,uc1,'--b',
h=legend('$u_{c,h1}$','Location','Northeast');%'${u1 FTC}$'
set(h,'interpreter','latex')
xlabel('$Time[s]$','interpreter','latex');
ylabel('Value','interpreter','latex');
axis([0.5 100 -30 30]);
subplot(2,1,2)
plot(t,u_2,'m','LineWidth',1);%t,uc2,'--b',
h=legend('$u_{c,h2}$','Location','Northeast');%'${u2 FTC}$',
set(h,'interpreter','latex')
xlabel('$Time[s]$','interpreter','latex');
ylabel('Value','interpreter','latex');
axis([0.5 100 -30 30]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%