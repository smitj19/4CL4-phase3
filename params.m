clear;

x0_up_1 = [0; pi; 0; 0];
x0_up_2 = [0; 0.95*pi; 0; 0];
x0_up_3 = [0.1; pi; 0; 0];
x0_up_4 = [0; pi; 0; 0];

x0_down_1 = [0; 0; 0; 0];
x0_down_2 = [10; 10; 0; 0];
x0_down_3 = [0; 0.1; 0; 0];
x0_down_4 = [-0.1; 0; 0; 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X0= x0_up_1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q20 = X0(2); 

kp=4;
Ts = 0.001;
m1=0.095; 
m2=0.024;
L1= 0.1;
L2= 0.129;
l1 =L1/2;
l2 =L2/2;
g= 9.81;
J1 = (1/12)*m1*L1^2;
J2 = (1/12)*m2*L2^2;
Rm=8.4;
km= 0.042;
b1= 0.001;
b2 = 1e-5;
Jb1 = J1+m1*l1^2+m2*L1^2;
Jb2 = J2+m2*l2^2;
%%
%--------PHASE 3 STUFF----------
P = load('ph.mat').Ph;
M = [P(1) P(6)*cos(q20); P(6)*cos(q20) P(2)];
F = [P(3) 0; 0 P(4)];
K = [0 0; 0 P(5)*cos(q20)];
T = [1; 0];

a11 = [0 0; 0 0];
a12 = [1 0; 0 1];
a21 = -inv(M)*K;
a22 = -inv(M)*F;
A = [a11 a12; a21 a22];
b1 = [0; 0];
b2 = inv(M)*T;
B = [b1; b2];
sys_poles = eig(A);
disp(sys_poles) %TASK 1

%If real -ve poles -> stable, if real +ve poles -> unstable
Pc = zeros(size(A,1), size(A,2));
for N = 1:size(A,1)
Pc(N,:) = A^(N-1)*B;
end
disp(rank(Pc)) %TASK 2

%Sys is controllable if rank = 4

ps_1 = [-10+10*1i -10-10*1i -15 -18];
ps_2 = 0.5*ps_1;
ps_3 = 1.3*ps_1;

%%%%%%MODIFY FOR TASK 5%%%%%%%
ps = ps_3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kf = place(A, B, ps);
sys_poles2 = eig(A - B*Kf);
disp(sys_poles2) %Confirm stability