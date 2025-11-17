P =Ph; %= [0.0341; 0.0150; 0.0618; 0.0012; 1.8664; 0.0158];

p1 = P(1);  % J1_bar
p2 = P(2);  % J2_bar
p3 = P(3);  % b1 + kt^2/R
p4 = P(4);  % b2
p5 = P(5);  % m2*g*l2
p6 = P(6);  % m2*L1

% Choose equilibrium: q20 = 0 for DOWN, q20 = pi for UP
q20 = pi;                     % <-- set to 0 or pi
c   = cos(q20);

% M, F, K from the lab handout
M = [p1,    p6*c;
     p6*c,  p2   ];

F = [p3, 0;
     0,  p4];

K = [0,     0;
     0,  p5*c];

detM = p1*p2 - (p6*c)^2;
Minv = (1/detM) * [ p2,     -p6*c;
                   -p6*c,    p1   ];

A = [ 0 0  1 0;
      0 0  0 1;
     -(Minv*K),  -(Minv*F) ];

B = [ 0;
      0;
      Minv*[1;0] ];

% --- Controllability matrix Pc = [B, A*B, A^2*B, A^3*B] ---
AB   = A*B;
A2B  = A*AB;
A3B  = A*A2B;
Pc   = [B, AB, A2B, A3B];    % 4×4 for this SISO 4-state model
%%
% Parameters
P =Ph; %= [0.0341; 0.0150; 0.0618; 0.0012; 1.8664; 0.0158];
p1=P(1); p2=P(2); p3=P(3); p4=P(4); p5=P(5); p6=P(6);


c = cos(0);
M = [p1, p6*c; p6*c, p2];
F = [p3, 0;    0,    p4];
K = [0,  0;    0, p5*c];

Minv = inv(M);
A_down = [zeros(2), eye(2); -Minv*K, -Minv*F];
B_down = [zeros(2,1); Minv*[1;0]];

% Desired poles
ps = [-10+10j, -10-10j, -15, -18];


% Gain and verification
K_down = place(A_down, B_down, ps);
eig_cl_down = eig(A_down - B_down*K_down)   % should equal poles (order may vary) 

% ---------- UP equilibrium (q20 = pi) ----------
c = cos(pi);
M = [p1, p6*c; p6*c, p2];
F = [p3, 0;    0,    p4];
K = [0,  0;    0, p5*c];

Minv = inv(M);
A_up = [zeros(2), eye(2); -Minv*K, -Minv*F];
B_up = [zeros(2,1); Minv*[1;0]];

K_up = place(A_up, B_up, ps);
eig_cl_up = eig(A_up - B_up*K_up) % should equal the poles given
