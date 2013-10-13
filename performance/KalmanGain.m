function [flag,K]=KalmanGain(input_val)
% Calculate Kalman Gain

% parameter
A11 = input_val(1);
A12 = input_val(2);
A22 = input_val(3);
C12 = input_val(4);
W11 = 0.005;
W22 = 0.005;
V11 = 0.001;

data_amount=400;

state=2;
out=1;
in=2;
r0=0.95; 
r=0.91;  
P=eye(state*1)*10^(-1); %#ok<*MFAMB>

% convert the system(X) into the observe form(Xo)
A=[A11 A12;-0.1 A22];
B=eye(2);
C=[0 C12];

% To1=([C*A;C])\[1;0];
% To=[A*To1 To1];
% invTo=inv(To);
% Ao=invTo*A*To;
% Bo=invTo*B;
% Co=C*To;

Cov_W=diag([W11 W22]);
Cov_V=diag(V11);

% self tunning
A1=-A(1,1);
A2=-A(2,1);
B1=B(1,:);
B2=B(2,:);
% A1=-Ao(1,1);
% A2=-Ao(2,1);
% B1=Bo(1,:);
% B2=Bo(2,:);

D1=0;
D2=0;
% d(:,1)=[D1;D2];

sprev = rng;
rng(0, 'twister');
% rstream = RandStream('mcg16807', 'Seed',4262);
W=sqrt(Cov_W)*randn(state,data_amount);
% W=invTo*sqrt(Cov_W)*invTo'*randn(rstream,state,data_amount);
% W=invTo*sqrt(Cov_W)*invTo'*randn(state,data_amount);
rng(3571, 'twister');
% rstream = RandStream('mcg16807', 'Seed',91320);
V=sqrt(Cov_V)*randn(out, data_amount);
rng(sprev);
% V=sqrt(Cov_V)*randn(out, data_amount);

X			= zeros(state, data_amount-1);
X_h			= zeros(state, data_amount-1);
Y			= zeros(1, data_amount-1);
Y(1)		= C * X(:, 1)+V(:, 1);
epsilon		= zeros(1, data_amount-1);
epsilon(1)	= Y(1);
e			= zeros(1, data_amount-1);
e(1)		= Y(1);
U			= zeros(in, data_amount-1);
cita_k=[D1 D2]';

for k=2:data_amount-1
   %k
   r=r0*r + (1-r0);
   if k==2
      fi=[epsilon(:,k-1); zeros(out,1)];
   else
      fi=[epsilon(:,k-1); epsilon(:,k-2)];
   end;
   X(:,k)=[-A1  eye(out); -A2  zeros(out)]*X(:,k-1)+[B1; B2]*U(:,k-1)+W(:,k-1);
   Y(:,k)=C*X(:,k)+V(:,k);
   if k==2
      epsilon(:,k)=Y(:,k)+A1*Y(:,k-1)-B1*U(:,k-1)-cita_k'*fi;
   else
      epsilon(:,k)=Y(:,k)+A1*Y(:,k-1)+A2*Y(:,k-2)-B1*U(:,k-1)-B2*U(:,k-2)-cita_k'*fi;   
   end;

   cita=cita_k;
   G=(P*fi)/(r+fi'*P*fi);
   cita_k=cita+G*epsilon(1,k);
   P=(P-((P*fi)*(P*fi)'/(r+fi'*P*fi)))/r;
   
   tmp=cita_k';
   D1=tmp(1,1);
   D2=tmp(1,2);
%    d(:,k)=[D1;D2];
   
   K=[D1-A1; D2-A2];
   
   X_h(:,k)=A*X_h(:,k-1)+B*U(:,k-1)+K*e(:,k-1);
   e(:,k)=Y(:,k)-C*X_h(:,k);
   U(:,k)=0;
end;

% test the Ko
test_D=[-(A1+K(1,1)) eye(1);-(A2+K(2,1)) zeros(1)];
s=1;
flag=0;
if (max(abs(eig(test_D))))>=s
   flag=1;
end;