function [cost,flag]=cost_zeta(AC,K,individual,zeta)

A11 = individual(1);
A12 = individual(2);
A22 = individual(3);
C12 = individual(4);

W11 = 0.005;
W22 = 0.005;
V11 = 0.001;

%%%%% initialize
data_amount=400;
state=2;
out=1;
in=2;

AA=[AC(1)  AC(2); -0.1  AC(3)];
CC=[0 AC(4)];
BB = eye(2);
% TTo1=inv([CC*AA;  CC])*[1; 0];
% TTo=[AA*TTo1  TTo1];
% invTTo=inv(TTo);
% AAo=invTTo*AA*TTo;
% BBo=invTTo*eye(2);
% CCo=CC*TTo;

%%%%% convert the system(X) into the observe form(Xo)
A=[A11 A12; -0.1  A22];
B=eye(2);
C=[0 C12];

% To1=inv([C*A;  C])*[1; 0];
% To=[A*To1  To1];
% invTo=inv(To);
% Ao=invTo*A*To;
% Bo=invTo*B;
% Co=C*To;

%%%%% self tunning
% A1=-Ao(1,1);
% A2=-Ao(2,1);
% B1=Bo(1,:);
% B2=Bo(2,:);
A1=-A(1,1);
A2=-A(2,1);
B1=B(1,:);
B2=B(2,:);

%%% test the Ko
test_D=[-(A1+K(1,1))  eye(1)
        -(A2+K(2,1))  zeros(1)];
s=1;
if (isinf(-(A1+K(1,1))) || isnan(-(A1+K(1,1))) || isinf(-(A2+K(2,1))) || isnan(-(A2+K(2,1))))
    flag=1;
    cost=0;
    return;
end
if (max(abs(eig(test_D))))>=s
   flag=1;
   cost=0;
   return;
end;

Cov_W=diag([W11  W22]);
Cov_V=diag(V11);

sprev = rng;
rng(0, 'twister');
% rstream = RandStream('mcg16807', 'Seed',4262);
W=sqrt(Cov_W)*randn(state,data_amount);
% W=invTo*sqrt(Cov_W)*invTo'*randn(rstream,state,data_amount);
% W=invTo*sqrt(Cov_W)*invTo'*randn(state,data_amount);
rng(3571, 'twister');
% rstream = RandStream('mcg16807', 'Seed',91320);
V=sqrt(Cov_V)*randn(out,data_amount);
rng(sprev);
% V=sqrt(Cov_V)*randn(out,data_amount);


X			= zeros(state, data_amount);
X(:,1)		= [0.1;-0.1];
X_h			= zeros(state, data_amount);
% Xo(:,1)=zeros(state,1);
% Xo_h(:,1)=zeros(state,1);
Y			= zeros(1, data_amount);
Y(1)		= C * X(:, 1) + V(:, 1);
U			= zeros(in, data_amount);
e			= zeros(1, data_amount);
e(1)		= Y(1);
sum=(X(:,1)-X_h(:,1))'*(X(:,1)-X_h(:,1));
for k=2:data_amount
   X(:,k)=[-A1  eye(out); -A2  zeros(out)]*X(:,k-1)+[B1; B2]*U(:,k-1)+W(:,k-1);
   Y(:,k)=C*X(:,k)+V(:,k);
   
   X_h(:,k) = AA * X_h(:,k-1) + BB * U(:, k-1) + K * zeta * e(:,k-1);
   e(:,k)=Y(:,k)-CC*X_h(:,k);
   U(:,k)=0;
   sum=sum+(X(:,k)-X_h(:,k))'*(X(:,k)-X_h(:,k));
end;
cost=sum/data_amount;
flag=0;