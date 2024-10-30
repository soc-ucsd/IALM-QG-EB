clc; clear; close all
yalmip('clear')
%ALM for Maxcut
%Package requirement: Yalmip and Mosek

addpath("Packages\");

datapath = 'data\LASSO_';
savepath = 'results/LASSO_';
name     = {'1','2','3'}; %lp_afiro %lp_blend
idx      =3;


load([datapath,name{idx},'.mat']);

At  = full(model.A');
b   = full(model.b);
c   = full(model.C);
K   = model.K;
[xx,yy,info] = sedumi(At,b,c,K);
Optimal.Cost = c.'*xx;


m        = height(At);
n        = width(At);
nf       = K.f;
nl       = K.l;
nq       = K.q;
r        = 10;               %augmetned term

x        = sdpvar(n,1); %Yalmip variable
u        = sdpvar(nl+nq,1);


Xtrue                     = xx;

ytrue                    = yy;
Ztrue                    = c -At.'*ytrue;

%nonnegative and quadratic part
A2   = At(:,nf+1:end);
c2   = c(nf+1:end);

%initialization for ALM
yk       = zeros(m,1);
z        = c2 - A2.'*yk;
zk       = z;
C        = c;
% zk       = reshape(z,n,n);
% C        = reshape(z,n,n);


Max_iter     = 20;
cost         = [];
Dcost        = [];
xstar        = [];
Affinefeasi  = [];
DAffinefeasi = [];
Conefeasi    = [];
DualGap      = [];
Dist        = [];
DDist        = [];
normb        = norm(b);
normC        = norm(c,'fro');
normXtrue    = norm(Xtrue,'fro');
normytrue    = norm(ytrue,'fro');
normZtrue    = norm(Ztrue,'fro');


%main algorithm here

for k = 1:Max_iter
    %solve the subproblem
    bAx       = b-At*x;
    %xu        = x-u; %trace(C*x)
    lk        =  c(:).'*x(:) + 1/(2*r)*norm(r*x(nf+1:end)-zk-u,'fro')^2 + 1/(2*r)*norm(yk+r*bAx)^2;
    cons      = [u(1:nl)>=0,u(nl+1)>=norm(u(nl+2:end))];

    ops = sdpsettings('solver','mosek');
    optimize(cons,lk);

    xk        = value(x);
    cost(k)   = c.'*x(:);
    xstar{k}  = xk;
    
    %dual update
    tempzk    = zk - r*xk(nf+1:end);

    %tempzk_f  = tempzk(1:nf);
    tempzk_l  = tempzk(1:nl);
    tempzk_q  = tempzk(nl+1:end);
    
    tempzk_l(tempzk_l<= 0) = 0;

    tempzk_q  = projSOCP(tempzk_q);
    zk        = [tempzk_l;tempzk_q];


    yk        = yk+r*value(bAx);
    
    Ztest     = c2 - A2.'*yk;
    
    %residual
    Dcost        = [Dcost,b.'*yk]; 
    Affinefeasi  = [Affinefeasi,norm(b-At*vec(xk))/(1+normb)];
    DAffinefeasi = [DAffinefeasi,norm(Ztest-zk)/(1+normC)];

    xk_f         = xk(1:nf);  
    xk_l         = xk(nf+1:nf+nl);
    xk_q         = xk(nf+nl+1:end);
    residual     = sqrt((norm(xk_l(xk_l<0)))^2+(norm(xk_q -projSOCP(xk_q)))^2);
    Conefeasi    = [Conefeasi,residual/(1+norm(xk))];
    Dist         = [Dist,norm(xk-Xtrue,'fro')/(1+normXtrue)];
    %DDist        = [DDist,sqrt(norm(yk-ytrue)^2+norm(zk-Ztrue,'fro')^2)/(1+normytrue+normZtrue)];
end


Out.PCostgap     = abs(cost-Optimal.Cost)/(1+abs(Optimal.Cost));
Out.DCostgap     = abs(Dcost-Optimal.Cost)/(1+abs(Optimal.Cost));

%
Out.Conefeasi    = Conefeasi;
Out.PCost        = cost;
Out.DCost        = Dcost;
Out.Affinefeasi  = Affinefeasi;
Out.DAffinefeasi = DAffinefeasi;
Out.Dist         = Dist;
Out.DDist        = DDist;
Out.OptimalCost  = Optimal.Cost;


semilogy(Out.PCostgap);
hold on
semilogy(Out.Conefeasi);

save([savepath,name{idx},'_r',int2str(r),'_result'],'Out');


function y = projSOCP(x)
%     SOC = true;
%     threshold = x(1) - norm(x(2:end));
    normq = norm(x(2:end));
    if x(1) >= normq %x is already socp
        y = x;
    elseif -x(1)>= normq
        y = 0;
    else
        y = (x(1)+normq)/(2*normq)*[normq;x(2:end)];
    end
end
