clc; clear; close all
yalmip('clear')
%ALM for Maxcut
%Package requirement: Yalmip and Mosek

addpath("Packages\");
%addpath("Examples\")

datapath = 'data\';
savepath = 'results/';
name     = {'G1_n20'};
idx        = 1;


load([datapath,name{idx},'.mat']);


[info,mosektime]=SolveMosek(At,b,c,K);
Optimal.Cost = info.sol.itr.pobjval;


m        = height(At);
n        = sqrt(width(At));
r        = 1;               %augmetned term

x        = sdpvar(n); %Yalmip variable
u        = sdpvar(n);


%Extract the true solution from Mosek
Xtrue                    = zeros(n);
[IndSym,IndDiag,IndOffDiag,ShrinkIndDiag,ShrinkIndOffDiag,IndOffDiagCounter] = SymmetricIndices(n,false);
Xtrue(IndSym)            = info.sol.itr.barx;
Xtrue(IndOffDiagCounter) = Xtrue(IndOffDiag);

ytrue                    = info.sol.itr.y;
Ztrue                    = reshape(c -At.'*ytrue,n,n);


%initialization for ALM
yk       = zeros(m,1);
z        = c - At.'*yk;
zk       = reshape(z,n,n);
C        = reshape(z,n,n);


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
normC        = norm(C,'fro');
normXtrue    = norm(Xtrue,'fro');
normytrue    = norm(ytrue,'fro');
normZtrue    = norm(Ztrue,'fro');


%main algorithm here

for k = 1:Max_iter
    %solve the subproblem
    bAx       = b-At*vec(x);
    xu        = x-u; %trace(C*x)
    lk        =  c(:).'*x(:) + 1/(2*r)*norm(r*x-zk-u,'fro')^2 + 1/(2*r)*norm(yk+r*bAx)^2;
    cons      = [u>=0];

    ops = sdpsettings('solver','mosek');
    optimize(cons,lk);

    xk        = value(x);
    cost(k)   = c.'*x(:);
    xstar{k}  = xk;
    
    %dual update
    tempzk    = zk - r*xk;
    [V,D]     = eig(tempzk); 
    D(D<=0)   = 0;
    zk        = V*D*V';
    zk        = (zk+zk')/2;

    yk        = yk+r*value(bAx);
    
    Ztest     = reshape(c - At.'*yk,n,n);
    
    %residual
    Dcost        = [Dcost,b.'*yk]; 
    Affinefeasi  = [Affinefeasi,norm(b-At*vec(xk))/(1+normb)];
    DAffinefeasi = [DAffinefeasi,norm(Ztest-zk)/(1+normC)];
    %DualGap      = [DualGap,abs(cost(k) - b'*yk)];
    [V,D]        = eig(xk);
    Conefeasi    = [Conefeasi,norm(D(D<0))/(1+norm(xk))];
    Dist         = [Dist,norm(xk-Xtrue,'fro')/(1+normXtrue)];
    DDist        = [DDist,sqrt(norm(yk-ytrue)^2+norm(zk-Ztrue,'fro')^2)/(1+normytrue+normZtrue)];
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
Out.OptimalCost  = info.sol.itr.pobjval;


semilogy(Out.PCostgap);
hold on
semilogy(Out.Conefeasi);

%save([savepath,name{idx},'_r',int2str(r),'_result_new'],'Out');
