clc;clear;close all;

%An SDP that satisfies strict complementarity and locally QG
%The plot wants to show f(y) - f* satisfies QG
%2 x 2 example
%unique dual solution y* = [0,0]


%Problem data
n   = 2;
rho = 4; %exact penalty parameter
b   = ones(n,1);
C = [1,-1;
     -1,1];
c  = reshape(C,[],1);
A1      = zeros(n,n);
A1(1,1) = 1;
A2      = zeros(n,n);
A2(2,2) = 1;


ystar   = [0;0];
resolution = 500;
dy      = 0.9; %interval
y1      = linspace(ystar(1) - dy,ystar(1)+dy,resolution);
y2      = y1;

[Y1,Y2] = meshgrid(y1,y2);
obj_nonlinear = zeros(length(y1),length(y2));
obj_linear = zeros(length(y1),length(y2));
dist    = zeros(length(y1),length(y2));
%grad    = zeros(length(y1),length(y2));
%EIGS    = zeros(length(y1),length(y2));
for i =1:length(y1)
    for j =1:length(y2)
        lambda_max = max([0;eig(-C+A1*y1(i)+A2*y2(j))]);
        if lambda_max == 0
            obj_linear(i,j)    = -b(1)*y1(i)-b(2)*y2(j)+rho*lambda_max + b.'*ystar;
            obj_nonlinear(i,j) = nan;
        else
            obj_nonlinear(i,j) = -b(1)*y1(i)-b(2)*y2(j)+rho*lambda_max + b.'*ystar;
            obj_linear(i,j)    = nan;
        end
        dist(i,j) = sqrt((y1(i)-ystar(1))^2+(y2(j)-ystar(2))^2);
        % [V,D]     = eig(C-A1*y1(i)-A2*y2(j));
        % EIGS(i,j) = min(diag(D));
        % if EIGS(i,j) >= 0
        %     grad(i,j) = 2;
        % elseif EIGS(i,j) < 0
        %     %check the multiplicity
        %     repeat = sum(D == EIGS(i,j),'all');
        %     if repeat>1
        %         disp('wait');
        %     end
        %     [d,ind] = sort(diag(D));
        %     Vs = V(:,ind);
        %     v = Vs(:,1);
        %     grad(i,j) = norm((-b-rho*v.^2));
        % end
    end
end


f1 = figure;
%firstax = axes (f1, 'FontSize', 16); 
h1 = surf(Y1,Y2,obj_linear,'EdgeColor', 'none','FaceColor','#EDB120','FaceAlpha',0.5);

% 
hold on;
h2 = surf(Y1,Y2,obj_nonlinear,'EdgeColor', 'none','FaceColor','#4DBEEE','FaceAlpha',0.5);


%optimal solution 
plot3(0,0,0,'-o','Color','none','MarkerSize',6,'MarkerFaceColor','r');

%Draw the boundary where eig = 0
d1 = linspace(-dy,0);
d2 = d1./(d1-1);
QG_dir = zeros(length(d1),1);
for i = 1:length(d1)
    QG_dir(i) = -b(1)*d1(i)-b(2)*d2(i)+rho*max([0;eig(-C+A1*d1(i)+A2*d2(i))]) + b.'*ystar;
end
plot3(d1,d2,QG_dir,'Color','m','LineWidth', 2);


d2 = linspace(-dy,0);
d1 = d2./(d2-1);
QG_dir = zeros(length(d1),1);
for i = 1:length(d1)
    QG_dir(i) = -b(1)*d1(i)-b(2)*d2(i)+rho*max([0;eig(-C+A1*d1(i)+A2*d2(i))]) + b.'*ystar;
end
plot3(d1,d2,QG_dir,'Color','m','LineWidth', 2);

h3 = surf(Y1,Y2,0.3*dist.^2,'EdgeColor', 'none','FaceColor','#77AC30','FaceAlpha',0.7);

hd =  legend([h1,h2,h3],'$-b^{\top}y$','$-b^{\top}y + \rho\max\{ \lambda_{\max}(\mathcal{A}^*y - C),0\}$','$\kappa \cdot\mathrm{Dist}^2(y,S)$','interpreter','latex',...
    'Location','none','Position',[0.07,0.93,0.87,0.05],'Box','off','FontSize', 12, 'NumColumns',3);

xlabel('$y_1$','interpreter','latex');
ylabel('$y_2$','interpreter','latex');
grid on


view([47 30]);



width  = 5.2;     % Width in inches
height = 4;    % Height in inches
set(gcf, 'Position', [300 100  width*100, height*100]); %<- Set size
set(gca, 'FontSize', 12); %<- Set properties

%print(gcf,'ExactPenalty_QG_3D.eps','-depsc2','-r300');


%%
%Sectional view

figure();
%Draw the boundary where eig = 0
d1 = linspace(-dy,0);
d2 = d1./(d1-1);
QG_dir = zeros(length(d1),1);
for i = 1:length(d1)
    QG_dir(i) = -b(1)*d1(i)-b(2)*d2(i)+rho*max([0;eig(-C+A1*d1(i)+A2*d2(i))]) + b.'*ystar;
end
plot(d1,QG_dir,'Color','m','LineWidth', 2);


hold on
d2 = linspace(-dy,0);
d1 = d1./(d2-1);
QG_dir = zeros(length(d1),1);
for i = 1:length(d1)
    QG_dir(i) = -b(1)*d1(i)-b(2)*d2(i)+rho*max([0;eig(-C+A1*d1(i)+A2*d2(i))]) + b.'*ystar;
end
plot(d1,QG_dir,'Color','m','LineWidth', 2);
xlabel('$y_1$','interpreter','latex');
ylabel('$f(y)$','interpreter','latex');
width  = 4;
height = 3.5;
set(gcf, 'Position', [300 100  width*100, height*100]); %<- Set size
set(gca, 'FontSize', 12); %<- Set properties
plot(0,0,'-o','Color','none','MarkerSize',6,'MarkerFaceColor','r');%#D9FFFF
grid on;
%print(gcf,'ExactPenalty_QG_Sectional.eps','-depsc2','-r300');
