close all;clc; clear;

%description

%parameter setting 
width  = 10;     % Width in inches
height = 3;    % Height in inches
alw    = 0.75;    % AxesLineWidth
fsz    = 12;      % Fontsize
lw     = 2;      % LineWidth
msz    = 8;       % MarkerSize

colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.4940 0.1840 0.5560],...
          [0.9290 0.6940 0.1250],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],...
          [1 0 0],[0 0 1],[0 1 0]}; 


datapath = "G1_n20";
name     = {'1','2','3'};

range = 1:20;

for idx = 1:1
    for R = [1,5,10]
        load(datapath+"_r"+int2str(R)+"_result_new");
        subplot(1,3,1);
        %semilogy(Out.Dist(range),'-o','LineWidth',1);
        semilogy(Out.PCostgap(range),'-o','LineWidth',1);
        hold on 
        %
        
        subplot(1,3,2);
        semilogy(Out.DCostgap(range),'-o','LineWidth',1);
        hold on 
    
        subplot(1,3,3);
        Feasibility = max([Out.Affinefeasi;Out.DAffinefeasi;Out.Conefeasi;abs(Out.PCost - Out.DCost)/(1+abs(Out.OptimalCost))]);
        semilogy(Feasibility(range),'-o','LineWidth',1);
        hold on 
    end
end

set(gcf, 'Position', [300 100  width*100, height*100]); %<- Set size

subplot(1,3,1);
xlim([range(1),range(end)]);
xlabel('Iteration','interpreter','latex');
set(gca,'TickLabelInterpreter','latex' ,'FontSize', fsz, 'LineWidth', alw); %<- Set propertiesylabel('$\epsilon_1$','interpreter','latex', 'LineWidth', alw,'FontSize', 14)

set(gca, 'Position', [0.1 0.3 0.225 0.6]); %<- Set properties
title('Primal cost value gap','interpreter','latex','FontSize', fsz);
hold on 

subplot(1,3,2);
xlim([range(1),range(end)]);
xlabel('Iteration','interpreter','latex');
set(gca,'TickLabelInterpreter','latex' ,'FontSize', fsz, 'LineWidth', alw); %<- Set properties
ylabel('$\epsilon_2$','interpreter','latex', 'LineWidth', alw,'FontSize', 14)

set(gca, 'Position', [0.42 0.3 0.225 0.6]); %<- Set properties
title('Dual cost value gap','interpreter','latex','FontSize', fsz);
hold on 


subplot(1,3,3);
xlim([range(1),range(end)]);
xlabel('Iteration','interpreter','latex');
set(gca,'TickLabelInterpreter','latex' ,'FontSize', fsz, 'LineWidth', alw); %<- Set properties
ylabel('$\epsilon_3$','interpreter','latex', 'LineWidth', alw,'FontSize', 14)
set(gca, 'Position', [0.73 0.3 0.225 0.6]); %<- Set properties
title('KKT Residuals','interpreter','latex','FontSize', fsz);
hold on 
    
lg = legend( 'r=1',...
    'r=5', ...
    'r=10','interpreter','latex','NumColumns',4,'Position',[0.115,0.03,0.8,0.1],'FontSize', fsz);

legend('boxoff');

%print("Numerical-result-SDP-different-para",'-depsc','-tiff');