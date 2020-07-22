clc;clear;
%% Environment Conditions.
global enviParams simTime T0 hub_height alpha
load('Norway5EC3D.mat')
load('Norway5params');

%%
T0          = 50;%(10:10:50); % return period
p3          = [0.5,];%[0.5,1e-1,1e-2,1e-3,1e-4,1e-5];
simTime     = 1; % simulation duration in hour
num_samples = 3600;
beta        = abs(norminv(simTime./(T0*365.25*24)));
theta       = linspace(1,360,num_samples)';
EC_U        = cell(size(beta));
EC_X        = cell(size(beta));
% phyVars = struct();
%% 10-meter wind speed U10 ~ truncated Weibull distribution

% turbine parameters
load('NREL5MW_OC4Semi')
% Uhub_in  = Uhub_cut(1);                                 % cut in wind speed
% Uhub_out = Uhub_cut(2);                                 % cut out wind speed
% hub_height = hub_height;                                % 133.5;% hub height
% Power-law wind shear profile
alpha = 0.1;                                            % exponential coeff
Uw_in   = Uhub_cutin /((hub_height / 10)^alpha);
Uw_out  = Uhub_cutout/((hub_height / 10)^alpha);



%% EC samples
U_hub = [3:11,11.3,12:25];
[EC_samples_U, EC_samples_X] = get_EC_samples(U_hub);

% %% Wave height, Hs|Uw ~ Weibull distribution
% % distribution parameters
% % alpha_u  = 2.029;                        % shape parameter
% % beta_u   = 9.409;                        % scale parameter
% U10_dist = makedist('Weibull','a', enviParams.Uw.beta, 'b', enviParams.Uw.alpha);
% % Uw_dist = truncate(Uw_dist, Uw_in,Uw_out);
% % q = linspace(0,1,10);
% % Uwq = icdf(Uw_dist,q);
% 
% % phyVars(end+1) = struct('name', 'Hs', 'unit','m','dist','None','distName','None', 'distParas','None');
% % distribution parameters
% % enviParams.Hs.a = [2.136,0.013,1.709];
% % enviParams.Hs.b = [1.816,0.024,1.787];
% for ii = 1 : length(p3)
%     EC_U{ii}        = sqrt(beta^2 - norminv(1-p3(ii))^2) * [cosd(theta), sind(theta)];
%     EC_X{ii}        = zeros (num_samples, 3); % initialize with (Uw, Hs, Tp);
%     EC_X{ii}(:,1)   = icdf(U10_dist,normcdf(EC_U{ii}(:,1))); % Rosenblatt transform for Uw
%     for i = 1 : num_samples
%         U10         = EC_X{ii}(i,1); % wind speed at 10 m height
%         HsUw_dist   = cal_dist_Hs(U10);
%         Hs          = icdf(HsUw_dist,normcdf(EC_U{ii}(i,2)));
%         U_hub       = U10 * ((hub_height / 10)^alpha);
%         Tp_HsUw_dist= cal_dist_Tp(U10,Hs);
%         Tp_median   = icdf(Tp_HsUw_dist,0.5);
%         % change to hub height wind speed
%         EC_X{ii}(i,1) = U_hub;
%         EC_X{ii}(i,2) = Hs;
%         EC_X{ii}(i,3) = Tp_median;
%     end
% end



% %%
% close all
% % figure
% % for i = 1 : size(U,2)
% %     plot(U{i}(:,1),U{i}(:,2))
% %     hold on
% % end
% % legend('p_3=0.5','p_3=1e-1','p_3=1e-2','p_3=1e-3','p_3=1e-4','p_3=1e-5')
% % hold off
% textpos = [5,7;25,11.3;5,5.8;25,9.78;5,4.2;21,4.9]; 
% textcon = {'$p_3=0.5$','$p_3=1e-1$','$p_3=1e-2$','$p_3=1e-3$','$p_3=1e-4$','$p_3=1e-5$'};
% figure
% for i = 1 : size(EC_U,2)
%     plot(EC_X{i}(:,1),EC_X{i}(:,2),'k','LineWidth',1)
%     text(textpos(i,1),textpos(i,2),textcon{i},'Interpreter','latex','FontSize',15);
%     hold on
% end
% 
% insdieInx = [9,11:17];
% 
% plot(phyRvs(insdieInx,1),phyRvs(insdieInx,2)-0.5,'o','MarkerSize',8,'MarkerEdgeColor','b','MarkerFaceColor','w','LineWidth',2)
% plot(phyRvs(insdieInx,1),phyRvs(insdieInx,2)-1,'o','MarkerSize',8,'MarkerEdgeColor','b','MarkerFaceColor','w','LineWidth',2)
% plot(phyRvs(insdieInx,1),phyRvs(insdieInx,2)-2,'o','MarkerSize',8,'MarkerEdgeColor','b','MarkerFaceColor','w','LineWidth',2)
% plot(phyRvs(insdieInx,1),phyRvs(insdieInx,2)-3,'o','MarkerSize',8,'MarkerEdgeColor','b','MarkerFaceColor','w','LineWidth',2)
% plot(phyRvs(insdieInx,1),phyRvs(insdieInx,2)-4,'o','MarkerSize',8,'MarkerEdgeColor','b','MarkerFaceColor','w','LineWidth',2)
% plot(phyRvs(insdieInx,1),phyRvs(insdieInx,2)-5,'o','MarkerSize',8,'MarkerEdgeColor','b','MarkerFaceColor','w','LineWidth',2)
% insdieInx = [9,11:17];
% plot(phyRvs(:,1),phyRvs(:,2),'o','MarkerSize',8,'MarkerEdgeColor','r','MarkerFaceColor','r')
% % plot(phyRvs(insdieInx,1),phyRvs(insdieInx,2)-1,'o','MarkerSize',7,'MarkerEdgeColor','b','MarkerFaceColor','g')
% % plot(phyRvs(insdieInx,1),phyRvs(insdieInx,2)-2,'o','MarkerSize',7,'MarkerEdgeColor','b','MarkerFaceColor','g')
% 
% xlim([2,28])
% xlabel('Hub height mean wind speed, $V(m/s)$','Interpreter','latex')
% ylabel('Significant wave height , $H_s(m)$','Interpreter','latex')
% title('Environment Contour, $T = 50$ years ','Interpreter','latex')
% % lines = findobj(gcf,'type','line');
% % lgd=legend('$p_3=0.5$','$p_3=1e-1$','$p_3=1e-2$','$p_3=1e-3$','$p_3=1e-4$','$p_3=1e-5$','Location','best');
% % set(lgd,'Interpreter','latex')
% % legend(lines([1,3,4]),'3D Inverse FORM','2D Inverse FORM',' Response surfaces ','Interpreter','latex')
% % set(gcf, 'Position',[50,50,1250,950])
% % set(gca,'FontSize',20,'FontName','Times New Roman')
% grid on
% hold off
%%
function [EC_U, EC_X] = get_EC_samples(U_hub)
    global enviParams simTime T0 hub_height alpha
    U10  = U_hub * ((hub_height / 10)^(-alpha));
    beta = abs(norminv(simTime ./(T0*365.25*24)));
    EC_X = zeros(length(U10), 3);
    EC_X(:,1)= U_hub;
    U10_dist = makedist('Weibull','a', enviParams.Uw.beta, 'b', enviParams.Uw.alpha);
    U1 = norminv(cdf(U10_dist,U10));
    U2 = sqrt(beta^2 - U1.^2); %% U3 assume to be 0, median value
    U3 = 0 * U2;
    EC_U = [U1;U2;U3]';
    for ii = 1 : length(U10)
        dist_Hs = cal_dist_Hs(U10(ii));
        Hs      = icdf(dist_Hs,normcdf(U2(ii)));
        dist_Tp = cal_dist_Tp(U10(ii),Hs);
        Tp      = icdf(dist_Tp,normcdf(U3(ii)));   
        EC_X(ii,2:3) = [Hs, Tp];
    end

end

function HsUw_dist = cal_dist_Hs(U10)
    global enviParams
    alpha_hu    = enviParams.Hs.a(1) + enviParams.Hs.a(2) * U10 ^(enviParams.Hs.a(3));
    beta_hu     = enviParams.Hs.b(1) + enviParams.Hs.b(2) * U10 ^(enviParams.Hs.b(3));
    HsUw_dist   = makedist('Weibull','a', beta_hu, 'b', alpha_hu);
end

function dist_Tp = cal_dist_Tp(U10, Hs)
    global enviParams
    theta = enviParams.Tp.theta;
    gamma = enviParams.Tp.gamma;
    e = enviParams.Tp.e;
    f = enviParams.Tp.f;
    k = enviParams.Tp.k;
    Tp_niu    = k(1) + k(2) * exp(Hs * k(3));  % eqn 21
    U10_bar   = f(1) + f(2) * Hs^f(3);
    Tp_bar    = e(1) + e(2) * Hs^e(3);
    
    Tp_mu = Tp_bar * (1 + theta * ( (U10 - U10_bar)/U10_bar )^gamma);
    Tp_sigma = Tp_niu * Tp_mu;
    
    log_Tp_mu = log(Tp_mu / sqrt(1 + Tp_niu^2));
    log_Tp_sigma = log(Tp_niu^2 + 1);
    dist_Tp = makedist('Lognormal', 'mu', log_Tp_mu, 'sigma', log_Tp_sigma);

end
