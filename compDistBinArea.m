function phyArea = compDistBinArea(enviParams, xgrids, returnPeriod, plotParams, simDuration, turbineModel, pwrLawAlpha )
% Plot 2D environmental contour and sea states 
% Arguments:
%   enviParams:     struct contains environment distribution paramesters
%   xgrids:         sea states executed
%   returnPeriod:   list of targeted return periods
%   simDuration:    length of simulation time for one sample
%   innerFractiles: targeted fractiles inside the contour 
%   turbineModel:   SNL or NREL, used to define hub height
%   pwrLawAlpha:    exponent index of power law wind profile
% Return:
%   enviContourPhyVars: array of shape (No. of betas, 
%                                       No. of innerFractiles, 
%                                       num_precision, 
%                                       size of physical variables (Uw,Hs));
% 
%   Example: compDistBinArea('Norway5params','Norway5EC2D');
fprintf('Calculating and plotting environment contour...\n')
load(enviParams, 'enviParams');

if nargin < 2 || isempty(xgrids)
    phyRvs = [];
else
    load(xgrids,'phyRvs');
    if contains(xgrids,'3D')
        markerProp = {'b.'};
    elseif contains(xgrids, '2D')
        markerProp = {'ro', 'MarkerSize',5, 'MarkerEdgeColor','r', 'MarkerFaceColor','r'};
    end
end

if nargin < 3 || isempty(returnPeriod)
    returnPeriod = 50;  % return period
end

if nargin < 4 || isempty(plotParams)
    plotParams.isPlot = 1; % plot figures in default
    plotParams.xlabel = 'Hub height mean wind speed, $V(m/s)$';
    plotParams.ylabel = 'Significant wave height , $H_s(m)$';
    plotParams.title  = 'Environment Contour, $T = 50$ years';
    plotParams.xlim   = ([2,28]);

end

if nargin < 5
    simDuration = 1;
end

if nargin < 6
    turbineModel = 'SNL';
end
if contains(turbineModel,'SNL')

    hubHeight = 133.5;
elseif contains(turbineModel,'NREL')
    hubHeight = 90;
%     hubHeight = 133.5;
else
    fprintf('Error: wind turbine model hub height is not defined')
end

if nargin < 7
    pwrLawAlpha = 0.1;
end


fprintf('\t%-25s:  %-30s\n', 'Wind turbine model',          turbineModel);
fprintf('\t%-25s:  %-30s\n', 'Return Period',               returnPeriod);
fprintf('\t%-25s:  %-30s\n', 'Simulation duration',         num2str(simDuration));
fprintf('\t%-25s:  %-30s\n', 'Turbine Model',               turbineModel);
fprintf('\t%-25s:  %-30s\n', 'Turbine hub height',          num2str(hubHeight));
fprintf('\t%-25s:  %-30s\n', 'Power law expoent',           num2str(pwrLawAlpha));






% Load environment (Uw, Hs, Tp) conditional distributions 
% load(enviParams,'enviParams');
Uw_dist = makedist('Weibull','a', enviParams.Uw.beta, 'b', enviParams.Uw.alpha);
u_dist  = makedist('Normal','mu',0,'sigma',1);
% num_precision   = 3600;
% theta = linspace(1,360,num_precision)';
% Calculate reliability indices, one for each return period: beta
beta = abs(norminv(simDuration./(returnPeriod*365.25*24)));

dUw = 1;
dHs = 0.5;
phyArea = 0*phyRvs;
for iss = 1 : size(phyRvs,1)
    
    iUw = phyRvs(iss,1);
    iHs = phyRvs(iss,2);
    meshX = [max(iUw - dUw/2,3), min(iUw + dUw/2,25)];
    meshY = [max(iHs - dHs/2,0), iHs + dHs/2];

    alpha_hu    = enviParams.Hs.a(1) + enviParams.Hs.a(2) * iUw^(enviParams.Hs.a(3));
    beta_hu     = enviParams.Hs.b(1) + enviParams.Hs.b(2) * iUw^(enviParams.Hs.b(3));
    HsUw_dist   = makedist('Weibull','a', beta_hu, 'b', alpha_hu);
    u1 = icdf(u_dist,cdf(Uw_dist,iUw));
    u2 = icdf(u_dist,cdf(HsUw_dist,iHs));
    if u1^2 + u2^2 > beta^2
        u3 = 0;
        iss
    else
        u3 = sqrt(beta^2 - u2^2 - u1^2);
    end
    p3 = 1-cdf(u_dist,u3);
    phyArea(iss,1) = p3;
    nsubsize = 1000;
    area = 0;

    for ii = linspace(meshX(1),meshX(2),nsubsize)
       for jj = linspace(meshY(1),meshY(2),nsubsize)
           area = area + pdf(Uw_dist, ii)*pdf(HsUw_dist, jj);
       end
    end
    area = area * ((meshX(2)-meshX(1))/(nsubsize-1)) * ((meshY(2)-meshY(1))/(nsubsize-1));
    phyArea(iss,2) = area;
    
end
phyArea(:,3) = phyArea(:,2) * (50*365.25*24);
phyArea(:,4) = phyArea(:,1) .* phyArea(:,3);
%     function A = compArea(X,Y)
% enviContourPhyVars = zeros( size(beta,2),...
%                             num_precision,...
%                             size(phyRvs,2));

% for ibeta = 1:size(beta,2)
%     b = beta(1,ibeta);
%     fprintf(['\tCalculating contour for return period: ' num2str(returnPeriod(ibeta)),' years\n']);
%     % U vector in standard normal space
%     U = sqrt(b^2 - norminv(1 - q)^2) * [cosd(theta), sind(theta)];
% %         wind speed in physical space
%     enviContourPhyVars(ibeta,iFractile,:,1) =  icdf(Uw_dist,normcdf(U(:,1)));
%     for itheta = 1 : num_precision
%         Uw          = enviContourPhyVars(ibeta,iFractile,itheta,1);
%         alpha_hu    = enviParams.Hs.a(1) + enviParams.Hs.a(2) * Uw^(enviParams.Hs.a(3));
%         beta_hu     = enviParams.Hs.b(1) + enviParams.Hs.b(2) * Uw^(enviParams.Hs.b(3));
%         HsUw_dist   = makedist('Weibull','a', beta_hu, 'b', alpha_hu);
% 
%         enviContourPhyVars(ibeta, iFractile, itheta, 2) = icdf(HsUw_dist,normcdf(U(itheta,2)));
%         % change Uw to hub height wind speed, Power law profile
%         enviContourPhyVars(ibeta, iFractile, itheta, 1) = enviContourPhyVars(ibeta, iFractile, itheta, 1) * ((hubHeight / 10)^pwrLawAlpha);
%     end
% end


%%
% if plotParams.isPlot
%     fprintf('\tMaking plots...\n')
%     for ibeta = 1:size(beta,2)
% %         close all
% %         figure
% 
%         hold on
%         for iFractile = 1:size(innerFractiles,2)
%             data2plot = squeeze(enviContourPhyVars(ibeta, iFractile, :, :));
%             plot(data2plot(:, 1),data2plot(:, 2),'k','LineWidth',1);
%         end
%         if ~isempty(phyRvs)
% %             plot(phyRvs(:,1),phyRvs(:,2),'b.')
%             plot(phyRvs(:,1),phyRvs(:,2),markerProp{:});
% %             xlim([floor(min(phyRvs(:,1)))-1, ceil(max(phyRvs(:,1)))+1])
%         end
% 
%         xlabel(plotParams.xlabel);
%         ylabel(plotParams.ylabel);
%         title(plotParams.title);
%         xlim(plotParams.xlim);
%         grid on
%         hold off
%         saveas(gcf,['Figures/EnviContour', num2str(returnPeriod(ibeta)),'.png'])
%         saveas(gcf,['Figures/EnviContour', num2str(returnPeriod(ibeta)),'.fig'])
%         saveas(gcf,['Figures/EnviContour', num2str(returnPeriod(ibeta)),'.eps'], 'epsc')
%     end
% 
%     %%
% 
%     
% end


end