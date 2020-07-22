% mu = [0 2];
% Sigma = [.25 .3; .3 1];
% x1 = -3:.2:3; x2 = -3:.2:3;
% [X1,X2] = meshgrid(x1,x2);
% F = mvnpdf([X1(:) X2(:)],mu,Sigma);
% F = reshape(F,length(x2),length(x1));
% surf(x1,x2,F);
% caxis([min(F(:))-.5*range(F(:)),max(F(:))]);
% axis([-3 3 -3 3 0 .4])
% xlabel('x1'); ylabel('x2'); zlabel('Probability Density');
close all
% cmap = jet(20);
% cmap = flipud(cmap(1:10,:));
% cmap(1,:) = [1,1,1];
% colormap(cmap);



load('Norway5params');
Uw_dist = makedist('Weibull','a', enviParams.Uw.beta, 'b', enviParams.Uw.alpha);

ngrid = 200;
Uw = linspace(0,45,ngrid*2);
Hs = linspace(0,16,ngrid);

[HsGrids,UwGrids] = meshgrid(Hs,Uw);
HsUwJointPdf = 0*UwGrids;
Uw_pdf = pdf(Uw_dist, Uw);
for igrid = 1 : length(Uw)
    iUw = Uw(igrid);
    alpha_hu    = enviParams.Hs.a(1) + enviParams.Hs.a(2) * iUw^(enviParams.Hs.a(3));
    beta_hu     = enviParams.Hs.b(1) + enviParams.Hs.b(2) * iUw^(enviParams.Hs.b(3));
    HsUw_dist   = makedist('Weibull','a', beta_hu, 'b', alpha_hu);
    HsUw_pdf    = pdf(HsUw_dist,Hs) * Uw_pdf(igrid);
    HsUwJointPdf(igrid,:) = HsUw_pdf;
end



% jointPdfSurf = surf(UwGrids*(90 / 10)^0.1, HsGrids,HsUwJointPdf);
% jointPdfSurf.EdgeColor = 'none';
% jointPdfSurf.LineStyle = 'none';

% Remove first row and column, pdf cannot start with 0 when apply log
UwGrids = UwGrids(2:end, 2:end);
HsGrids = HsGrids(2:end, 2:end);
HsUwJointPdf = HsUwJointPdf(2:end, 2:end);
HsUwJointPdf(HsUwJointPdf<1e-10) = 1e-10;
contourf(UwGrids*(90 / 10)^0.1, HsGrids,log10(HsUwJointPdf));


% contourf(UwGrids*(90 / 10)^0.1, HsGrids,HsUwJointPdf);
colormap(bluewhitered)
% caxis([min(F(:))-.5*range(F(:)),max(F(:))]);
% axis([-3 3 -3 3 0 .4])
xlabel('Uw'); ylabel('Hs'); zlabel('Probability Density');view(2);
hold on;