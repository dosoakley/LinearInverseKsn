%This is a heavily modified form of the INCA code from Smith and Fox (2024).
%It is faster, has some adjustments to scale factors, and also contains several new options.
%It may not give exactly the same results as INCA, but it should be close.

clear all
close all

%%% Options. Change these as needed %%%
tt_dir = 'PATH_TO_TOPOTOOLBOX'; %Path for the directory in which topotoolbox is located.
DEM_path = 'PATH_TO_DEM/clearwater_2.tif'; %Path for the DEM.
mn_vec = 0.3:0.02:0.8; %Vector of m/n values to try, in the form min:step:max
alpha = 1; %The regularization parameter. Note: 1 here is roughly equivalent to 1e4 in INCA, because I also removed the scaling of chi by 1e4.
minarea = 1e6; %Minimum area to define a stream (m^2).
algorithm = 'interior-point'; %Algorithm for least-squares minimization.
display = 'none'; %Display option for lsqlin.
project_to_UTM = false; %Whether to reproject the DEM to UTM. Note: This uses reproject2utm, which expects that the input is in a WGS84 geographic coordinate system, which the Clearwater DEM isn't.
use_full_resnorm = true; %If true, use the full residual norm to choose the best m/n . If false, use only the residual norm for the stream channel elevations, not including the regularization.
DEM_cellsize = 0; %In meters. Set to 0 to determine it from the DEM.
soln_cellsize = 2000; %In meters, size of the grid on which to solve for ksn.
boundary_file = 'PATH_TO_Boundary/outlet.kml'; %Optional, set to '' to not use. Path to a shape or kml file defining a boundary. Only streams upstream of the boundary will be considered.
klargest = 1; %Number of largest connected stream networks to use. Set to 0 to use all.

%%% Run the analsysis. %%%
%The results will be stored in a structure called results.
tic

%Add the path to topotoolbox
addpath([tt_dir,'/topotoolbox-master'])
addpath([tt_dir,'/topotoolbox-master/utilities'])

%Load the DEM, and extract the stream network.
[S,z,Area,DEM] = GetStreamsFromDEM(DEM_path,minarea,'boundary_file',boundary_file,'cellsize',DEM_cellsize, ...
    'klargest',klargest,'reproject',project_to_UTM);

%Loop through Ksn values, do the linear inverse Ksn for each one, and find
%the best one.
results = ConstrainConcavity(S,z,Area,mn_vec,soln_cellsize,'alpha',alpha,...
    'A0',minarea,'algorithm',algorithm,'display',display,'use_full_resnorm',use_full_resnorm);

%%% Make some plot. %%%

%Plot the misfit.
figure(1)
ax = gca;
ax.FontSize = 16;
hold on
box on
plot(mn_vec, results.misfit_matrix, 'LineWidth',2)
xlim([mn_vec(1),mn_vec(end)])
xlabel("\theta", 'FontSize',20)
ylabel("RMS Misfit [m]")
% savefig(gcf,'Misfit.fig')

%Plot the roughness.
%It makes more sense to plot the scaled roughness. The unscaled roughness
%is in terms of Ksn, and therefore it increases as Ksn increases.
figure(2)
ax = gca;
ax.FontSize = 16;
hold on
box on
plot(mn_vec, results.roughness_matrix_scaled, 'LineWidth',2)
xlim([mn_vec(1),mn_vec(end)])
xlabel("Roughness", 'FontSize',20)
ylabel("RMS Roughness")
% savefig(gcf,'Roughness.fig')

%Plot ksn for the best theta value.
figure(3)
ksn_grid=reshape(results.ksn_best(1:results.nx_grid*results.ny_grid),results.nx_grid,results.ny_grid);
imageschs(DEM,[],'colormap','gray') %It is possible to give this a second argument to plot something directly on top.
hold on
x_ksn = results.x_grid+(results.x_grid(2,1)-results.x_grid(1,1))/2; %x_grid and y_grid give the coordinates of the bottom left column. Get the center coordinates instead.
y_ksn = results.y_grid+(results.y_grid(1,2)-results.y_grid(1,1))/2;
s = pcolor(x_ksn,y_ksn,ksn_grid);
colormap('parula') %After running imageschs, the colormap and its limits get set by that, so we need to reset them for Ksn.
clim([min(ksn_grid(:)),max(ksn_grid(:))])
s.FaceAlpha = 0.75;
plot(S,'k');  %Plot the stream network.
axis equal
if strcmp(DEM.georef.SpatialRef.CoordinateSystemType,'geographic')
    xlabel('Longitude')
    ylabel('Latitude')
else
    xlabel('Easting (m)')
    ylabel('Northing (m)')
end
hold off
colorbar

t = toc;
disp(' ')
disp(['Best m/n = ',num2str(results.theta_best)])
disp(['Total Run Time = ',num2str(toc),' seconds.'])
