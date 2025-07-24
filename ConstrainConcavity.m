function [results] = ConstrainConcavity(S,h,Area,mn_values,grid_cellsize,options)
%ConstrainConcavity Map Ksn and find best theta.
%   Following Smith et al. (2024), use a linear inverse method to produce a
%   map of Ksn values, and test different concavity (theta) values to see
%   which gives the best fit.

arguments
    S {mustBeA(S,'STREAMobj')} %A topotoolbox stream object.
    h {mustBeNAL(S,h)} %A node-ordered list of elevation corresponding to the nodes in S. Must be in either a geographic coordinate system or a projected coordinate system with units of meters.
    Area {mustBeNAL(S,Area)} %A node-ordered list of drainage areas corresponding to the nodes in S. (In pixels, as returned by floacc.)
    mn_values (:,1) {mustBeReal} %A 1D array of the values of m/n (theta) to test.
    grid_cellsize (1,1) {mustBeReal} %The cellsize of the Ksn grid.
    options.alpha (1,1) {mustBeReal} = 1 %The regularization parameter.
    options.A0 (1,1) {mustBeReal} = 1e6; %The reference area for the chi transformation.
    options.algorithm {mustBeMember(options.algorithm,['interior-point',...
        'trust-region-reflective','active-set'])} = 'interior-point'; %Algorithm for least-squares minimization.
    options.display {mustBeMember(options.display,['off','none','final','iter','iter-detailed','final-detailed'])} = 'none'; %Display option for lsqlin.
    options.use_full_resnorm (1,1) logical = true; %If true, use the full residual norm to choose the best m/n . If true, use only the residual norm for the stream channel elevations, not the regularization.
end

%Get the options.
[alpha,A0,algorithm,display] = deal(options.alpha,options.A0,options.algorithm,options.display);


%If S uses a geographic coordinate system, convert degrees to meters.
%Otherwise, assume that the units of S are in meters.
if strcmp(S.georef.SpatialRef.CoordinateSystemType,'geographic')
    %degree to metre calculation
    %... Mean radius of the Earth (m)
    radiusEarth = 6378e3;
    %... Meters per arc degree for the Earth's surface
    conversion = pi*radiusEarth/180;
    disp(['Degrees will be converted to meters by multiplying by ',num2str(conversion)])
else
    conversion = 1; %No need to convert.
end


%define elevations
z = zerobaselevel(S, h); %Adjust elevations so that all base levels are 0.

%produce a linear list of IDs
id = (1:length(S.ix))';
number_of_nodes=size(id,1);
disp(['number of nodes = ',num2str(number_of_nodes)])

misfit_matrix = zeros(1, numel(mn_values));
misfit_matrix_scaled = zeros(1, numel(mn_values));
roughness_matrix = zeros(1, numel(mn_values));
roughness_matrix_scaled = zeros(1, numel(mn_values));

%Build the grid on which to solve for u* (ksn).
grid_size = grid_cellsize/conversion;
[x1,x2] = deal(min(S.x),max(S.x));
[y1,y2] = deal(min(S.y),max(S.y));
xL=x2-x1;
yL=y2-y1;
nx_grid=1 + floor((xL)/grid_size);
ny_grid=1 + floor((yL)/grid_size);
nn_grid=nx_grid*ny_grid;
[x_grid,y_grid] = ndgrid(0:grid_size:xL,0:grid_size:yL); %x_grid and y_grid give the bottom left coordinates of each cell.

%Build the W matrix since we only need to do this once.
W = MakeW(nx_grid,ny_grid,alpha);

%%%%% inverse scheme begins here %%%%%%%%%%%%%%%%%%%
for n_mn = 1:numel(mn_values)
    clear C
    clear Model_damp

    mn = mn_values(n_mn);
    C = chitransform(S, Area, 'mn', mn,'a0',A0);
    C = C*conversion; %Convert degrees to meters, if necessary.
    chi_scale_factor = 1/max(C);
    C = C*chi_scale_factor;
    %perform palaeotopography calculation
    chi = C;
    delta_chi = zeros(size(chi));
    delta_chi(S.ix) = chi(S.ix) - chi(S.ixc); %Difference between chi at each node and its receiver. Base level nodes will stay at 0.

    %%% BUILD MODEL MATRIX AND DAMPING MATRIX!!!%%%%%%

    %First loop through to determine structure
    A = MakeA(S,delta_chi,grid_size,nn_grid,nx_grid,x1,y1);

    fprintf(1,'A matrix built\n');

    %Scale z and pad it with zeros for the regularization terms.
    z_scale_factor = 1/max(z);
    z_scaled = z(S.ix) * z_scale_factor; %By taking z(S.ix), this excludes the base-level nodes.
    elev_padded = [z_scaled;zeros(nn_grid,1)];

    Model_damp=[A;W];

    fprintf(1,'Solving system');
    lb = zeros([nn_grid,1]); %Lower bound. Ksn cannot go negative.
    lsqlin_options = optimoptions('lsqlin','Algorithm',algorithm,'Display',display);
    [ksn_scaled,resnorm]=lsqlin(Model_damp,elev_padded,[],[],[],[],lb,[],ones(nn_grid,1),lsqlin_options);

    %Calculate the un-scaled ksn.
    %The elevation is z = sum(delta_chi*ksn/a0^mn).
    %Because we used chi * chi_scale_factor and z * z_scale_factor, what we fit for is z_scale_factor * (ksn / (a0^mn * chi_scale_factor)).
    ksn_scale_factor = z_scale_factor / (chi_scale_factor * A0^mn); %The ksn we fit for is actually ksn * ksn_scale_factor.
    ksn = ksn_scaled / ksn_scale_factor;

    %Calculate the RMS misfit in z.
    misfit_scaled = z_scaled - A*ksn_scaled; %Misfit in terms of scaled elevation.
    if ~options.use_full_resnorm
        resnorm = sum(misfit_scaled.^2); %Sum of squares for only the (scaled) z misfit part.
    end
    misfit = misfit_scaled / z_scale_factor; %Misfit in terms of absolute elevation.
    RMS_misfit_scaled = sqrt(mean(misfit_scaled.^2));
    RMS_misfit = sqrt(mean(misfit.^2));
    misfit_matrix_scaled(n_mn) = RMS_misfit_scaled;
    misfit_matrix(n_mn) = RMS_misfit;
    
    %See if this is a new best fit.
    if n_mn == 1 || resnorm < resnorm_best
        ksn_best_scaled = ksn_scaled;
        ksn_best = ksn;
        theta_best = mn;
        resnorm_best = resnorm;
    end

    roughness_ksn_scaled= W*ksn_scaled;
    roughness_ksn = roughness_ksn_scaled / ksn_scale_factor / grid_cellsize.^2; %Convert the roughness to units of ksn/m^2.
    RMS_roughness_scaled = sqrt(mean(roughness_ksn_scaled.^2));
    RMS_roughness = sqrt(mean(roughness_ksn.^2));
    roughness_matrix_scaled(n_mn) = RMS_roughness_scaled;
    roughness_matrix(n_mn) = RMS_roughness;

end %this is end for mn

%Assemble the results into a structure to return.
results.ksn_best_scaled = ksn_best_scaled;
results.ksn_best = ksn_best;
results.theta_best = theta_best;
results.resnorm_best = resnorm_best;
results.misfit_matrix_scaled = misfit_matrix_scaled;
results.misfit_matrix = misfit_matrix;
results.roughness_matrix_scaled = roughness_matrix_scaled;
results.roughness_matrix = roughness_matrix;
results.nx_grid = nx_grid;
results.ny_grid = ny_grid;
results.x_grid = x1+x_grid; %Add x1 back in so it is shifted to the correct real-world position.
results.y_grid = y1+y_grid; %Add y1 back in so it is shifted to the correct real-world position.
end

%Define a validation function for ensuring that an input is a node
%attribute list of the stream object.
function mustBeNAL(S,a)
    if ~(isnal(S,a))
        eidType = 'mustBeNAL:notNAL';
        msgType = 'Input must be a node attribute list of the stream object.';
        error(eidType,msgType)
    end
end