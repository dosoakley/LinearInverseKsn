function [S,z,Area,varargout] = GetStreamsFromDEM(DEM_path,minarea,options)

%Warning: If boundary_file is specified it must have the same coordinate
%system as the DEM. This will not be checked.

%The optional output in varargout is the DEM itself. Use [S,z,Area,DEM] as
%the outputs if you need the DEM or just [S,z,Area] if you only need the
%stream network and don't want to have to keep storing the DEM in memory.

arguments
    DEM_path {mustBeFile} %The path to the DEM file.
    minarea (1,1) {mustBePositive} %The minimum drainage area to use to define a stream (in square meters).
    options.boundary_file char = '' %The path to a shapefile or kml file that defines a boundary. Only streams above the boundary will be used.
    options.reproject (1,1) logical = false %Whether to reproject the map to UTM. (Note: This may require the DEM to be in a WGS84 coordinate system.)
    options.cellsize (1,1) {mustBeNonnegative} = 0 %The cellsize in meters for the output DEM. Set to 0 to determine automatically from the DEM.
    options.klargest (1,1) {mustBeInteger} = 0 %Extract only the klargest largest connected components in the stream network. Set to 0 to ignore and use all.
end
cellsize_from_DEM = options.cellsize == 0; %Whether to determine the cellsize automatically from the DEM.

%Load the DEM and process it to get streams.
DEM = GRIDobj(DEM_path);
%Get the boundary if one is specified.
if ~isempty(options.boundary_file)
    %Load the boundary file before reprojecting DEM so we can use
    %line2grid with DEM and boundary both in geographic coordinates.
    if strcmp(options.boundary_file(end-2:end),'shp')
        boundary = shaperead(options.boundary_file);
    elseif strcmp(options.boundary_file(end-2:end),'kml')
        boundary = kml_shapefile(options.boundary_file);
    else
        error('boundary_file type is unrecognized.')
    end
    boundary_line = line2GRIDobj(DEM, boundary);
    boundary_line = dilate(boundary_line,ones(2)); %Make the boundary line twice as wide - prevents streams that flow diagonally being missed
end
%Reproject the DEM to UTM, if desired.
%degree to metre calculation
radiusEarth = 6378e3; % Mean radius of the Earth (m)
mPerDegree = pi*radiusEarth/180; % Meters per arc degree for the Earth's surface
input_DEM_geographic = strcmp(DEM.georef.SpatialRef.CoordinateSystemType,'geographic');
if options.reproject
    if ~input_DEM_geographic
        disp('Warning: The DEM will not be reprojected. Only geographic coordinate systems can be reprojected.')
    else
        if ~strcmp(DEM.georef.SpatialRef.GeographicCRS.Name,'WGS 84')
            disp('Warning: The DEM is not in WGS 84. I am not sure if reproject2utm will work correctly.')
        end
        if cellsize_from_DEM
            DEM_cellsize = DEM.cellsize * mPerDegree;
            disp(['  DEM cellsize = ',num2str(DEM_cellsize),' meters.'])
        else
            DEM_cellsize = options.cellsize;
        end
        [DEM,zone] = reproject2utm(DEM,DEM_cellsize); %TopoToolbox prefers a projected coordinate system, and I think snap2stream may need it.
        disp(['DEM has been reprojected to zone',zone])
        if ~isempty(options.boundary_file)
            boundary_line = reproject2utm(boundary_line,DEM); %Would it be better to the dilate before or after the reprojection?
        end
    end
end
%Set the cellsize for the DEM.
if ~options.reproject || ~input_DEM_geographic %If we reprojected, cell size was already set while doing that.
    DEM_cellsize = DEM.cellsize;
    if input_DEM_geographic
        DEM_cellsize = DEM_cellsize * mPerDegree; %Convert to meters.
    end
    if cellsize_from_DEM
        disp(['  DEM cellsize = ',num2str(DEM_cellsize),' meters.'])
        DEM.cellsize = DEM_cellsize; %Set it in meters rather than degrees.
    else
        if DEM_cellsize ~= options.cellsize
            %This is mostly in here because it is the way that the original INCA code worked.
            disp('Warning: The specified cellsize does not match the cellsize of the DEM.')
            disp(['  DEM cellsize = ',num2str(DEM_cellsize),' meters.'])
            disp(['  Specified cellsize = ',num2str(options.cellsize),' meters.'])
            disp('  The cellsize for the DEM GRIDobj will be set to the specified value anyway, but this is likely to be wrong.')
        end
        DEM.cellsize = options.cellsize; %Set it in meters rather than degrees.
    end
end
%Set areas below sea level to NaN.
DEM.Z(DEM.Z <= 0) = NaN;
%Get the drainage area.
FD = FLOWobj(DEM,'preprocess','carve');
DEM = imposemin(FD,DEM,0.0001);
AS = flowacc(FD);
%extract stream network
S = STREAMobj(FD, 'minarea', minarea, 'unit', 'map');
if ~isempty(options.boundary_file)
    S = modify(S, 'upstreamto', boundary_line);
end
if options.klargest > 0
    S = klargestconncomps(S, options.klargest);
end
z = getnal(S,DEM);
Area = getnal(S,AS);
varargout = {DEM};
end