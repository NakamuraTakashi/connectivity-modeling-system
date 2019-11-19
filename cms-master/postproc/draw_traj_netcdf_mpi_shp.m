close all
%% parameters
YEAR = '1993'

% nest_filename='../expt/expt_cots_2011/nests/nest_1_20110401000000.nc';
nest_filename=['../expt/expt_cots_' YEAR '/nests/nest_1_' YEAR '0401000000.nc'];

np = 10; % Number of CPU

% output_dir = '../expt/expt_cots_2011/output';
output_dir = ['../expt/expt_cots_' YEAR '/output'];

polyfile = 'D:/Documents/GIS_data/Yasuda_COTS/habitat.shp';

%% read in the nest data

% Get the values of Longitude from the nestfile
lonAxis = ncread(nest_filename,'Longitude');
% Get the values of Latitude from the nestfile
latAxis = ncread(nest_filename,'Latitude');
% Get the values of U-velocity from the nestfile
uvel = ncread(nest_filename,'zu',[1 1 1 1], [Inf Inf 1 1]);

%% draw the land

mask=squeeze(uvel).';
mask(mask<2^100)= 1;
mask(isnan(mask))=0;

xsize=800; ysize=530;
xmin=115;xmax=155;
ymin=15;ymax=40;

%color of the land in rgb color model
mymap = [1 1 1; 0.75 0.8 1];

f1=figure;
f1.Color=[1 1 1]; f1.Position=[0 0 xsize ysize];
f1.GraphicsSmoothing='off';
axes1 = axes('Parent',f1,...
    'FontSize',9,...
    'FontName','Arial',...
    'Box','on');
xlim(axes1,[xmin xmax]);
ylim(axes1,[ymin ymax]);

%draw the land and water
% surface(lonAxis,latAxis,mask','linestyle','none')
h_contour=contour(lonAxis, latAxis, mask,...
    'LineColor',[0.48 0.06 0.92],...
    'LevelList',[-1 1],...
    'Parent',axes1,...
    'ShowText','off');

xlabel('Longitude','FontName','Arial');
ylabel('Latitude','FontName','Arial');
hold on;

%% read the polygon shapefile
S = shaperead(polyfile);
n_poly = size(S,1);
for i=1:n_poly
    S(i).Lat = rmmissing(S(i).Y);
    S(i).Lon = rmmissing(S(i).X);
end
% Sort by id field
S2=S;
for i=1:n_poly-1
    for j=i+1:n_poly
        if S(i).id >S(j).id
            S2(i)=S(i);
            S(i)=S(j);
            S(j)=S2(i);
        end
    end
end

%different color for each trajectory
colors=jet(n_poly);

for j = 1:np
    if np >= 100      
        str_file_num = num2str(j,'%03d');
    elseif np >= 10
        str_file_num = num2str(j,'%02d');
    else
%         str_file_num = num2str(j,'%01d');
        str_file_num = num2str(j,'%02d');
    end
    
    traj_filename=[output_dir,'/traj_file_',str_file_num,'.nc'];

    %% read in the trajectory data

    % open trajectory file
    % Get the values of time
    time = ncread(traj_filename,'time');
    % Get the values of Longitude
    lon = ncread(traj_filename,'lon');
    % Get the values of Latitude
    lat = ncread(traj_filename,'lat');
    % Get the values of depths
    depth = ncread(traj_filename,'depth');
    % Get the values of status
    status = ncread(traj_filename,'exitcode');
    % Get the values of release date (in julian)
    release = ncread(traj_filename,'releasedate');
    %close nestfile
    polygon = ncread(traj_filename,'releasepolygon');

    lat(lat>999) = NaN;
    lon(lon>999) = NaN;


    %number of trajectories
    num_traj = size(lat,2);



    %% plot all the trajectories


    for i=1:num_traj
    %     plot(lon(:,i), lat(:,i),'color','green');
        if status(i) == -4
%            plot(lon(:,i), lat(:,i),'color',colors(polygon(i),:));
            plot(lon(:,i), lat(:,i),'color',colors(polygon(i),:));
%             plot(lon(1,i), lat(1,i),'Marker','o','color', colors(polygon(i),:));
%             plot(lon(size(lon,1),i), lat(size(lat,1),i),'Marker','*','color', colors(polygon(i),:));
        end
    end
end

%% Draw the polygon

for i=1:n_poly
    pgon = polyshape(S(i).Lon, S(i).Lat);
    plot(pgon,'FaceColor',colors(S(i).id,:),'FaceAlpha',0.5);
    text(mean(S(i).BoundingBox(:,1))-0.7, mean(S(i).BoundingBox(:,2)+0.2), num2str(S(i).id),'FontSize',11);
end


hold off;
drawnow
hgexport(figure(1), ['output/traj_' YEAR '.png'], hgexport('factorystyle'),'Format','png');

