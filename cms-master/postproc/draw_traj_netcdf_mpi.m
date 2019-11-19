%% parameters
nest_filename='../expt/expt_cots_2018/nests/nest_1_20180401000000.nc';
np = 10; % Number of CPU
output_dir = '../expt/expt_cots_2018/output';
polyfile = 'D:/Documents/GIS_data/Yasuda_COTS/habitat.shp';

n_polygon = 29;

%% read in the nest data

% Get the values of Longitude from the nestfile
lonAxis = ncread(nest_filename,'Longitude');
% Get the values of Latitude from the nestfile
latAxis = ncread(nest_filename,'Latitude');
% Get the values of U-velocity from the nestfile
uvel = ncread(nest_filename,'zu',[1 1 1 1], [Inf Inf 1 1]);

%% draw the land

%layer of depth. layer = 1 is the surface.
%     layer=1;
%divide velocities in land(value=1) and water(value=0)
%velocity of land is 2^100
mask=squeeze(uvel);
mask(mask<2^100)=0;
mask(mask>=2^100)= 1;
%draw the land and water
figure;
% contourf(lonAxis,latAxis,mask',[0.5 0.5],'linestyle','none');
surface(lonAxis,latAxis,mask','linestyle','none');
%color of the land in rgb color model
mymap = [1 1 1; 0.75 0.8 1];
colormap(mymap);

hold on;

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

    %different color for each trajectory
    colors=jet(n_polygon);

    %% plot all the trajectories


    for i=1:num_traj
    %     plot(lon(:,i), lat(:,i),'color','green');
        if status(i) == -4
            plot(lon(:,i), lat(:,i),'color',colors(polygon(i),:));
%             plot(lon(1,i), lat(1,i),'Marker','o','color', colors(polygon(i),:));
%             plot(lon(size(lon,1),i), lat(size(lat,1),i),'Marker','*','color', colors(polygon(i),:));
        end
    end
end
hold off;

%print title
title(['Plotted ',num2str(num_traj),' trajectories']);
axis equal;
