clear all
%% parameters
nest_filename='../expt/expt_cots_2018/nests/nest_1_20180401000000.nc';
ncpu = 10; % Number of CPU
output_dir = '../expt/expt_cots_2018/output';
polyfile = 'D:/Documents/GIS_data/Yasuda_COTS/habitat.shp';

%% read in the nest data

% Get the values of Longitude from the nestfile
lonAxis = ncread(nest_filename,'Longitude');
% Get the values of Latitude from the nestfile
latAxis = ncread(nest_filename,'Latitude');
% Get the values of U-velocity from the nestfile
uvel = ncread(nest_filename,'zu',[1 1 1 1], [Inf Inf 1 1]);
nlon=size(lonAxis);
nlat=size(latAxis);
for i=1:nlat
   lon_map(:,i)=lonAxis; 
end
for i=1:nlon
   lat_map(i,:)=latAxis; 
end

%% draw the land
close all

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
f1.Colormap=colors;
%% Draw scatter plot
h_scatter=scatter(0,0,3,0,'fill'); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Point size

%% Draw the polygon

for i=1:n_poly
    pgon = polyshape(S(i).Lon, S(i).Lat);
    plot(pgon,'FaceColor',colors(S(i).id,:),'FaceAlpha',0.5);
    text(mean(S(i).BoundingBox(:,1))-0.7, mean(S(i).BoundingBox(:,2)+0.2), num2str(S(i).id),'FontSize',11);
end

drawnow;

%% Read data
polygon2=[];
status2=[];
release2=[];
    
for j = 1:ncpu
    if ncpu >= 100      
        str_file_num = num2str(j,'%03d');
    elseif ncpu >= 10
        str_file_num = num2str(j,'%02d');
    else
        str_file_num = num2str(j,'%02d');
    end  
    traj_filename=[output_dir,'/traj_file_',str_file_num,'.nc'];
    
    status = ncread(traj_filename,'exitcode');
    release = ncread(traj_filename,'releasedate');
    polygon = ncread(traj_filename,'releasepolygon');
    
    polygon2=[polygon2, polygon'];
    release2=[release2, release'];
    status2=[status2, status'];
    np(j)=size(release,1);
end
time = ncread(traj_filename,'time');
nstep = (max(release2)-min(release2))*86400/time(2)+35;
start_day = min(release2);
dt = time(2)/86400;
npt = sum(np);% Total number of particles

plon=NaN(1,npt);
plat=NaN(1,npt);
cnt=zeros(1,npt);
skip=24;
%% Plot particles
for i=1:skip:nstep

    t = start_day+dt*i;
    
    for j = 1:ncpu

        if ncpu >= 100      
            str_file_num = num2str(j,'%03d');
        elseif ncpu >= 10
            str_file_num = num2str(j,'%02d');
        else
            str_file_num = num2str(j,'%02d');
        end  
        traj_filename=[output_dir,'/traj_file_',str_file_num,'.nc'];
        
        for k=1:np(j)
            ip=sum(np(1:j-1))+k;
            if release2(ip) <= t && cnt(ip)+skip<=size(time,1)
                cnt(ip)=cnt(ip)+skip;
                plon(ip) = ncread(traj_filename,'lon',[cnt(ip) k], [1 1]);
                plat(ip) = ncread(traj_filename,'lat',[cnt(ip) k], [1 1]);
            else
                plon(ip)=NaN;
                plat(ip)=NaN;
            end
        end
    end

    set(h_scatter,'XData',plon,'YData',plat,'CData',polygon2)

%     title(datestr(datetime(t,'ConvertFrom','juliandate'),'yyyy-mm-dd HH:MM'), 'FontSize' , 12)
    title(datestr(datetime(t,'ConvertFrom','juliandate'),'yyyy-mm-dd'), 'FontSize' , 12)
   
    drawnow
    hgexport(figure(1), strcat('output/figs_png\flt01_',num2str(i,'%0.4u'),'.png'),hgexport('factorystyle'),'Format','png');

end
hold off;
