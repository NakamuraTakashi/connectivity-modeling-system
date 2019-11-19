clear all
%% parameters
YEAR = '1996'
polyfile = 'D:/Documents/GIS_data/Yasuda_COTS/habitat.shp';

% nest_dir='../expt/expt_cots_2002/nests';

% start_date = '2002/4/1';
% end_date = '2002/11/1';

% polygonfile='polygonFile';
% relfile='releaseFile';

nest_dir=['../expt/expt_cots_' YEAR '/nests'];

start_date = [YEAR '/4/1'];
end_date = [YEAR '/11/1'];

polygonfile=['../expt/input_cots_' YEAR '/polygonFile'];
relfile=['../expt/input_cots_' YEAR '/releaseFile'];

n_particle_per_polygon = 100;
temp_threshold = 28; % (oC) Temperature threshold of COST spawing
tot_days = 30; % Total days for spawning
depth = 0;

%% 
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
%% read in the nest data
t = datetime(start_date);

nest_filename=[nest_dir, '/nest_1_', datestr(t,'yyyymmddHHMMss'), '.nc'];

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

%% draw polygons
colors=jet(n_poly);

for i=1:n_poly
    pgon = polyshape(S(i).Lon, S(i).Lat);
    plot(pgon,'FaceColor',colors(S(i).id,:),'FaceAlpha',0.5);
    text(mean(S(i).BoundingBox(:,1))-0.7, mean(S(i).BoundingBox(:,2)+0.2), num2str(S(i).id),'FontSize',11);
end
drawnow();

%% Write polygon file
fid = fopen(polygonfile,'w');
format = '%f %f %i\n';
for i=1:n_poly
    for j=1:size(S(i).Lat,2)
        fprintf(fid,format, S(i).Lon(j), S(i).Lat(j),S(i).id );
    end
end

fclose(fid);
    
%% Pick up the particle release points

for i=1:n_poly
    ip = 0;
    iloc = 0;
    lonb = S(i).BoundingBox(:,1);
    latb = S(i).BoundingBox(:,2);   
  
    while ip < n_particle_per_polygon

        lon = lonb(1) + (lonb(2)-lonb(1))*rand;
        lat = latb(1) + (latb(2)-latb(1))*rand;
        
         [M,id_lon] = min(abs(lon-lonAxis));
         [M,id_lat] = min(abs(lat-latAxis));
        for j=2:size(lonAxis,1)
            if lonAxis(j)>lon
                id_lon1=j-1;
                id_lon2=j;
                break
            end
        end
        for j=2:size(latAxis,1)
            if latAxis(j)>lat
                id_lat1=j-1;
                id_lat2=j;
                break
            end
        end
        
      
        if inpolygon(lon, lat, S(i).Lon, S(i).Lat) && ...
            uvel(id_lon1, id_lat1)<2^100 && uvel(id_lon2, id_lat1)<2^100 && ...
            uvel(id_lon1, id_lat2)<2^100 && uvel(id_lon2, id_lat2)<2^100
            
            ip = ip+1;
            Poly(i).ID(ip)=ip;
            Poly(i).pLat(ip)=lat;
            Poly(i).pLon(ip)=lon;
            if ip == 1
                iloc = iloc+1;
                Poly(i).N_inpoly=iloc;
                Poly(i).iLat_inpoly(iloc)=id_lat;
                Poly(i).iLon_inpoly(iloc)=id_lon;
            elseif isempty( find(Poly(i).iLat_inpoly==id_lat)) && isempty( find(Poly(i).iLon_inpoly==id_lon))
                iloc = iloc+1;
                Poly(i).N_inpoly=iloc;
                Poly(i).iLat_inpoly(iloc)=id_lat;
                Poly(i).iLon_inpoly(iloc)=id_lon;           
            end       
        end        
    end
end

%% Extract time-series temperature data at each polygon

for i=1:n_poly
    it = 1;
    st = 0;
    t = datetime(start_date);    

    while t <= datetime(end_date)
        nest_filename=[nest_dir, '/nest_1_', datestr(t,'yyyymmddHHMMss'), '.nc'];
        
        for iloc=1:Poly(i).N_inpoly
            try
                Poly(i).temp(iloc) = ncread(nest_filename,'zt',[Poly(i).iLon_inpoly(iloc) Poly(i).iLat_inpoly(iloc) 1 1], [1 1 1 1]);
            catch
                st = 1;
                disp(i)
                disp(nest_filename)
                break
            end
        end
        if st == 0
            Poly(i).mean_temp(it) = mean(Poly(i).temp);
        else
            Poly(i).mean_temp(it) = Poly(i).mean_temp(it-1);
        end
        time(it) = t;

        it = it+1; % update time counter
        t = t+1; % update date
    end
    [Poly(i).max_temp, Poly(i).id_max]  = max(Poly(i).mean_temp);
end

%% 

fid = fopen(relfile,'w');
        %Polygon Longitude Latitude Depth Number Year Month Day Second
format = '%i %f %f %i %i %i %i %i %i\n';

n_day=size(Poly(1).mean_temp,2);

for i=1:n_poly
    idate=0;
    Poly(i).reldays=0;
    fspawn = 0;
    if Poly(i).max_temp >= temp_threshold
        for j=1:n_day       
            if Poly(i).mean_temp(j) >= temp_threshold
                fspawn = 1;
            end
            if Poly(i).mean_temp(j) >= temp_threshold-2 && fspawn == 1
                for ip=1:n_particle_per_polygon
                    fprintf(fid,format, S(i).id,Poly(i).pLon(ip), Poly(i).pLat(ip),depth,1,year(time(j)),month(time(j)),day(time(j)),0 );
                end
                idate=idate+1;
                Poly(i).reldays = idate;
                if idate >= tot_days
                    break
                end
            end      
        end
    else
        for j=Poly(i).id_max:n_day       
            if Poly(i).mean_temp(j) >= temp_threshold-2
                for ip=1:n_particle_per_polygon
                    fprintf(fid,format, S(i).id,Poly(i).pLon(ip), Poly(i).pLat(ip),depth,1,year(time(j)),month(time(j)),day(time(j)),0 );
                end
                idate=idate+1;
                Poly(i).reldays = idate;
                if idate >= tot_days
                    break
                end
            end
        end       
    end   
end

fclose(fid);


