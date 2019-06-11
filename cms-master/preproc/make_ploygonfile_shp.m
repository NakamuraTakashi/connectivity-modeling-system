%% parameters
nest_filename='../expt/expt_japan/nests/nest_1_20020601000000.nc';
polyfile = 'D:/Documents/GIS_data/Yasuda_COTS/habitat.shp';

polygonfile='polygonFile';
relfile='releaseFile';
n_particle_per_polygon = 100;
depth = 0;

%% 
S = shaperead(polyfile);
n_poly = size(S,1);
for i=1:n_poly
    S(i).Lat = rmmissing(S(i).Y);
    S(i).Lon = rmmissing(S(i).X);
end


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
%% 
colors=jet(n);

for i=1:n_poly
    pgon = polyshape(S(i).Lon, S(i).Lat);
    plot(pgon,'FaceColor',colors(i,:),'FaceAlpha',0.3);
    text(mean(S(i).BoundingBox(:,1))-0.1, mean(S(i).BoundingBox(:,2)+0.1), num2str(S(i).id));
end


%%
fid = fopen(polygonfile,'w');
format = '%f %f %i\n';
for i=1:n_poly
    for j=1:size(S(i).Lat,2)
        fprintf(fid,format, S(i).Lon(j), S(i).Lat(j),i );
    end
end

fclose(fid);
    
%% 
fid = fopen(relfile,'w');
        %Polygon Longitude Latitude Depth Number Year Month Day Second
format = '%i %f %f %i %i %i %i %i %i\n';

for i=1:n_poly
    ip = 0;
    lonb = S(i).BoundingBox(:,1);
    latb = S(i).BoundingBox(:,2);   
    while ip < n_particle_per_polygon

        lon = lonb(1) + (lonb(2)-lonb(1))*rand;
        lat = latb(1) + (latb(2)-latb(1))*rand;
        
%         [M,id_lon] = min(abs(lon-lonAxis));
%         [M,id_lat] = min(abs(lat-latAxis));
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
                    
            fprintf(fid,format, i,lon, lat,depth,1,2002,6,1,0 );
            ip = ip+1;
        end
        
    end
end
fclose(fid);