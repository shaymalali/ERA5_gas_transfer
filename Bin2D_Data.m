function [kmean_2d]=Bin2D_Data(kw,output,plot_flag)
%% Script to bin transfer velocity data on a 2x2 degree grid 
% First, calculate the climatological means for each bin
%Second, flip the data onto a 361x720 grid with the latitudes going
%northward (-90 to 90) and the longitudes going eastward (0-360)
%current data currently on a 720x360 grid 

%Inputs: kw=transfer velocity
%        output=OCIM output that provides OCIM grid
%        plot_flag= flag to create plots or not

%% calculate the climatological mean
[kmean]=calcClimateMean(kw,kw.lat,kw.lon,kw.dt);
%kmean=kw.kw;
%% Transform the data from 720x361 to 361x720
%Flip the data so latitudes range from [-90 90]

kw.lat2=flip(kw.lat);

kmean2=flip(kmean,2);
kmean2=pagetranspose(kmean2);


if plot_flag == 1
    figure(1)
    subplot(2,1,1)
    worldmap([-90 90], [-180 180])
    load coastlines
    % pcolorm(double(era5_monthly.lat),double(era5_monthly.lon),era5_monthly.swh(:,:,2)); hold on
    pcolorm(double(kw.lat),double(kw.lon),kmean(:,:,1)'); hold on
    geoshow(coastlat,coastlon,'Color','k','LineWidth',1); hold off
    h=colorbar;
    clim([0 60])
    ylabel(h,'Transfer Velocity (cm/hr)','FontSize',10);
    title('Climatological Mean BEFORE Transposing and Flipping')
    subplot(2,1,2)
    worldmap([-90 90], [-180 180])
    load coastlines
    pcolorm(double(kw.lat),double(kw.lon),kw.kw(:,:,2)); hold on
    % pcolorm(double(grd.yt),double(grd.xt),squeeze(kmean2_2d(1,:,:))); hold on
    geoshow(coastlat,coastlon,'Color','k','LineWidth',1); hold off
    h=colorbar;
    clim([0 60])
    ylabel(h,'Transfer Velocity (cm/hr)','FontSize',10);
    title('Climatological Mean AFTER Transposing and Flipping')
end


%% Interpolate Data so its on a 91x180 grid (2x2deg)
% use OCIM grid to grid the data
%load CTL.mat            % load the OCIM1 matrix transport model
grd = output.grid;      % OCIM1 grid definition
msk = output.M3d;  

lons=kw.lon;
lats=kw.lat2;

%dx=0.5;
%copy first and last longitude to top and end to respect peridocity
lons2=[];
lons=kw.lon;
lons2=NaN(722,1);
lons2(1)=-0.5; % x(1) - dx 
lons2(end)=360.5; % x(end) + dx
lons2(2:721)=lons;

[X,Y]=meshgrid(lons2,lats);
for i=1:size(kmean2,3)

    d=kmean2(:,:,i);
    % interpolate the data
    X=double(X);
    Y=double(Y);
    lats=double(lats);
    lons=double(lons);
    d2_interp=interp2(lons,lats,d,X,Y);

    % Use inpaint_nans to get rid of the nans
    test=inpaint_nans(d2_interp,2);
    
    %interpolate data onto 2deg grid
    
    mu1=interp2(X,Y,test,grd.XT,grd.YT);

    %apply land mask
    idx=find(msk(:,:,1)==1);

    mu2=NaN(size(grd.XT));
    mu2(idx)=mu1(idx);

    kmean_2d(:,:,i)=mu2;
    
    clear mu1 mu2 var1 n1 Q1 test d d2_interp
end

   
