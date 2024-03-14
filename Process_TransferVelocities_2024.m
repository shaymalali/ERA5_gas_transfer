%% Transfer Velocity Calculation 
% Use script to download ERA5 reanalysis data to calculate the transfer
% velocity (Kw) using 2 different parameterizations 

% Transfer Velocity 1: Reichle and Deike (2020) Parameterization that
% includes bubbles and breaking waves and has two components: The
% Non-Bubble and Bubble Components (Kw=KwNB+KwB)
%Inputs: Significant Wave Height, Friction Velocity, Sea Surface
%Temperature, and Salinity

% Transfer Velocity 2: Wanninkhof (2014) Parameterization that includes the
% a wind speed dependence 
%Inputs: Wind Speed and Sea Surface Temperature
%

%% STEP ONE: Download ERA5 DATA
% Data is downloaded using a bash (.sh) for each variable
% (U10,V10,SWH,U*,SST)
% Only the U and V components of wind speed are available with ERA5 Hourly
% Data 

%Step One: Create an account with ECMWF
%Step Two: Create a cdsapirc file with Account # and Key 
%Step Three: chmod +x file.sh 
%Step Four: run file using ./file.sh 

%Each variable should have a netcdf file for each month of every year
%specified in the bash file 
%this script uses the years 1990-2022

%In this project, each variable has its own file for more download
%efficiency and speed, however you can download the variables to be stored
%in the same file
%However, the ERA5 atmospheric data is on a different grid than the ocean
%data and downloading any ocean data in the same file as atmospheric data
%will just push the ocean grid on the atmospheric grid without
%interpolation which will mess up the data so try to group the atmospheric
%and ocean data seperately

%% STEP TWO: Load in Salinity Data
% Salinity Data is a climatological mean of the Sea Surface Salinity from
% the NASA OISSS mission 
% The salinity is then interpolated onto the same grid as the ERA5 Ocean
% grid

load ('SSS.mat')

%% STEP 3: Read in the DATA 

%Load in the ERA5 files by creating a directory for each variable
direc_swh=dir('era_swh*');
direc_ustar=dir('era_ustar');
direc_sst=dir('era_sst');
direc_uwind=dir('era_u10*');
direc_vwind=dir('era_vwind');

%use a for loop to process each file separately 
for i=1:length(direc_swh)
    %The ERA5 atmospheric data is on a separate grid than the ERA5 ocean
    %data therefore we will group all of the ocean data into the same
    %structure while grouping the atmospheric data into its own structure

    %We will then interpolate the atmospheric data to be on the same grid
    %as the ocean data

    %load in the ocean data
    fname_swh=direc_swh(i).name;
    fname_sst=direc_sst(i).name;
    era_ocean.swh=ncread(fname_swh,'swh');
    era_ocean.lat=ncread(fname_swh,'latitude');
    era_ocean.lon=ncread(fname_swh,'longitude');
    era_atmo.sst=ncread(fname_sst,'sst'); %ocean product but data is in atmospheric grid
    era_ocean.time=ncread(fname_swh,'time');
    %load in the atmospheric data
    fname_ustar=direc_ustar(i).name;
    fname_uwind=direc_uwind(i).name;
    fname_vwind=direc_vwind(i).name;
    era_atmo.ustar=ncread(fname_ustar,'zust');
    era_atmo.u=ncread(fname_uwind,'u10');
    era_atmo.v=ncread(fname_vwind,'v10');
    era_atmo.lat=ncread(fname_ustar,'latitude');
    era_atmo.lon=ncread(fname_ustar,'longitude');
    era_atmo.time=ncread(fname_ustar,'time');

    %calculate the 10 m wind speed using the U and V components
    %Where Wspd=sqrt(U^2+V^2);
    era5_atmo.wspd=sqrt((era5_atmo.u).^2+(era5_atmo.v).^2);

    %interpolate the atmosphere data so its on the same grid as the ocean
    %data
    [X,Y]=meshgrid(era_ocean.lat,era_ocean.lon); %create a meshgrid of the data
    % ustar
    for j=1:size(era_atmo.ustar,3)
        V=era_atmo.ustar(:,:,j);
        V2=era_atmo.wspd(:,:,j);
        V3=era_atmo.sst(:,:,j);
        era_ocean.ustar(:,:,j)=interp2(era_atmo.lat,era_atmo.lon,V,X,Y);
        era_ocean.wspd(:,:,j)=interp2(era_atmo.lat,era_atmo.lon,V2,X,Y);
        era_ocean.sst(:,:,j)=interp2(era_atmo.lat,era_atmo.lon,V3,X,Y);
    end

    %apply land mask to ustar and u10
    era_ocean.ustar(isnan(era_ocean.swh))=NaN;
    era_ocean.wspd(isnan(era_ocean.swh))=NaN;
    era_ocean.sst(isnan(era_ocean.swh))=NaN;

    %calculate the solubility for the Reichl and Deike (2020)
    %Parameterizations
    [K0]=calcSolubility3(era_ocean,SSS);

    %Calculate the Schmidt Number (Sc) for both parameterizations
    %Sc is calculated using the Wanninkhof(2014) calculations 
    %Sc=A+Bt+Ct^2+Dt^3+Et^4 where t is SST in C
    for j=1:size(era_ocean.swh,3)
       t=era_ocean.sst(:,:,j);

       A= 2116.8;
       B=(-136.25).*(t);
       C=4.7353*((t).^2);
       D=(-0.092307)*((t).^3);
       E=0.000755*((t).^4);
       
       Sc(:,:,j)=A+B+C+D+E;
    end

    %calculate the bubble and non-bubble transfer velocities using the
    %Reichl & Deike (2020) Parameterizations

    %These calculations will NOT include the dimensional fitting
    %coefficients found in R&D (2020) paper but will be included in
    %further codes

    [Kw]=calcTransferVelocity3(era_ocean,K0,Sc);
    %Kw is a structure where:
    %Kw.NB is the non-bubble portion
    %Kw.B is the bubble portion

    %calculate the Wanninkhof (2014) transfer velocities
    %These calculations will NOT include the fitting coefficient 
    [Kw_Wan]=calcWankw(era_ocean,Sc);
     
    %calculate the monthly mean of the transfer velocities
    kw_monthmean.kwnb(:,:,i)=mean(kw.NB,3,"omitnan");
    kw_monthmean.kwb(:,:,i)=mean(kw.B,3,"omitnan");
    kw_monthmean.wan(:,:,i)=mean(Kw_Wan,3,"omitnan");
    kw_monthmean.lat=era_ocean.lat;
    kw_monthmean.lon=era_ocean.lon;
    kw_monthmean.dt(i)=datetime(year(era_ocean.dt(1)),month(era_ocean(1)),1);
  
    waveparams.meanswh(:,:,i)=mean(era_ocean.swh,3,"omitnan");
    waveparams.stdswh(:,:,i)=std(era_ocean.swh,0,3);
    waveparams.meanu(:,:,i)=mean(era_ocean.ustar,3,"omitnan");
    waveparams.stdu(:,:,i)=std(era_ocean.ustar,0,3);
    waveparams.meanu10(:,:,i)=mean(era_ocean.u10,3,"omitnan");
    waveparams.meanSc(:,:,i)=mean(Sc,3,"omitnan");

    clear era_ocea era_atmo Sc
end
save('ERA_Hour_Data1990.mat','kw_monthmean','waveparams');
