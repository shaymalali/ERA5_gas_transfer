%% Script to concatenate and bin calculated Transfer Velocity (kw) values 
% current kw values calculated using ERA hourly SWH and friction velocity values between 1990 and 2023
% Current grid is a 0.5x0.5 grid (361x721) 
% Binning to OCIM grid (2x2)

%load OCIM grid 
load CTL.mat 
grd=output.grid;
M3d=output.M3d;


%concatenate the kw values from each decade 
flag=0; %0 to concatenate the transfer velocities using Reichl and Deike (2020) parametrization %1 to concatenate Wanninkhof(2014) parameterization

if flag == 0
	direc_90s=dir('ERA_Hour_Data1990*');
	direc_2000s=dir('ERA_Hour_Data2000s*');
	direc_2010s=dir('ERA_Hour_Data2010s*');
      
       
        %keyboard;  %different fraction # used in each kw file
	for i=14:length(direc_90s)
   
   	load(direc_90s(i).name)
    	kw_90=kw_monthmean;
        a=kw_90.kwnb;
        a2=kw_90.kwb;
        a3=kw_90.wan;
       

    	load(direc_2000s(i).name)
    	kw_2000=kw_monthmean;
    	b=kw_2000.kwnb;
        b2=kw_2000.kwb;
        b3=kw_2000.wan;
       

    	load(direc_2010s(i).name)
    	kw_2010=kw_monthmean;
    	c=kw_2010.kwnb;
        c2=kw_2010.kwb;
        c3=kw_2010.wan;

    %concatenate the transfer velocities
    %kw.kw=cat(3,kw_90.kw,kw_2000.kw,kw_2010.kw);
    	kw.NB=cat(3,a,b,c);
        kw.B=cat(3,a2,b2,c2);
        kw.wan=cat(3,a3,b3,c3);
    	kw.lat=kw_90.lat;
    	kw.lon=kw_90.lon;
    	kw.dt=cat(2,kw_90.dt,kw_2000.dt,kw_2010.dt);

    	output.grid=grd;
    	output.M3d=M3d;

    	[kmean_2d_NB]=Bin2D_Data(kw.NB,output,0);
        [kmean_2d_B]=Bin2D_Data(kw.B,outout,0);
        [kmean_2d_wan]=Bin2D_Data(kw.wan,outout,0);
    
    %d=frac(i)
    	outname=sprintf('kmean2D_RD.mat');
    	save(outname,'kmean_2d_NB','kmean_2d_B','kmean_2d_wan','-v7.3');
    	%print('Done processing kw with fraction num %1.1f',frac(i));
        
            if i == 14
               sprintf('breaking now')
               break
            end 
	end


end  
