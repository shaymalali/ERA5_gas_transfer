function data = getglodap_c14(t0,tf)
% function data = getglodap_c14(t0,tf)
%     data: required, struct  
%       vn: required, variable name for extract
%       t0: optional, starting year for data extraction
%       tf: optional, ending year for data extraction
% extract variable from GLODAPv2.2021 data
% interpolate data to OCIM grid
% build an operator to map OCIM data to where
% GLODAP data available
    
%load /DFS-L/DATA/primeau/salali/CTL_OCIM_gas_exchange_OCT2023/DATA/M3d91x180x24.mat M3d grd;
    
    load M3d91x180x24.mat M3d grd
    iwet = find(M3d==1); nwet = length(iwet);
    
    % load the GLODAPv2 data
    disp(sprintf('now loading %s from GLODAPv2.2022','G2c14')); 
    
    %dpath= '/DFS-L/DATA/primeau/salali/CTL_OCIM_gas_exchange_OCT2023/DATA';
    dpath = '.';
    fpath=sprintf('%s/GLODAPv2.2022_Merged_Master_File.mat',dpath);
    load(fpath,'G2latitude','G2longitude','G2pressure','G2year','G2pressure','G2year','G2month','G2depth','G2c14','G2c14f');
    
    lat=G2latitude;
    lon=G2longitude;
    pressure=G2pressure;
    year=G2year;
    month=G2month;
    depth=G2depth;
    %fix negative longitudes
    ineg=find(lon<0);
    lon(ineg)=360+lon(ineg);
    
    
    disp(sprintf('extracting %s from %i to %i','G2c14',t0,tf));
    c14  = G2c14;
    c14f = G2c14f;
    
    c14err = max(0.1*c14,0.05);% NOT CORRECT
    idat = find( ( c14 ~= -9999 ) & ( depth >= 0 ) & ( c14f==2 ) & ( year >= t0 ) & ( year <= tf) );
    c14 = c14(idat);
    c14err = c14err(idat);
    
    % remove outliers of GLODAPv2 data
    %[~,TF] = rmoutliers(c14,'percentiles',[0.1 99.9]);
    %ibadV = find(TF==1);
    
    
    % outlier removed
    lon = lon(idat);
    lat = lat(idat);
    z = depth(idat);
    yr = year(idat);
    mo = month(idat);
    idx = yr * 12 + mo; 
    ID = unique(idx); 
    
    % bin to grd (year,month)
    for i = 1:length(ID)
        id = find( idx == ID(i) );
        
        [mu,vard,n] = bin3d(lon(id), lat(id), z(id), c14(id), c14err(id), grd.XT3d, grd.YT3d, grd.ZT3d);  
        
        data.c14star(:,i) = mu(iocn);
        data.varc14(:,i) = vard(iocn);
        data.nc14(:,i) = n(iocn);
    end
    data.c14id = ID;
    data.c14z = z;

    %% create operator to compare the monthly GLODAPv2 data and model result
    id1 = data.c14id;
    II1 = ones(12,1); % months
    II2 = ones(tf-t0+1,1); % years 
    idd1 = kron(II1,[t0:tf]'*12) + kron([1:12]',II2);
    l1 = length(id1);
    l2 = length(idd1);
    H = sparse(l1,l2); 
    for i = 1:l1
        j = find(idd1 == id1(i));
        H(i,j) = 1;
    end  
    H1 = kron(H,speye(nwet)); 

    % step2: find grd
    c14star = data.c14star;
    d0 = @(x) spdiags(x(:),[0],length(x(:)),length(x(:)));
    H2 = d0(c14star(:)~=-9999);
    ikeep = find(c14star(:)~=-9999);
    
    H2 = H2(:,ikeep)'; 
    data.C14h1 = H1;
    data.C14h2 = H2;
    
end
