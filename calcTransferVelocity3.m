function [kw,Sc]=calcTransferVelocity3(era,K0)
%Use function to calculate the gas transfer velocity 
%Calculation is based on Reichl and Deike (2020)

%Kw is the sum of non-bubble (Kwnb) and bubble (Kwb) components:

%Kwnb=Anb*ustar*(Sc/660)^-(1/2)
%where Anb=1.55x10^-4 (Fairral et all (2011))
%      Sc=Schmidts number 

%and Kwb=(Ab/K0*R*T0)*ustar^(5/3)*(gHs)^(2/3)*(Sc/660)^(-1/2)
%where Ab=1 +/- 2 x10^-5 
%      R=ideal gas constant
%      T0=SST
%      g=gravity

%Anb=1.55e-04; Ab= 6.55315e-5; g=9.8;
%R=0.08206; %units=L.atm/K.mol
R=8.206e-5; %units=m3.atm/K.mol
g=9.8;

%find the valid SST and SSS for each SWH point
kw.kw=[]; kw.NB=[];  kw.B=[]; 
for i=1:size(era.swh,3)

    k=era.sst(:,:,i);
    kk=Sc(:,:,i);

    %Convert T0 from Celcius to Kelvin
    t0k=k + 273;


    %% calculate kwnb

   l=era.ustar(:,:,i).*((kk./660).^(-1/2));
    nb=l*360000; %convert m/s to cm/hr
   
    kw.NB(:,:,i)=nb;

    %% calculate kwb
    ko=K0.k0(:,:,i).*(1/0.001); %convert to mol.(m3.atm)-1
   
    % alpha=(ko.*R.*t0k)/101325;
    l=(1/(ko.*R.*t0k)); 
    % l=Ab./alpha;
    ll=(era.ustar(:,:,i)).^(5/3);
    l2=(g.*era.swh(:,:,i)).^(2/3);
    ll2=(kk./660).^(-1/2);

    j=l.*ll.*l2.*ll2;
  
    j=j*360000; %convert m/s to cm/hr
 

    kw.B(:,:,i)=j;
    
    %Calculate the total transfer velocity 
    %kw=kw.nb+kw.b

    % j3=j+nb;
    % kw.kw(:,:,i)=j3;

end


end




