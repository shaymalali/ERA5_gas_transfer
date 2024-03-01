function [K0]=calcSolubility3(era,SSS)
%use function to calculate the solubilty at every significant wave height
%from ERA-5 reanalysis data

%Solubility Equation is derived from Reichle and Deike (2020)

%ln(K0)=a1+a2(100/T0)+a3*ln(T0/100)+S0(b1+b2(T0/100)+b3(t0/100)^2)
%a1 a2 a3 b1 b2 b2 are coefficients
%T0 is the sea surface salinity in Kelvin 
%S0 is the surface salinity in g/kg

%write out coefficients
a1=-58.0931; a2=90.5069; a3=22.2940;
b1=0.027766; b2=-0.025888; b3=0.0050578; 

% a1=-60.2409; a2=93.4517; a3=23.3585;
% b1=0.023517; b2=-0.023656; b3=0.0047036;

K0.k0=[]; K0.lat=[]; K0.lon=[];  K0.SST=[]; K0.sss=[];

for i=1:size(era.swh,3)
    
    %get data for the month and year of each data set


    month_era=month(era.dt(i));

    sss_e=sss.sss(:,:,month_era);

% convert SST from celcius to Kelvin 

     sst_k=era.sst(:,:,i) +273;

     j1=a1;
     j2=(a2.*(100./sst_k));
     j3=(a3.*log(sst_k./100));
     j4=(b1+(b2.*(sst_k./100)))+(b3.*((sst_k./100).^2));
     j5=sss_e.*j4;
     %jj=a1+(a2.*(100./T0))+(a3.*log(T0./100))+(S0.*(b1+b2.*(T0./100))+(b3.*((T0./100).^2)));
     jj=j1+j2+j3+j5;
     jj=exp(jj);

     K0.k0(:,:,i)=jj;
end





