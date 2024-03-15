%% Process_Results.m 
%Author: Shayma Al Ali 3/14/2024

%Use script to process the results from the OCIM model, compare with
%GLODAP observations, calculate the Root Mean Square Error, 
% and optimize the model

%Step 1: Read in GLODAP Observations

%Step 2: Read in the results and calculate D14C

%Step 3: Calculate the RMSE 

%Step 4: Optimize the model

%% Step 1: Read in GLODAP Observations
t0=1750;
tf=2021;

data = getglodap_c14(t0,tf);

c14star=data.c14star;
H1=data.C14h1;
H2=data.C14h2;

Obs= H2*c14star(:); %find the observations that are in the OCIM grid
%% Step 2: Read in the results and calculate D14C
%Model results will be stored in the output directory with a file for each
%year calculated
%so if your model runs from 1750 to 2021, expect 272 files in the output
%directory 
%This script will parse through each file, read in the results, calculate
%the D14C for each year, and store all the results in one array 

%

% Set path where the output is stored
% each model run with a different dimensional fitting coefficient has its
% own model output directory
Path_02='./ModelOutput/Wan0.20'; %coefficient 0.20
Path_025='./ModelOutput/Wan0.25'; %coefficient 0.25
Path_030='./ModelOutput/Wan0.30'; %coefficient 0.30
Path_039='./ModelOutput/Wan0.39'; %coefficient 0.39
Path_045='./ModelOutput/Wan0.45'; %coefficient 0.45
Path_06='./ModelOutput/Wan0.60'; %coefficient 0.60
Path_09='./ModelOutput/Wan0.9'; %coefficient 0.9
Path_12='./ModelOutput/Wan1.2' %coefficient 1.2

[d14c_02]=load_results(Path_02);
Model_02=H2*H1*d14c_10; %find the data that exists in the same time and place as the GLODAP observations

[d14c_025]=load_results(Path_025);
Model_025=H2*H1*d14c_025; %find the data that exists in the same time and place as the GLODAP observations

[d14c_030]=load_results(Path_030);
Model_030=H2*H1*d14c_030;

[d14c_039]=load_results(Path_039);
Model_039=H2*H1*d14c_039;

[d14c_045]=load_results(Path_045);
Model_045=H2*H1*d14c_045;

[d14c_06]=load_results(Path_06);
Model_06=H2*H1*d14c_06;

[d14c_09]=load_results(Path_09);
Model_09=H2*H1*d14c_09;

[Path_12]=load_results(Path_12);
Model_12=H2*H1*d14c_12;

%% Step 3: Calculate the RMSE
RMSE_02=rmse(Model_02,Obs);
RMSE_025=rmse(Model_025,Obs);
RMSE_030=rmse(Model_030,Obs);
RMSE_039=rmse(Model_039,Obs);
RMSE_045=rmse(Model_045,Obs);
RMSE_06=rmse(Model_06,Obs);
RMSE_09=rmse(Model_09,Obs);
RMSE_12=rmse(Model_12,Obs);

%% Optimize the model

x=[0.20;0.251;0.30;0.39;0.45;0.6;0.9;1.2];
y=[RMSE_02;RMSE_025;RMSE_030;RMSE_039;RMSE_045;RMSE_06;RMSE_09;RMSE_12];

[alpha,min_x]=solveparabola(x,y,0);
x2=linspace(1e-6,1e-4,1000);
y2=polyval(alpha,x2);

%find min point of parabola
f=@(x) (alpha(1)*x^2)+(alpha(2)*x)+alpha(3)
initial=0;
min_x=fminsearch(f,initial)
min_y=f(min_x)

figure(1)
plot(x,y,'or','MarkerSize',10); hold on
plot(min_x,min_y,'ob','MarkerSize',10); hold on
plot(x2,y2,'-.k'); hold off
xlim([1e-6 1e-4]);
%ylim([50 250])
grid on
xlabel('Dimensional Fitting Coefficients')
ylabel('RMSE')
legend('Root Mean Square Error','Minimum Point','Parabola')
title('Model Optimization')





