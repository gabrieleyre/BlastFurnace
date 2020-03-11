% Read initial Chemical and Reaction Data from .csv files:

ChemData = readtable('/Users/gabrieleyre/Dropbox/Warwick/ES327/MATLAB/chemdata.csv');
% Columns: 1:Element,2:Name,3:Density,4:MMass,5:Mvol,6:stdHForm,7:Cpn
% Units: Density:[kg/m^3],MMass:[kg/mol],Mvol:[m^3/mol],stdHForm:[J/mol],Cpn:[J/mol·K]
% Rows: 1:C,2:O2,3:CO,4:CO2,5:Fe2O3,6:Fe3O4,7:FeO,8:Fe,9:N2,10:CaCO3,11:CaO,12:SiO2,13:CaSiO3,14:Air

ReactionData = readtable('/Users/gabrieleyre/Dropbox/Warwick/ES327/MATLAB/reactdata.csv');
% Columns: 1:Reaction,2:Description,3:f_factor,4:a_energy
% Units: f_factor:[s^-1],a_energy:[J]

% Set initial parameters:

tspan = 50000; % Time span for simulation
CokeOreRatio = 3; % Coke : Ore molar ratios to investigate
HearthDiameter = 20; % Diameter of furnace [m]
HearthHeight = 10; % Height of furnace [m]
WallThickness = 2; % Thickness of BF walls [m]
ThermalConductivity = 30; % W/mK. Assuming Corundrum (Al2O3).
HotBlastRate = 60; % Volumetric flow rate of air [m^3/s]
HotBlastTemp = 1200; % Hot Blast Air Temperature [°C]
OxygenEnrichment = 3; % Assumes 20.95% of O2 in dry air. [%]
Temp_ext = 40; % Temperature outside Blast Furnace, assumed constant. [°C]
Temp_initial = 40; % Initial Blast Furnace temperature at start of model [°C]
SilicaImpurity = 2; % Silica content in Fe2O3[%]
nZones = 10; % Number of zones to model. Must be at least 3.
Res = 1; % Size of each time step for  modelling [s]. Lower is more accurate but requires more power.

% Pre-define arrays for efficiency:

[RateCoef,Rate,Hr,SpeciesVolume] = deal(zeros([tspan/Res height(ReactionData) nZones]));
[nTotal,CpnAv,VolUsed,VolUsedFraction,VolAvailable,VolAvailableFraction,HrNet] = deal(zeros([tspan/Res nZones]));
[dTdt,dhWall,dhProducts,dhBurdenDescentNet,dhGasAscentNet,dhHotBlast,HNet] = deal(zeros([tspan/Res nZones]));
[FeedRate,dnOut] = deal(zeros([tspan/Res 13]));
[dnBurdenDescent,dnGasAscent] = deal(zeros([tspan/Res 13 nZones]));
[VolUsedNew,VolAvailableBurden] = deal(zeros([(tspan/Res) 1]));
n = zeros([(tspan/Res)+1 13 nZones]);
dn = zeros([tspan/Res 13 nZones]);
Temp = zeros([(tspan/Res)+1 nZones]);

% Initialise Blast Furnace simulation:

for j = CokeOreRatio
    [RateCoef,Rate,Hr,SpeciesVolume,nTotal,CpnAv,VolUsed,VolUsedFraction,VolAvailable,VolAvailableFraction,HrNet,dTdt,dhWall,dhProducts,dhBurdenDescentNet,dhGasAscentNet,dhHotBlast,HNet,FeedRate,dnOut,dnBurdenDescent,dnGasAscent,VolUsedNew,VolAvailableBurden,n,dn,Temp,t] = blastfurnace(ChemData,ReactionData,tspan,j,HearthDiameter,HearthHeight,WallThickness,ThermalConductivity,HotBlastRate,HotBlastTemp,OxygenEnrichment,Temp_ext,Temp_initial,SilicaImpurity,nZones,Res);
end

% Tata charge in terms of tonnes. Express feedrate in terms of that.

% For reduction rate, use grey colormap.

% Plot stuff

%imagesc(y,x,Temp');hold on;colorbar;colormap(hot);xlabel('Time (s)'),ylabel('Zone');title('Blast Furnace zone temperature vs. time'),ylabel(colorbar, 'Temperature (°C)');

% FINAL CALCULATIONS FOR PLOTTING

% Calculate Feed Rate in terms of mass:

FeedRateMass = (ChemData.Mmass(1:13).*FeedRate(:,:)')';

% Calculate percentage of reduction in each zone:

%PercentReduction(:,:) = (n(:,8,:)./(n(:,5,:)+n(:,6,:)+n(:,7,:)+n(:,8,:)))*100;

% Total percentage reduction is calculated in final column of PercentReduction:

%PercentReduction(:,nZones+1) = sum(PercentReduction,2)/nZones;

% Calculate dhTopGas = dhGasAscent from top zone = 1:

dhTopGas = dhGasAscentNet(:,1);

% Calculate average blast furnace temperature:

TempAvg = sum(Temp,2)/nZones;

% Cumulative Heat calculations:

HrNetCumulative = cumsum(HrNet);
dhWallCumulative = cumsum(dhWall);
dhProductsCumulative = cumsum(dhProducts(:,nZones));
dhTopGasCumulative = cumsum(dhTopGas(:,1));
dhHotBlastCumulative = cumsum(dhHotBlast(:,nZones));
HNetCumulative = cumsum(HNet);
FeoutCumulative = cumsum(dnOut(:,8));
SlagoutCumulative = cumsum(dnOut(:,13));
