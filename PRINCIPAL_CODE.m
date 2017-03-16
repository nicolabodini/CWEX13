% _______________________________________________
% MULTIPLE WAKES DETECTION ALGORITHM
% based on wake detection for a row of 4 turbines
% (version for southwesterly wind)

% Author: Nicola Bodini - Sept 2016
% nicola.bodini@colorado.edu
% _______________________________________________

clear; clc; close all
warning off


% _______________________________________________
% INPUT DATA

% Load Windcube vertical profiling lidar data 
% represent upwind conditions
load WC49;

% Load Obukhov length data calculated from surface flux station measurements
% (a data - a row - for each second from the start of the day)
load Ldef;

% Load data from Windcube 200S scanning lidar
% and transform variables to double type
% Year UTC
year=ncread('WLS200S_2013-08-26.nc','year');
year=double(year);
% Month UTC
month=ncread('WLS200S_2013-08-26.nc','month');
month=double(month);
% Day UTC
day=ncread('WLS200S_2013-08-26.nc','day');
day=double(day);
% Hour UTC
hour=ncread('WLS200S_2013-08-26.nc','hour');
hour=double(hour);
% Minute UTC
minute=ncread('WLS200S_2013-08-26.nc','minute');
minute=double(minute);
% Second UTC
second=ncread('WLS200S_2013-08-26.nc','second');
second=double(second);
% Time since start of the day at the first instant of current scan (s)
stime=ncread('WLS200S_2013-08-26.nc','stime');
stime=double(stime);
% Scan type (1=PPI, 2=RHI)
stype=ncread('WLS200S_2013-08-26.nc','stype');
% Radial wind speed, i.e. line-of-sight velocity (m/s)
losv=ncread('WLS200S_2013-08-26.nc','rspd');
losv=double(losv);
losv=abs(losv); % positive values
losv=transpose(losv);
% Radial wind speed dispersion (m/s)
rdisp=ncread('WLS200S_2013-08-26.nc','rdisp');
rdisp=double(rdisp);
rdisp=abs(rdisp); % positive values
rdisp=transpose(rdisp);
% (Slant) range gates
range=ncread('WLS200S_2013-08-26.nc','rgate');
range=double(range);
range=transpose(range);
% Azimuth angle (deg)
az=ncread('WLS200S_2013-08-26.nc','cAngA');
az=double(az);
az = az*(pi/180); % convert azimuth angle from deg to rad
% Elevation angle (deg)
el=ncread('WLS200S_2013-08-26.nc','cAngE');
el=double(el);
el = el*(pi/180); % convert elevation angle from deg to rad
% Carrier-to-noise ratio (dB)
cnr=ncread('WLS200S_2013-08-26.nc','rCNR');
cnr=double(cnr);
cnr=transpose(cnr);

% _______________________________________________
% PRE-PROCESSING OF INPUT DATA

% Create a single date+time variable for WC49 profiling lidar data
tL=datenum(WC49(:,1),WC49(:,2),WC49(:,3),WC49(:,4),WC49(:,5),WC49(:,6));
% Create a single date+time variable for Windcube 200S scanning lidar data
decTime = datenum(year,month,day,hour,minute,second);

% Turbine rotor diameter [m]
D = 80;
% Express range gates in terms of rotor diameters D
range=range/D;

% Indices of range gates that will be used for wakes detection 
% (tot 98 horizontal ranges, from 100 m to 4950 m, every 50m)
rInd = 15:44; % from 800m to 2250m 
% Number of range gates that will be used for wakes detection
nR = numel(rInd);
% Number of models used for wake detection
nModels = 2; % no-wake model, 4 wakes model
% Initialize output table row index
iTable = 1;
% Plot settings
lws = 'LineWidth';
lw = 2;

% Eliminate RHI scans        
PPI = find(stype == 1);
decTime = decTime(PPI);
losv = losv(PPI,:);
rdisp=rdisp(PPI,:);
az = az(PPI,:);
el = el(PPI,:);
cnr = cnr(PPI,:);
stime = stime(PPI,:);

% Eliminate rows for which consecutive azimuth angles are the same
good = find(diff(az(:,1)) ~= 0);
decTime = decTime(good);
losv = losv(good,:);
rdisp=rdisp(good,:);
az = az(good,:);
el = el(good,:);
cnr = cnr(good,:);
stime = stime(good,:);

% Remove losv measurements whose corresponding CNR < -27 dB
losv(cnr < -27) = NaN;
rdisp(cnr < -27) = NaN;
cnr(cnr < -27) = NaN;

% Number of timestamps (m) and ranges (n)
[m,n] = size(losv);

% Calculate xG and yG coordinates
xG=zeros(m,n);
yG=zeros(m,n);
for i=1:m
    for j=1:n
        % longitudinal coordinate [D]
        xG(i,j) = range(1,j).*cos(el(i,1)).*cos(az(i,1));
        % transverse coordinate [D]
        yG(i,j) = range(1,j).*cos(el(i,1)).*sin(az(i,1));
    end
end

% Differences between adjacent elements of azimuth angle
dAz = diff(az(:,1));
% Flag to keep track of PPI scan direction
if dAz(1) > 0
   % clockwise scan
   flag = 1;
else
   % counterclockwise scan
   flag = 0;
end


% _______________________________________________
% FROM HERE: CONSIDER SINGLE PPI SCANS

% Initialize index that keeps track of time of beginning of scan
i1 = 1;

% Loop through each complete PPI scan
while i1 < m % m % m = number of timesteps in the considered day
     % i1 beginning of the scan; let's find the index i2 for the end
     % of the scan (i.e. when the scan changes its rotation sense)
     switch flag
           case 1 % if flag=1.. (i.e. if clockwise)
               i2 = find(dAz(i1:end) < 0,1,'first') + i1 - 1; % finds the first index where dAz<0
               flag = 0;
           case 0 % if flag=0.. (i.e. if counterclockwise)
               i2 = find(dAz(i1:end) > 0,1,'first') + i1 - 1;
               flag = 1;
      end

      % set i2 = m at the end of the file, for the last scan in the file
      if isempty(i2)
         i2 = m;
      end
      
      % skip VAD scans - scans to get wind profiles, large azimuth ranges
      if (max(az(i1:i2,1)) - min(az(i1:i2,1))) > 3.14
         i1=i2+1;
         continue % go to next iteration
      end
      
      % only proceed if in this PPI scan there are at least 25 azimuth angles
      if (i2 - i1) > 25
          
          % _______________________________________________
          % PROCESSING OF WINDCUBE49 PROFILING LIDAR DATA
          % WC49 measurements are taken every 2 minutes, 
          % while a WC200S PPI scan usually lasts for something like 100s.
          % So let's look for the WC49 data temporally closest to the scan.
          tDiff = abs(tL - decTime(i1));
          WCInd = find(tDiff == min(tDiff));
          % median WC49-measured wind direction during the scan @80m
          meanDir = WC49(WCInd,8);
          % median WC49-measured wind speed during the scan @80m
          meanSpd = WC49(WCInd,7);
          % _______________________________________________
          
          % _______________________________________________
          % PROCESSING OF OBUKHOV LENGTH DATA - for atmospheric stability
          % Obukhov length closer to the start of the actual PPI scan
          LInd = round(stime(i1));
          LObuk = Ldef(LInd);
          % _______________________________________________
          
          
          % Number of losv measurements in this PPI scan, per each range gate
          nT = i2 - i1 + 1;
          % Total number of measurements in the PPI scan
          nM = nT*nR;
          % Reshape predictors
          x = reshape(xG(i1:i2,rInd),nM,1);
          y = reshape(yG(i1:i2,rInd),nM,1);
          % Predictor array
          X = [x y];
          
          % Select losv, rdisp and cnr for the single PPI scan
          losvSelect = losv(i1:i2,rInd);
          rdispSelect = rdisp(i1:i2,rInd);
          rdispSelect = inpaint_nans(rdispSelect); % interpolates NaN values
          cnrSelect = cnr(i1:i2,rInd);
          % Reshape losv and cnr matrices
          losvReshape = reshape(losvSelect,nM,1);
          rdispReshape = reshape(rdispSelect,nM,1);
          cnrReshape = reshape(cnrSelect,nM,1);

          % Eliminate (NaNs) losv outliers using MAD method
          mu = nanmedian(losvReshape);
          sigma = 1.4826*mad(losvReshape,1); % http://it.mathworks.com/help/stats/mad.html
          losvSelect(losvSelect < mu-3*sigma | losvSelect > mu+3*sigma) = NaN;
          losvInterp = inpaint_nans(losvSelect); % interpolates NaN values
          losvReshape = reshape(losvInterp,nM,1);
                
          % Eliminate (NaNs) SNR outliers using MAD method
          muCNR = nanmedian(cnrReshape);
          sigmaCNR = 1.4826*mad(cnrReshape,1); % http://it.mathworks.com/help/stats/mad.html
          cnrSelect(cnrSelect < muCNR-3*sigmaCNR | cnrSelect > muCNR+3*sigmaCNR) = NaN;
          cnrInterp = inpaint_nans(cnrSelect); % interpolates NaN values
          cnrReshape = reshape(cnrInterp,nM,1);

          
          % _______________________________________________
          % ANALYSIS OF FREE-STREAM VARIABLES FOR SPATIAL TURBOLENCE INTENSITY
          % Freestream azimuth angles
          freeInd = 58:64; % freestream range indexes (upwind of the turbine)
          azFree = (180/pi)*az(i1:i2,1)';
          freeIndDim=freeInd(end)-freeInd(1)+1;
          azFree = repmat(azFree,freeIndDim,1);
          nFree = numel(azFree); % number of elements in azFree
          azFree = reshape(azFree,nFree,1);

          % Select and reshape freestream losv and snr matrices
          losvFree2D = losv(i1:i2,freeInd)';
          losvFree = reshape(losvFree2D,nFree,1);
          rdispFree2D = rdisp(i1:i2,freeInd)';
          cnrFree2D = cnr(i1:i2,freeInd)';
          cnrFree = reshape(cnrFree2D,nFree,1);

          % Eliminate freestream losv outliers using MAD method
          mFree = nanmedian(losvFree);
          sFree = 1.4826*mad(losvFree,1);
          losvFree2D(losvFree2D < mFree-3*sFree | losvFree2D > mFree+3*sFree) = NaN;
          losvInterpFree = inpaint_nans(losvFree2D);
          losvFree = reshape(losvInterpFree,nFree,1);
          rdispFree2D = 1./abs(inpaint_nans(rdispFree2D)); % weigths are inverse of losv dispersion
          rdispFree = reshape(rdispFree2D,nFree,1);

          % Eliminate freestream SNR outliers using MAD method
          mSNRfree = nanmedian(cnrFree);
          sSNRfree = 1.4826*mad(cnrFree,1);
          cnrFree2D(cnrFree2D < mSNRfree-3*sSNRfree | cnrFree2D > mSNRfree+3*sSNRfree) = NaN;
          snrInterpFree = inpaint_nans(cnrFree2D);
          cnrFree = reshape(snrInterpFree,nFree,1);

          % Ambient flow model
          % x indipendent variable
          % b(1) e b(2) model parameters
          mdlFree = @(b,x)b(1)*cosd(x - b(2));

          % Initial guess for parameters:
          % (median free-stream wind speed, turbine yaw angle)
          b0Free = [median(losvFree) (meanDir)];
          % Model fit, using azFree as independent variable, 
          % losvFree as dependent variable, 
          % b0Free as initial parameters, wFree as data weights
          fitFree = NonLinearModel.fit(azFree,losvFree,mdlFree,b0Free,'Weights',rdispFree);

          % Best-fit parameter estimates
          bHatFree = fitFree.Coefficients.Estimate';
          % Best-fit uLOS
          losvHatFree = mdlFree(bHatFree,azFree);
          % Root mean squared error
          rmse = fitFree.RMSE;
          
          % Spatial turbulence intensity [%]          
          spaceTI = rmse/bHatFree(1)*100;
          % _______________________________________________
          
          
          % _______________________________________________
          % HEART OF THE CODE: WAKES DETECTION
          % The wake detection is performed at each range gate
          
          % Number of fit tries (with different first-guess values, see below)
          ntent=4;
          
          for iR = rInd % loop across different range gates
              % initialize fit outputs
              % best-fit parameters
              bHat = cell(ntent,nModels);  % nModels=2 (no-wake, 4 wakes)
              % residuals
              res = NaN(nT,nModels); % nT = number of measurments per range gate
              % Jacobian
              J = cell(ntent,nModels);
              % covariance
              cov = cell(ntent,nModels);
              % mean squared error
              mse = NaN(ntent,nModels);
              % best-fit wind speeds
              losvHat = NaN(nT,nModels);
              % number of parameters in each model
              nP = [2 14];
              % sum of squares
              ss = NaN(ntent,nModels);
              % parameter confidence intervals
              ci = cell(ntent,nModels);
              % linear correlation coefficient (Pearson) 
              % between measured and modeled values of losv
              rho = NaN(ntent,nModels);

              % Select losv measurements at THIS particular range gate
              losvGate = losvInterp(:,iR - min(rInd) + 1);
              rdispGate = rdispSelect(:,iR - min(rInd) + 1);
              weightsGate = 1./rdispGate;  % weights: inverse of losv dispersion

              % Models for line-of-sight velocity
              % b wake model parameters
              % y transverse distance
              mdl = @(b,y)weightsGate.*PPImodel(b,y,range(1,iR));
              % model without weights - will be useful for some stathistics
              mdlnowgts = @(b,y)PPImodel(b,y,range(1,iR));
              
              % Transverse coordinates where measurements where performed 
              % at this range gate 
              yorig = yG(i1:i2,iR);
              
              % Identify indexes for first-guess peaks locations
              for itent=1:ntent % 4 fitting tries - with different first-guess
                  % ______________
                  if itent==1 % FIRST TRY
                      % First-guess transverse locations for the peaks of the wakes
                      % (this expressions come from some first fits)
                      ymu1=(-0.00000674868394838631*meanDir^2 + 0.00248035793126009000*meanDir - 0.228404567497937)*iR*D+ 0.01722404640348190000*meanDir^2 - 6.44726906449759*meanDir + 605.390889305966;
                      ymu2=(-0.0000156009515098093*meanDir^2+0.00611219086065158*meanDir-0.601592219288617)*iR*D+0.0328518719379387*meanDir^2-12.8743558019012*meanDir+1271.02403470688;
                      ymu3=(-0.00002114070267917320*meanDir^2+0.00831930228590515*meanDir-0.821484084356202)*iR*D+0.04228286106955580000*meanDir^2-16.5576639553901*meanDir+1634.78291290128;
                      % 1st wake
                      %Diff = abs(yorig-(ymu1-2.5));
                      %yind1 = find(Diff == min(Diff));
                      %Diff = abs(yorig-(ymu1+1));
                      %yind2 = find(Diff == min(Diff));
                      [mintemp,yind1] = (min(yorig));
                      [maxtemp,yind8] = max(yorig);
                      if flag==1 % clockwise PPI scan
                          yind2 = yind1 + floor(0.25*abs(yind8-yind1));
                          yind3 = yind2 + 1;
                          yind4 = yind3 + floor(0.25*abs(yind8-yind1));
                          yind5 = yind4 + 1;
                          yind6 = yind5 + floor(0.25*abs(yind8-yind1));
                          yind7 = yind6 + 1;
                      else % counter-clockwise PPI scan
                          yind2 = yind1 - floor(0.25*abs(yind8-yind1));
                          yind3 = yind2 - 1;
                          yind4 = yind3 - floor(0.25*abs(yind8-yind1));
                          yind5 = yind4 - 1;
                          yind6 = yind5 - floor(0.25*abs(yind8-yind1));
                          yind7 = yind6 - 1;
                      end
                      % 2nd wake
                      %Diff = abs(yorig-(ymu2-2));
                      %if flag==1 % clockwise PPI scan
                      %    yind3 = min(max(find(Diff == min(Diff)),yind2+4),length(yorig));
                      %else % counter-clockwise PPI scan
                      %    yind3 = max(min(find(Diff == min(Diff)),yind2-4),1);
                      %end
                      %Diff = abs(yorig-(ymu2+2));
                      %yind4 = find(Diff == min(Diff));
                      % 3d wake
                      %Diff = abs(yorig-(ymu3-2));
                      %if flag==1 % clockwise PPI scan
                      %    yind5 = min(max(find(Diff == min(Diff)),yind4+4),length(yorig));
                      %else % counter-clockwise PPI scan
                      %    yind5 = max(min(find(Diff == min(Diff)),yind4-4),1);
                      %end
                      %Diff = abs(yorig-(ymu3+2));
                      %yind6 = find(Diff == min(Diff));
                      % 4th wake
                      %Diff = abs(yorig - 10);
                      %if flag==1 % clockwise PPI scan
                      %    yind7 = min(max(find(Diff == min(Diff)),yind6+4),length(yorig));
                      %else % counter-clockwise PPI scan
                      %    yind7 = max(min(find(Diff == min(Diff)),yind6-4),1);
                      %end
                      %Diff = abs(yorig - 15.5);
                      %yind8 = find(Diff == min(Diff));
                      %[maxtemp,yind8] = max(yorig);
                  % ______________
                  elseif itent==2  % SECOND TRY
                      % 1st wake
                      Diff = abs(yorig-(ymu1+0.5));
                      yind1 = find(Diff == min(Diff));
                      Diff = abs(yorig-(ymu1+3.5));
                      yind2 = find(Diff == min(Diff));
                      % 2nd wake
                      Diff = abs(yorig-(ymu2+0.5));
                      if flag==1
                          yind3 = min(max(find(Diff == min(Diff)),yind2+4),length(yorig));
                      else
                          yind3 = max(min(find(Diff == min(Diff)),yind2-4),1);
                      end
                      Diff = abs(yorig-(ymu2+4.5));
                      yind4 = find(Diff == min(Diff));
                      % 3d wake
                      Diff = abs(yorig-(ymu3+0.5));
                      if flag==1
                          yind5 = min(max(find(Diff == min(Diff)),yind4+4),length(yorig));
                      else
                          yind5 = max(min(find(Diff == min(Diff)),yind4-4),1);
                      end 
                      Diff = abs(yorig-(ymu3+4.5));
                      yind6 = find(Diff == min(Diff));
                      % 4th wake
                      Diff = abs(yorig - 10);
                      if flag==1
                          yind7 = min(max(find(Diff == min(Diff)),yind6+4),length(yorig));
                      else
                          yind7 = max(min(find(Diff == min(Diff)),yind6-4),1);
                      end
                      Diff = abs(yorig - 15.5);
                      yind8 = find(Diff == min(Diff));
                  % ______________
                  elseif itent==3  % THIRD TRY
                      % 1st wake 
                      Diff = abs(yorig-(ymu1+3));
                      yind1 = find(Diff == min(Diff));
                      Diff = abs(yorig-(ymu1+6.5));
                      yind2 = find(Diff == min(Diff));
                      % 2nd wake
                      Diff = abs(yorig-(ymu2+3));
                      if flag==1
                          yind3 = min(max(find(Diff == min(Diff)),yind2+4),length(yorig));
                      else
                          yind3 = max(min(find(Diff == min(Diff)),yind2-4),1);
                      end
                      Diff = abs(yorig-(ymu2+7));
                      yind4 = find(Diff == min(Diff));
                      % 3d wake
                      Diff = abs(yorig-(ymu3+4));
                      if flag==1
                          yind5 = min(max(find(Diff == min(Diff)),yind4+4),length(yorig));
                      else
                          yind5 = max(min(find(Diff == min(Diff)),yind4-4),1);
                      end
                      Diff = abs(yorig-(ymu3+8));
                      yind6 = find(Diff == min(Diff));
                      % 4th wake
                      Diff = abs(yorig - 10);
                      if flag==1
                          yind7 = min(max(find(Diff == min(Diff)),yind6+4),length(yorig));
                      else
                          yind7 = max(min(find(Diff == min(Diff)),yind6-4),1);
                      end
                      Diff = abs(yorig - 15.5);
                      yind8 = find(Diff == min(Diff));
                  % ______________
                  else % FOURTH TRY
                      % 1st wake  
                      Diff = abs(yorig-(ymu1+6));
                      yind1 = find(Diff == min(Diff));
                      Diff = abs(yorig-(ymu1+10));
                      yind2 = find(Diff == min(Diff));
                      % 2nd wake
                      Diff = abs(yorig-(ymu2+8));
                      if flag==1
                          yind3 = min(max(find(Diff == min(Diff)),yind2+4),length(yorig));
                      else
                          yind3 = max(min(find(Diff == min(Diff)),yind2-4),1);
                      end
                      Diff = abs(yorig-(ymu2+12));
                      yind4 = find(Diff == min(Diff));
                      % 3d wake
                      Diff = abs(yorig-(ymu3+10));
                      if flag==1
                          yind5 = min(max(find(Diff == min(Diff)),yind4+4),length(yorig));
                      else
                          yind5 = max(min(find(Diff == min(Diff)),yind4-4),1);
                      end
                      Diff = abs(yorig-(ymu3+14));
                      yind6 = find(Diff == min(Diff));
                      % 4th wake
                      Diff = abs(yorig - 12);
                      if flag==1
                          yind7 = min(max(find(Diff == min(Diff)),yind6+4),length(yorig));
                      else
                          yind7 = max(min(find(Diff == min(Diff)),yind6-4),1);
                      end
                      Diff = abs(yorig - 17);
                      yind8 = find(Diff == min(Diff));
                  end
              
                  % First-guess estimates for model parameters
                  % Turbine yaw angle
                  phi0 = (meanDir-180)*(pi/180);
                  % Ambient wind speed
                  u0 = max(mu,b0Free(1));
                  % Velocity deficits amplitude
                  a0 = zeros(1,4);
                  a0(1,1) = u0 - min(losvGate(min(yind1,yind2):max(yind1,yind2),1));
                  a0(1,2) = u0 - min(losvGate(min(yind3,yind4):max(yind3,yind4),1));
                  a0(1,3) = u0 - min(losvGate(min(yind5,yind6):max(yind5,yind6),1));
                  a0(1,4) = u0 - min(losvGate(min(yind7,yind8):max(yind7,yind8),1));
                  % Positions of the center of the wakes
                  mu0 = [yorig(losvGate==min(losvGate(min(yind1,yind2):max(yind1,yind2),1)))'...
                         yorig(losvGate==min(losvGate(min(yind3,yind4):max(yind3,yind4),1)))'...
                         yorig(losvGate==min(losvGate(min(yind5,yind6):max(yind5,yind6),1)))'...
                         yorig(losvGate==min(losvGate(min(yind7,yind8):max(yind7,yind8),1)))'];
                  % b0 contains the initial estimates for model paramaters
                  b0 = cell(nModels);
                  % No-wake model
                  b0{1} = [phi0 u0];
                  % 4 wakes model
                  b0{2} = [b0{1} a0(1) a0(2) a0(3) a0(4) mu0(1) mu0(2) mu0(3) mu0(4) 0.3 0.3 0.3 0.3];


                  % Nonlinear regression
                  for iMdl = 1:2 % loop across the 2 models (no-wake, 4 wakes)
                      % lower boundary for the paramters to be estimated
                      lb = [-360 0 0 0 0 0 min(yG(i1:i2,iR)) min(yG(i1:i2,iR)) min(yG(i1:i2,iR)) min(yG(i1:i2,iR)) 0 0 0 0];
                      % upper boundary for the paramters to be estimated
                      ub = [360 +Inf u0 u0 u0 u0 +Inf +Inf +Inf +Inf +Inf +Inf +Inf +Inf];
                      % Regression
                      [bHat{itent,iMdl},ss(itent,iMdl),res(:,iMdl),exitflag{itent,iMdl},output{itent,iMdl},lambda{itent,iMdl},J{itent,iMdl}] = lsqcurvefit(mdl,b0{iMdl},yG(i1:i2,iR),losvGate.*weightsGate,lb,ub);

                      if iMdl == 2 % 4 wakes model
                         bHat{itent,iMdl}(7:10) = sort(bHat{itent,iMdl}(7:10)); % sort peaks locations
                      end
                      % losv determination using the 2 models, with the
                      % estimates of the parameter just found (no weights
                      % issues..)
                      losvHat(:,iMdl) = mdlnowgts(bHat{itent,iMdl},yG(i1:i2,iR)); 
                      % Calculate residuals
                      res(:,iMdl)=losvHat(:,iMdl)-losvGate(:,1);
                      % Mean squared error
                      mse(itent,iMdl) = (sum(res(:,iMdl).^2))/(nT-nP(iMdl));
                      % Norm of the residuals
                      ss(itent,iMdl) = sum(res(:,iMdl).^2);
                      % Covariance matrix
                      %cov{itent,iMdl} = inv(J{itent,iMdl}'*J{itent,iMdl}).*mse(itent,iMdl);
                      % Nonlinear regression parameter confidence intervals
                      ci{itent,iMdl} = nlparci(bHat{itent,iMdl},res(:,iMdl),'jacobian',J{itent,iMdl});
                      % Linear correlation coefficient (Pearson)
                      rho(itent,iMdl) = corr(losvGate,losvHat(:,iMdl));
                  end
                    
                  % Degrees of freedom
                  df = nT - nP;
                  % F-test statistic
                  F12(itent) = ((ss(itent,1) - ss(itent,2))/(df(1) - df(2)))/(ss(itent,2)/df(2));
                  % p-value
                  p12(itent) = 1 - fcdf(F12(itent),df(1)-df(2),df(2));

                  % Which model to be chosen for each fit try? 
                  % Determine whether there are wakes with F-test
                  if p12(itent) < 0.05
                     % wake type
                     wakeTemp(itent) = 2; % 4 wakes model
                     % parameter estimate record
                     bHatRecTemp(itent,:) = bHat{itent,wakeTemp(itent)};
                     % parameter confidence interval record
                     ciRecTemp(itent,:) = reshape(ci{itent,wakeTemp(itent)}',1,nP(wakeTemp(itent))*2);
                  else
                     % wake type
                     wakeTemp(itent) = 1; % no-wake model
                     % parameter estimate record
                     bHatRecTemp(itent,:) = [bHat{itent,wakeTemp(itent)} NaN(1,12)];
                     % parameter confidence interval record
                     ciRecTemp(itent,:) = [reshape(ci{itent,wakeTemp(itent)}',1,nP(wakeTemp(itent))*2) NaN(1,24)];
                  end
              end % end 4 fitting tries
              
              % Choose between 4 fitting tries basing on R^2 coefficient
              if rho(1,wakeTemp(1))>=rho(2,wakeTemp(2)) && rho(1,wakeTemp(1))>=rho(3,wakeTemp(3)) && rho(1,wakeTemp(1))>=rho(4,wakeTemp(4))
                  seltent=1;
                  bHatRec = bHatRecTemp(1,:);
                  ciRec = ciRecTemp(1,:);
                  wake = wakeTemp(1);
              elseif rho(2,wakeTemp(2))>=rho(1,wakeTemp(1)) && rho(2,wakeTemp(2))>=rho(3,wakeTemp(3)) && rho(2,wakeTemp(2))>=rho(4,wakeTemp(4))
                  seltent=2;
                  bHatRec = bHatRecTemp(2,:);
                  ciRec = ciRecTemp(2,:);
                  wake = wakeTemp(2);
              elseif rho(3,wakeTemp(3))>=rho(1,wakeTemp(1)) && rho(3,wakeTemp(3))>=rho(2,wakeTemp(2)) && rho(3,wakeTemp(3))>=rho(4,wakeTemp(4))
                  seltent=3;
                  bHatRec = bHatRecTemp(3,:);
                  ciRec = ciRecTemp(3,:);
                  wake = wakeTemp(3);
              else
                  seltent=4;
                  bHatRec = bHatRecTemp(4,:);
                  ciRec = ciRecTemp(4,:);
                  wake = wakeTemp(4);
              end
              
              % Convert phi units from rad to deg
              bHatRec(1) = bHatRec(1)*(180/pi);
              ciRec(1:2) = ciRec(1:2)*(180/pi);
                 
              
              % ____________________________________________________
              % QUALITY CHECKS
              
              % Exclude wakes with peak locations outside the actual
              % transverse coordinate range of the scan
              if wake==2
                 for ind=7:10
                     if bHatRec(ind)>max(yorig)-0.2
                        bHatRec(ind)=NaN;
                        bHatRec(ind+4)=NaN;
                        bHatRec(ind-4)=NaN;
                     end
                 end
                 ind = 7; % 1st wake outside the range?
                     if bHatRec(ind)<min(yorig)+0.2
                        bHatRec(ind)=bHatRec(ind+1);
                        bHatRec(ind+1)=bHatRec(ind+2);
                        bHatRec(ind+2)=bHatRec(ind+3);
                        bHatRec(ind+3)=NaN;
                        bHatRec(ind+4)=bHatRec(ind+5);
                        bHatRec(ind+5)=bHatRec(ind+6);
                        bHatRec(ind+6)=bHatRec(ind+7);
                        bHatRec(ind+7)=NaN;
                        bHatRec(ind-4)=bHatRec(ind-3);
                        bHatRec(ind-3)=bHatRec(ind-2);
                        bHatRec(ind-2)=bHatRec(ind-1);
                        bHatRec(ind-1)=NaN;
                     end
              end
        
              % Exclude wakes with very low width, or with very high width,
              % or with very low velocity deficit, and re-order 
              % the remaining fitted wakes to match the real wakes.
              % NOTE: the re-order procedure here works for southwesterly
              % wind; for other directions it needs to be adapted!
              if wake==2
                 ind=14; % 4th detected wake (on the right edge of the row)
                 A = [bHatRec(ind-9), bHatRec(ind-10), bHatRec(ind-11)];
                 %if bHatRec(ind)<0.5/4 || bHatRec(ind-8)<iR*0.03 || 4*bHatRec(ind)>(max(yorig)-min(yorig))/3
                 if bHatRec(ind)<0.1 || bHatRec(ind-8)<0.5*nanmean(A) || 4*bHatRec(ind)>(max(yorig)-min(yorig))/4
                     bHatRec(ind)=NaN;
                     bHatRec(ind-4)=NaN;
                     bHatRec(ind-8)=NaN; 
                 end
                 ind=13; % 3rd detected wake
                 %if (bHatRec(ind)<0.5/4 || bHatRec(ind-8)<iR*0.04 || 4*bHatRec(ind)>(max(yorig)-min(yorig))/3) && iR<44-1
                 A = [bHatRec(ind-8), bHatRec(ind-9), bHatRec(ind-10)];
                 if (bHatRec(ind)<0.1 || bHatRec(ind-8)<0.5*nanmean(A) || 4*bHatRec(ind)>(max(yorig)-min(yorig))/4) && iR<44-1
                     bHatRec(ind)=bHatRec(ind+1);
                     bHatRec(ind-4)=bHatRec(ind-3);
                     bHatRec(ind-8)=bHatRec(ind-7);
                     bHatRec(ind+1)=NaN;
                     bHatRec(ind-3)=NaN;
                     bHatRec(ind-7)=NaN;
                 end
                 ind=12; % 2nd detected wake
                 A = [bHatRec(ind-7), bHatRec(ind-8), bHatRec(ind-9)];
                 %if (bHatRec(ind)<0.5/4 || bHatRec(ind-8)<iR*0.05 || 4*bHatRec(ind)>(max(yorig)-min(yorig))/3) && iR<44-2
                 if (bHatRec(ind)<0.1 || bHatRec(ind-8)<0.5*nanmean(A) || 4*bHatRec(ind)>(max(yorig)-min(yorig))/4) && iR<44-2
                     bHatRec(ind)=bHatRec(ind+1);
                     bHatRec(ind-4)=bHatRec(ind-3);
                     bHatRec(ind-8)=bHatRec(ind-7);
                     bHatRec(ind+1)=bHatRec(ind+2);
                     bHatRec(ind-3)=bHatRec(ind-2);
                     bHatRec(ind-7)=bHatRec(ind-6);
                     bHatRec(ind+2)=NaN;
                     bHatRec(ind-2)=NaN;
                     bHatRec(ind-6)=NaN;
                 end
                 ind=11; % 1st detected wake (on the left edge of the row)
                 A = [bHatRec(ind-6), bHatRec(ind-7), bHatRec(ind-8)];
                 %if (bHatRec(ind)<0.5/4 || bHatRec(ind-8)<iR*0.046 || 4*bHatRec(ind)>(max(yorig)-min(yorig))/3) && iR<44-2
                 if (bHatRec(ind)<0.1 || bHatRec(ind-8)<0.5*nanmean(A) || 4*bHatRec(ind)>(max(yorig)-min(yorig))/4) && iR<44-2
                     bHatRec(ind)=bHatRec(ind+1);
                     bHatRec(ind-4)=bHatRec(ind-3);
                     bHatRec(ind-8)=bHatRec(ind-7);
                     bHatRec(ind+1)=bHatRec(ind+2);
                     bHatRec(ind-3)=bHatRec(ind-2);
                     bHatRec(ind-7)=bHatRec(ind-6);
                     bHatRec(ind+2)=bHatRec(ind+3);
                     bHatRec(ind-2)=bHatRec(ind-1);
                     bHatRec(ind-6)=bHatRec(ind-5);
                     bHatRec(ind+3)=NaN;
                     bHatRec(ind-1)=NaN;
                     bHatRec(ind-5)=NaN;
                 end
            
              end
                
              % When not all the 4 wakes are present, the algorithm may
              % detect fake double-wakes which are actually a single
              % physical wake.
              % Let's check and try to solve this issue here.
              % NOTE: the re-order procedure here works for southwesterly
              % wind; for other directions it needs to be adapted!
              if wake==2 && iR<max(rInd)-1 && iR>min(rInd)+1 && iTable>1
                 % look at 1st and 2nd detected wakes
                 if (bHatRec(8)-bHatRec(7))<(2*bHatRec(11)+2*bHatRec(12))
                    % A: check if it is a double wake 
                    if (bHatRec(8)-bHatRec(7))<max(2*bHatRec(11),2*bHatRec(12)) && (bHatRec(8)-bHatRec(7))<1.5
                        bHatRec(3)=(bHatRec(4)+bHatRec(3))/2;
                        bHatRec(7)=(bHatRec(8)+bHatRec(7))/2;
                        bHatRec(11)=(abs(bHatRec(8)+2*bHatRec(12)-bHatRec(7)+2*bHatRec(11)))/4;
                        bHatRec(4)=bHatRec(5);
                        bHatRec(5)=bHatRec(6);
                        bHatRec(6)=NaN;
                        bHatRec(8)=bHatRec(9);
                        bHatRec(9)=bHatRec(10);
                        bHatRec(10)=NaN;
                        bHatRec(12)=bHatRec(13);
                        bHatRec(13)=bHatRec(14);
                        bHatRec(14)=NaN;
                    % B: if it is only a "width issue", check if the fit
                    %    at the previous range gate was great to get some
                    %    info from it
                    elseif tablePPI_1D(iTable,48,iR - min(rInd))<0.3
                        if abs(bHatRec(7)-tablePPI_1D(iTable,11,iR - min(rInd)))>abs(bHatRec(8)-tablePPI_1D(iTable,11,iR - min(rInd)))
                            bHatRec(3)=bHatRec(4);
                            bHatRec(7)=bHatRec(8);
                            bHatRec(11)=bHatRec(12);
                            bHatRec(4)=NaN;
                            bHatRec(5)=NaN;
                            bHatRec(6)=NaN;
                            bHatRec(8)=NaN;
                            bHatRec(9)=NaN;
                            bHatRec(10)=NaN;
                            bHatRec(12)=NaN;
                            bHatRec(13)=NaN;
                            bHatRec(14)=NaN;
                        end
                    % C: if it is only a "width issue", check if the fit
                    %    at 2 previous range gates was great to get some
                    %    info from it
                    elseif tablePPI_1D(iTable,48,iR - min(rInd)-1)<0.3 % if fit at 2-previous range gate was great
                        if abs(bHatRec(7)-tablePPI_1D(iTable,11,iR - min(rInd)-1))>abs(bHatRec(8)-tablePPI_1D(iTable,11,iR - min(rInd)-1))
                            bHatRec(3)=bHatRec(4);
                            bHatRec(7)=bHatRec(8);
                            bHatRec(11)=bHatRec(12);
                            bHatRec(4)=NaN;
                            bHatRec(5)=NaN;
                            bHatRec(6)=NaN;
                            bHatRec(8)=NaN;
                            bHatRec(9)=NaN;
                            bHatRec(10)=NaN;
                            bHatRec(12)=NaN;
                            bHatRec(13)=NaN;
                            bHatRec(14)=NaN;
                        end
                    % D: otherwise, just use NaNs, and don't replace anything
                    else 
                        bHatRec(4)=NaN;
                        bHatRec(5)=NaN;
                        bHatRec(6)=NaN;
                        bHatRec(8)=NaN;
                        bHatRec(9)=NaN;
                        bHatRec(10)=NaN;
                        bHatRec(12)=NaN;
                        bHatRec(13)=NaN;
                        bHatRec(14)=NaN;
                    end
                 end
                 % look at 2nd and 3rd detected wakes
                 if (bHatRec(9)-bHatRec(8))<2*bHatRec(12)+2*bHatRec(13) 
                    if (bHatRec(9)-bHatRec(8))<max(2*bHatRec(12),2*bHatRec(13)) && (bHatRec(9)-bHatRec(8))<1.5
                        bHatRec(4)=(bHatRec(4)+bHatRec(5))/2;
                        bHatRec(8)=(bHatRec(8)+bHatRec(9))/2;
                        bHatRec(12)=(abs(bHatRec(9)+2*bHatRec(13)-bHatRec(8)+2*bHatRec(12)))/4;
                        bHatRec(5)=bHatRec(6);
                        bHatRec(6)=NaN;
                        bHatRec(9)=bHatRec(10);
                        bHatRec(10)=NaN;
                        bHatRec(13)=bHatRec(14);
                        bHatRec(14)=NaN;
                    elseif tablePPI_1D(iTable,48,iR - min(rInd))<0.3 % if fit at previous range gate was great
                        if abs(bHatRec(8)-tablePPI_1D(iTable,12,iR - min(rInd)))>abs(bHatRec(9)-tablePPI_1D(iTable,12,iR - min(rInd)))
                            bHatRec(4)=bHatRec(5);
                            bHatRec(8)=bHatRec(9);
                            bHatRec(12)=bHatRec(13);
                            bHatRec(5)=NaN;
                            bHatRec(6)=NaN;
                            bHatRec(9)=NaN;
                            bHatRec(10)=NaN;
                            bHatRec(13)=NaN;
                            bHatRec(14)=NaN;
                        end
                    elseif tablePPI_1D(iTable,48,iR - min(rInd)-1)<0.3 % if fit at 2-previous range gate was great
                        if abs(bHatRec(8)-tablePPI_1D(iTable,12,iR - min(rInd)-1))>abs(bHatRec(9)-tablePPI_1D(iTable,12,iR - min(rInd)-1))
                            bHatRec(4)=bHatRec(5);
                            bHatRec(8)=bHatRec(9);
                            bHatRec(12)=bHatRec(13);
                            bHatRec(5)=NaN;
                            bHatRec(6)=NaN;
                            bHatRec(9)=NaN;
                            bHatRec(10)=NaN;
                            bHatRec(13)=NaN;
                            bHatRec(14)=NaN;
                        end
                    else % otherwise, just use NaNs, and don't replace anything
                        bHatRec(5)=NaN;
                        bHatRec(6)=NaN;
                        bHatRec(9)=NaN;
                        bHatRec(10)=NaN;
                        bHatRec(13)=NaN;
                        bHatRec(14)=NaN;
                    end    
                 end
                 % look at 3rd and 4th detected wakes
                 if (bHatRec(10)-bHatRec(9))<2*bHatRec(13)+2*bHatRec(14)
                    if (bHatRec(10)-bHatRec(9))<max(2*bHatRec(13),2*bHatRec(14)) && (bHatRec(10)-bHatRec(9))<1.5
                        bHatRec(5)=(bHatRec(6)+bHatRec(5))/2;
                        bHatRec(9)=(bHatRec(10)+bHatRec(9))/2;
                        bHatRec(13)=(abs(bHatRec(10)+2*bHatRec(14)-bHatRec(9)+2*bHatRec(13)))/4;
                    elseif tablePPI_1D(iTable,48,iR - min(rInd))<0.3 % if fit at 2-previous range gate was great
                        if  abs(bHatRec(9)-tablePPI_1D(iTable,13,iR - min(rInd)))>abs(bHatRec(10)-tablePPI_1D(iTable,13,iR - min(rInd)))
                            bHatRec(5)=bHatRec(6);
                            bHatRec(9)=bHatRec(10);
                            bHatRec(13)=bHatRec(14);
                        end
                    elseif tablePPI_1D(iTable,48,iR - min(rInd)-1)<0.3 % if fit at 2-previous range gate was great
                        if  abs(bHatRec(9)-tablePPI_1D(iTable,13,iR - min(rInd)-1))>abs(bHatRec(10)-tablePPI_1D(iTable,13,iR - min(rInd)-1))
                            bHatRec(5)=bHatRec(6);
                            bHatRec(9)=bHatRec(10);
                            bHatRec(13)=bHatRec(14);
                        end
                    end % otherwise, just use NaNs, and don't replace anything
                        bHatRec(6)=NaN;
                        bHatRec(10)=NaN;
                        bHatRec(14)=NaN;
                 end
              end

              % ______________________________________
              % OUTPUT: table (i,j,k)
              % i -> each single PPI scan
              % j -> different quantities stored
              % k -> each single range gate
              tablePPI_1D(iTable,:,iR - min(rInd) + 1) = [decTime(i1) meanDir meanSpd wake bHatRec ciRec rho(seltent,wake) mse(seltent,wake) min(yG(i1:i2,iR)) max(yG(i1:i2,iR)) el(i1,1)*(180/pi) bHatFree spaceTI LObuk];
              % ______________________________________
              
              % ______________________________________
              % PLOTS: measured losv and fit
              % yPlot = linspace(min(yG(i1:i2,iR)),max(yG(i1:i2,iR)),1000);
              % losvPlot = mdlnowgts(bHat{seltent,wake},yPlot);
              % figure('vis','off');
              % set(gcf,'Color','White','PaperPosition',[1 1 4 4]);
              % errorbar(yG(i1:i2,iR),losvGate,rdispGate,'xk')
              % hold on
              % plot(yPlot,losvPlot,'r',lws,lw);
              % ylim([0,12]);
              % xlabel('y'' [D]','FontSize',14);
              % ylabel('u_LOS [m/s]','FontSize',14);
              % str = sprintf('PPI scan %s UTC, Elevation angle: %0.2f%c',datestr(decTime(i1)),el(i1,1)*(180/pi),char(176));
              % str2 = sprintf('Range gate %0.0f m, Wind Dir %0.1f, MSE %0.2f, Rho %0.2f',range(iR)*D, meanDir, mse(seltent,wake),rho(seltent,wake));
              % title({str;str2},'FontSize',6);
              % if wake==2
              %     str3 = sprintf('a = (%0.1f %0.1f %0.1f %0.1f)',bHatRec(1,3),bHatRec(1,4),bHatRec(1,5),bHatRec(1,6));
              %     text(-1.5,0.8,str3,'FontSize',8);
              %     str4 = sprintf('mu = (%0.1f %0.1f %0.1f %0.1f)',bHatRec(1,7),bHatRec(1,8),bHatRec(1,9),bHatRec(1,10));
              %     text(-1.5,1.3,str4,'FontSize',8);
              %     str5 = sprintf('s = (%0.1f %0.1f %0.1f %0.1f)',4*bHatRec(1,11),4*bHatRec(1,12),4*bHatRec(1,13),4*bHatRec(1,14));
              %     text(-1.5,0.3,str5,'FontSize',8);
              % end      
              % saveas(gcf,sprintf('%s_%0.1f_%0.0f.jpg',datestr(decTime(i1)),el(i1,1)*(180/pi),range(iR)*D));
              % ______________________________________
          end
          % Reset index for next iteration of the loop
          iTable = iTable + 1;
      end
      % Reset temporal index for next iteration of the loop
      i1 = i2 + 1;
end

% Finally, save the output table
save tablePPI_1D_2017 tablePPI_1D