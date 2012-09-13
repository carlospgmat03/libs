% $Id$
%
% read4modeInCM.m
%
% Description:
%
% This Matlab-file reads in the data collected from the four homodyne
% detectors for the bound entanglement experiment.  There are six data
% files in total corresponding to the following measurement order:
%
%
% dataSetShot -- shota, shotc, shotb, shotb
%
% dataSet1 -- xa, xc, xb, xd
%
% dataSet2 -- pa, pc, pb, pd
%
% dataSet3 -- xa, pc, xb, xd
%
% dataSet4 -- xa, xc, pb, pd
%
% dataSet5 -- xa, xc, xb, pd
%
% dataSet6 -- pa, xc, xb, xd
%
% This order has been determined by the order in which the homodye
% detectors are read-in with the data acquisition system.
%
%
% Inputs:
%
% pathname -- A character string specifying the date and trial of the raw
%             data to be used.
%
% Outputs:
%
% dset1-sub -- A vector containing four columns of Ndata points
%	  
%
% Author: James DiGuglielmo
%
% Last Updated: May 19. 2010

function [xxxx, pppp, xpxx, xxpx, xxxp, pxxx, acbdsub] = read4modeInCM(Ndata)

% Open the files
fidshot = fopen('dataSetshot.dat','r');
fidset1 = fopen('dataSet1.dat','r');
fidset2 = fopen('dataSet2.dat','r');
fidset3 = fopen('dataSet3.dat','r');
fidset4 = fopen('dataSet4.dat','r');
fidset5 = fopen('dataSet5.dat','r');
fidset6 = fopen('dataSet6.dat','r');
fidsetsub = fopen('dataSetsub.dat','r');

% Read in the data
dshot = fread(fidshot, [Ndata, inf], 'float64');
dset1 = fread(fidset1, [Ndata, inf], 'float64');
dset2 = fread(fidset2, [Ndata, inf], 'float64');
dset3 = fread(fidset3, [Ndata, inf], 'float64');
dset4 = fread(fidset4, [Ndata, inf], 'float64');
dset5 = fread(fidset5, [Ndata, inf], 'float64');
dset6 = fread(fidset6, [Ndata, inf], 'float64');
dsetsub = fread(fidsetsub, [Ndata, inf], 'float64');

% Close the files
fclose(fidshot);
fclose(fidset1);    
fclose(fidset2);
fclose(fidset3);
fclose(fidset4);
fclose(fidset5);
fclose(fidset6);
fclose(fidsetsub);

% Subtract DC from the raw data
for cnt = 1:size(dshot, 2)
    dshot(:,cnt) = dshot(:,cnt) - mean(dshot(:,cnt));
    dset1(:,cnt) = dset1(:,cnt) - mean(dset1(:,cnt));
    dset2(:,cnt) = dset2(:,cnt) - mean(dset2(:,cnt));
    dset3(:,cnt) = dset3(:,cnt) - mean(dset3(:,cnt));
    dset4(:,cnt) = dset4(:,cnt) - mean(dset4(:,cnt));
    dset5(:,cnt) = dset5(:,cnt) - mean(dset5(:,cnt));
    dset6(:,cnt) = dset6(:,cnt) - mean(dset6(:,cnt));
    dsetsub(:,cnt) = dsetsub(:,cnt) - mean(dsetsub(:,cnt));
end

% Set normalization
shotstd = std(dshot);         % vacuum equals 1


% Initialize data matrix
xxxx = zeros(Ndata, 4);
pppp = zeros(Ndata, 4);
xpxx = zeros(Ndata, 4);
xxpx = zeros(Ndata, 4);
xxxp = zeros(Ndata, 4);
pxxx = zeros(Ndata, 4);
sub = zeros(Ndata, 4);

% Normalize to shot noise

% Measuring the amplitude quadratures
% size(xxxx)
% size(dset1)
xxxx(:,1) = dset1(:,1)/shotstd(1); % xa
xxxx(:,2) = dset1(:,2)/shotstd(2); % xc
xxxx(:,3) = dset1(:,3)/shotstd(3); % xb
xxxx(:,4) = dset1(:,4)/shotstd(4); % xd

% Measuring phase quadratures
pppp(:,1) = dset2(:,1)/shotstd(1); % pa
pppp(:,2) = dset2(:,2)/shotstd(2); % pc
pppp(:,3) = dset2(:,3)/shotstd(3); % pb
pppp(:,4) = dset2(:,4)/shotstd(4); % pd

% Measuring: xpxxx
xpxx(:,1) = dset3(:,1)/shotstd(1); % xa
xpxx(:,2) = dset3(:,2)/shotstd(2); % pc
xpxx(:,3) = dset3(:,3)/shotstd(3); % xb
xpxx(:,4) = dset3(:,4)/shotstd(4); % xd

% Measuring: xxpx
xxpx(:,1) = dset4(:,1)/shotstd(1); % xa
xxpx(:,2) = dset4(:,2)/shotstd(2); % xc
xxpx(:,3) = dset4(:,3)/shotstd(3); % pb
xxpx(:,4) = dset4(:,4)/shotstd(4); % xd

% Measuring: xxxp
xxxp(:,1) = dset5(:,1)/shotstd(1); % xa
xxxp(:,2) = dset5(:,2)/shotstd(2); % xc
xxxp(:,3) = dset5(:,3)/shotstd(3); % xb
xxxp(:,4) = dset5(:,4)/shotstd(4); % pd


% Measuring: pxxx
pxxx(:,1) = dset6(:,1)/shotstd(1); % pa
pxxx(:,2) = dset6(:,2)/shotstd(2); % xc
pxxx(:,3) = dset6(:,3)/shotstd(3); % xb
pxxx(:,4) = dset6(:,4)/shotstd(4); % xd

% Measuring 45 lock
acbdsub(:,1) = dsetsub(:,1)/shotstd(1); % xa [45 deg]
acbdsub(:,2) = dsetsub(:,2)/shotstd(2); % xc [45 deg]
acbdsub(:,3) = dsetsub(:,3)/shotstd(3); % xb [45 deg]
acbdsub(:,4) = dsetsub(:,4)/shotstd(4); % xd [45 deg]





