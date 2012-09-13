% $Id: read4modeIn.m 2887 2010-05-04 15:38:24Z jadigu $
%
% read4modeIn.m
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
% In total, 500,000 data values will be collected.
%
% Inputs:
%
% pathname -- A character string specifying the date and trial of the raw
%             data to be used.
%
% Outputs:
%
% data -- Ndata x 2 x 4 x 7 Matrix where the dimensions mean
%
%	  i.) Ndata -- The total number of raw data points
%
%	  ii.) 2 --  Specifies individual quadrature, e.g. x1, x2, ...etc
%
%	  iii.) 4 -- Specifies the 2x2 system, e.g. x1x2, p1p2, x1p2, p2x1
%
%	  iv.) 7 -- Specifies the 4x4 system, e.g. gammaAB, gammaAC, ...etc 
%               including the 45 degree lock.
%
% Author: James DiGuglielmo
%
% Last Updated: Dec. 6 2009
%               April 1, 2010
%               April 29, 2010

function [data, dset1, dset2, dset3, dset4, dset5,dset6] = read4modeIn(Ndata)

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
data = zeros(Ndata, 2, 4, 7);

% Build data matrix

% GammaAB
data(:,1,1,1) = dset1(:,1)/shotstd(1); % xa
data(:,2,1,1) = dset1(:,3)/shotstd(3); % xb

data(:,1,2,1) = dset2(:,1)/shotstd(1); % pa
data(:,2,2,1) = dset2(:,3)/shotstd(3); % pb

data(:,1,3,1) = dset4(:,1)/shotstd(1); % xa
data(:,2,3,1) = dset4(:,3)/shotstd(3); % pb

data(:,1,4,1) = dset6(:,1)/shotstd(1); % pa
data(:,2,4,1) = dset6(:,3)/shotstd(3); % xb


% GammaAC
data(:,1,1,2) = dset1(:,1)/shotstd(1); % xa
data(:,2,1,2) = dset1(:,2)/shotstd(2); % xc

data(:,1,2,2) = dset2(:,1)/shotstd(1); % pa
data(:,2,2,2) = dset2(:,2)/shotstd(2); % pc

data(:,1,3,2) = dset3(:,1)/shotstd(1); % xa
data(:,2,3,2) = dset3(:,2)/shotstd(2); % pc

data(:,1,4,2) = dset6(:,1)/shotstd(1); % pa
data(:,2,4,2) = dset6(:,2)/shotstd(2); % xc


% GammaAD
data(:,1,1,3) = dset1(:,1)/shotstd(1); % xa
data(:,2,1,3) = dset1(:,4)/shotstd(4); % xd

data(:,1,2,3) = dset2(:,1)/shotstd(1); % pa
data(:,2,2,3) = dset2(:,4)/shotstd(4); % pd

data(:,1,3,3) = dset5(:,1)/shotstd(1); % xa
data(:,2,3,3) = dset5(:,4)/shotstd(4); % pd

data(:,1,4,3) = dset6(:,1)/shotstd(1); % pa
data(:,2,4,3) = dset6(:,4)/shotstd(4); % xd


% GammaCB
data(:,1,1,4) = dset1(:,2)/shotstd(2); % xc
data(:,2,1,4) = dset1(:,3)/shotstd(3); % xb

data(:,1,2,4) = dset2(:,2)/shotstd(2); % pc
data(:,2,2,4) = dset2(:,3)/shotstd(3); % pb

data(:,1,3,4) = dset4(:,2)/shotstd(2); % xc
data(:,2,3,4) = dset4(:,3)/shotstd(3); % pb

data(:,1,4,4) = dset3(:,2)/shotstd(2); % pc
data(:,2,4,4) = dset3(:,3)/shotstd(3); % xb

% GammaBD
data(:,1,1,5) = dset1(:,3)/shotstd(3); % xb
data(:,2,1,5) = dset1(:,4)/shotstd(4); % xd

data(:,1,2,5) = dset2(:,3)/shotstd(3); % pb
data(:,2,2,5) = dset2(:,4)/shotstd(4); % pd

data(:,1,3,5) = dset5(:,3)/shotstd(3); % xb
data(:,2,3,5) = dset5(:,4)/shotstd(4); % pd

data(:,1,4,5) = dset4(:,3)/shotstd(3); % pb
data(:,2,4,5) = dset4(:,4)/shotstd(4); % xd

    
% GammaCD
data(:,1,1,6) = dset1(:,2)/shotstd(2); % xc
data(:,2,1,6) = dset1(:,4)/shotstd(4); % xd

data(:,1,2,6) = dset2(:,2)/shotstd(2); % pc
data(:,2,2,6) = dset2(:,4)/shotstd(4); % pd

data(:,1,3,6) = dset5(:,2)/shotstd(2); % xc
data(:,2,3,6) = dset5(:,4)/shotstd(4); % pd

data(:,1,4,6) = dset3(:,2)/shotstd(2); % pc
data(:,2,4,6) = dset3(:,4)/shotstd(4); % xd

% 45 degree lock
data(:,1,1,7) = dsetsub(:,1)/shotstd(1); % xa [45 deg]
data(:,1,2,7) = dsetsub(:,2)/shotstd(2); % xc [45 deg]
data(:,1,3,7) = dsetsub(:,3)/shotstd(3); % xb [45 deg]
data(:,1,4,7) = dsetsub(:,4)/shotstd(4); % xd [45 deg]





