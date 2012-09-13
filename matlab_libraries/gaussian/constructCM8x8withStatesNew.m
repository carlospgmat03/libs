% $Id$
%
% constructCM8x8RawData.m
%
% Description:
%
% This file reconstructs the 8x8 Covariance matrix from the raw data
% collected during a measurement.
%
%

% clear all; close all;

addpath '~/investigacion/matlab_libraries/directories';
% path
all_cfp_dirs;

addpath '~/investigacion/bound/JamesToolBox/CVQinfo'
addpath '~/investigacion/bound/JamesToolBox/CVQinfo/MathTools'

%% Read in the data

Ndata = 4e6; % That is 4'000,000
% Ndata = 5e5; % That is   500,000
% [dset1, dset2, dset3, dset4, dset5, dset6, dsetsub] = read4modeInCM(Ndata);


%% Determine number of samples

% Initialize the matrices
gamma = zeros(8,8);
Nmodes = 4;


Ndata = size(dset1(:,1),1);
%SizeInterval = 50000;
% SizeInterval = Ndata;
Nsamples = 5;
SizeInterval = Ndata/Nsamples;
results = zeros(Nsamples,3);

last = SizeInterval;
start = 1;

for idx = 1:Nsamples
    % Calculate main diagonal
    Vx1 = var(dset1(start:last,1)); Vp1 = var(dset2(start:last,1));
    Vx2 = var(dset1(start:last,2)); Vp2 = var(dset2(start:last,2));
    Vx3 = var(dset1(start:last,3)); Vp3 = var(dset2(start:last,3));
    Vx4 = var(dset1(start:last,4)); Vp4 = var(dset2(start:last,4));
    
    gamma(1,1) = Vx1; gamma(3,3) = Vx2; gamma(5,5) = Vx3; gamma(7,7) = Vx4; 
    gamma(2,2) = Vp1; gamma(4,4) = Vp2; gamma(6,6) = Vp3; gamma(8,8) = Vp4; 
    
    % Calculate the covariances
    Covx1x2 = cov(dset1(start:last,1), dset1(start:last,2)); gamma(1,3) = Covx1x2(1,2); % Covx1x2
    Covx1x3 = cov(dset1(start:last,1), dset1(start:last, 3)); gamma(1,5) = Covx1x3(1,2); % Covx1x3
    Covx1x4 = cov(dset1(start:last,1), dset1(start:last, 4)); gamma(1,7) = Covx1x4(1,2); % Covx1x4
    Covx2x3 = cov(dset1(start:last,2), dset1(start:last,3)); gamma(3,5) = Covx2x3(1,2); % Covx2x3
    Covx2x4 = cov(dset1(start:last,2), dset1(start:last,4)); gamma(3,7) = Covx2x4(1,2); % Covx2x4
    Covx3x4 = cov(dset1(start:last,3), dset1(start:last,4)); gamma(5,7) = Covx3x4(1,2); % Covx3x4

    Covp1p2 = cov(dset2(start:last,1), dset2(start:last,2)); gamma(2,4) = Covp1p2(1,2); % Covp1p2
    Covp1p3 = cov(dset2(start:last,1), dset2(start:last,3)); gamma(2,6) = Covp1p3(1,2); % Covp1p3
    Covp1p4 = cov(dset2(start:last,1), dset2(start:last,4)); gamma(2,8) = Covp1p4(1,2); % Covp1p4
    Covp2p3 = cov(dset2(start:last,2), dset2(start:last,3)); gamma(4,6) = Covp2p3(1,2); % Covp2p3
    Covp2p4 = cov(dset2(start:last,2), dset2(start:last,4)); gamma(4,8) = Covp2p4(1,2); % Covp2p4
    Covp3p4 = cov(dset2(start:last,3), dset2(start:last,4)); gamma(6,8) = Covp3p4(1,2); % Covp3p4
    
    Vl1 = var(dsetsub(start:last,1)); Vl2 = var(dsetsub(start:last,2));
    Vl3 = var(dsetsub(start:last,3)); Vl4 = var(dsetsub(start:last,4));

    % Now calculate the covariance
    Covx1p1 = Vl1 - 0.5*Vx1 - 0.5*Vp1; gamma(1,2) = Covx1p1;
    Covx2p2 = Vl2 - 0.5*Vx2 - 0.5*Vp2; gamma(3,4) = Covx2p2;
    Covx3p3 = Vl3 - 0.5*Vx3 - 0.5*Vp3; gamma(5,6) = Covx3p3;
    Covx4p4 = Vl4 - 0.5*Vx4 - 0.5*Vp4; gamma(7,8) = Covx4p4;
    
    % Calculate the cross variances
    
    Covx1p2 = cov(dset3(start:last, 1), dset3(start:last, 2)); gamma(1,4) = Covx1p2(1,2);
    Covx1p3 = cov(dset4(start:last, 1), dset4(start:last, 3)); gamma(1,6) = Covx1p3(1,2);
    Covx1p4 = cov(dset5(start:last, 1), dset5(start:last, 4)); gamma(1,8) = Covx1p4(1,2);

    Covp1x2 = cov(dset6(start:last, 1), dset6(start:last, 2)); gamma(2,3) = Covp1x2(1,2);
    Covp3x2 = cov(dset4(start:last, 3), dset4(start:last, 2)); gamma(3,6) = Covp3x2(1,2);
    Covp4x2 = cov(dset5(start:last, 4), dset5(start:last, 2)); gamma(3,8) = Covp4x2(1,2);

    Covx3p1 = cov(dset6(start:last, 3), dset6(start:last, 1)); gamma(2,5) = Covx3p1(1,2);
    Covx3p2 = cov(dset3(start:last, 3), dset3(start:last, 2)); gamma(4,5) = Covx3p2(1,2);
    Covx3p4 = cov(dset5(start:last, 3), dset5(start:last, 4)); gamma(5,8) = Covx3p4(1,2);

    Covx4p1 = cov(dset6(start:last, 4), dset6(start:last, 1)); gamma(2,7) = Covx4p1(1,2);
    Covx4p2 = cov(dset3(start:last, 4), dset3(start:last, 2)); gamma(4,7) = Covx4p2(1,2);
    Covx4p3 = cov(dset4(start:last, 4), dset4(start:last, 3)); gamma(6,7) = Covx4p3(1,2);
    
    gamma = gamma + triu(gamma,1)';
    [J] = symplecticForm(Nmodes);
    simon = min(real(eig(gamma + i*J)));

    if (simon<0)
	gamma =  physicalize(gamma, '2norm');
% 	gamma =  physicalize(gamma, '1norm');
    end
    simon = min(real(eig(gamma + i*J)));
    ppt = cmPeres(gamma, 'double', [2 2]);
    En = how_entangled(gamma, [2 2]);
    results(idx,1) = ppt; results(idx,2) = En; results(idx,3) = simon;
    start = last; last = last+SizeInterval; gamma = zeros(8,8);
    
end

results

Ndata = size(dset1(:,1),1);
%SizeInterval = 50000;
% SizeInterval = Ndata;
Nsamples = 1;
SizeInterval = Ndata/Nsamples;
results = zeros(Nsamples,3);

last = SizeInterval;
start = 1;

for idx = 1:Nsamples
    % Calculate main diagonal
    Vx1 = var(dset1(start:last,1)); Vp1 = var(dset2(start:last,1));
    Vx2 = var(dset1(start:last,2)); Vp2 = var(dset2(start:last,2));
    Vx3 = var(dset1(start:last,3)); Vp3 = var(dset2(start:last,3));
    Vx4 = var(dset1(start:last,4)); Vp4 = var(dset2(start:last,4));
    
    gamma(1,1) = Vx1; gamma(3,3) = Vx2; gamma(5,5) = Vx3; gamma(7,7) = Vx4; 
    gamma(2,2) = Vp1; gamma(4,4) = Vp2; gamma(6,6) = Vp3; gamma(8,8) = Vp4; 
    
    % Calculate the covariances
    Covx1x2 = cov(dset1(start:last,1), dset1(start:last,2)); gamma(1,3) = Covx1x2(1,2); % Covx1x2
    Covx1x3 = cov(dset1(start:last,1), dset1(start:last, 3)); gamma(1,5) = Covx1x3(1,2); % Covx1x3
    Covx1x4 = cov(dset1(start:last,1), dset1(start:last, 4)); gamma(1,7) = Covx1x4(1,2); % Covx1x4
    Covx2x3 = cov(dset1(start:last,2), dset1(start:last,3)); gamma(3,5) = Covx2x3(1,2); % Covx2x3
    Covx2x4 = cov(dset1(start:last,2), dset1(start:last,4)); gamma(3,7) = Covx2x4(1,2); % Covx2x4
    Covx3x4 = cov(dset1(start:last,3), dset1(start:last,4)); gamma(5,7) = Covx3x4(1,2); % Covx3x4

    Covp1p2 = cov(dset2(start:last,1), dset2(start:last,2)); gamma(2,4) = Covp1p2(1,2); % Covp1p2
    Covp1p3 = cov(dset2(start:last,1), dset2(start:last,3)); gamma(2,6) = Covp1p3(1,2); % Covp1p3
    Covp1p4 = cov(dset2(start:last,1), dset2(start:last,4)); gamma(2,8) = Covp1p4(1,2); % Covp1p4
    Covp2p3 = cov(dset2(start:last,2), dset2(start:last,3)); gamma(4,6) = Covp2p3(1,2); % Covp2p3
    Covp2p4 = cov(dset2(start:last,2), dset2(start:last,4)); gamma(4,8) = Covp2p4(1,2); % Covp2p4
    Covp3p4 = cov(dset2(start:last,3), dset2(start:last,4)); gamma(6,8) = Covp3p4(1,2); % Covp3p4
    
    Vl1 = var(dsetsub(start:last,1)); Vl2 = var(dsetsub(start:last,2));
    Vl3 = var(dsetsub(start:last,3)); Vl4 = var(dsetsub(start:last,4));

    % Now calculate the covariance
    Covx1p1 = Vl1 - 0.5*Vx1 - 0.5*Vp1; gamma(1,2) = Covx1p1;
    Covx2p2 = Vl2 - 0.5*Vx2 - 0.5*Vp2; gamma(3,4) = Covx2p2;
    Covx3p3 = Vl3 - 0.5*Vx3 - 0.5*Vp3; gamma(5,6) = Covx3p3;
    Covx4p4 = Vl4 - 0.5*Vx4 - 0.5*Vp4; gamma(7,8) = Covx4p4;
    
    % Calculate the cross variances
    
    Covx1p2 = cov(dset3(start:last, 1), dset3(start:last, 2)); gamma(1,4) = Covx1p2(1,2);
    Covx1p3 = cov(dset4(start:last, 1), dset4(start:last, 3)); gamma(1,6) = Covx1p3(1,2);
    Covx1p4 = cov(dset5(start:last, 1), dset5(start:last, 4)); gamma(1,8) = Covx1p4(1,2);

    Covp1x2 = cov(dset6(start:last, 1), dset6(start:last, 2)); gamma(2,3) = Covp1x2(1,2);
    Covp3x2 = cov(dset4(start:last, 3), dset4(start:last, 2)); gamma(3,6) = Covp3x2(1,2);
    Covp4x2 = cov(dset5(start:last, 4), dset5(start:last, 2)); gamma(3,8) = Covp4x2(1,2);

    Covx3p1 = cov(dset6(start:last, 3), dset6(start:last, 1)); gamma(2,5) = Covx3p1(1,2);
    Covx3p2 = cov(dset3(start:last, 3), dset3(start:last, 2)); gamma(4,5) = Covx3p2(1,2);
    Covx3p4 = cov(dset5(start:last, 3), dset5(start:last, 4)); gamma(5,8) = Covx3p4(1,2);

    Covx4p1 = cov(dset6(start:last, 4), dset6(start:last, 1)); gamma(2,7) = Covx4p1(1,2);
    Covx4p2 = cov(dset3(start:last, 4), dset3(start:last, 2)); gamma(4,7) = Covx4p2(1,2);
    Covx4p3 = cov(dset4(start:last, 4), dset4(start:last, 3)); gamma(6,7) = Covx4p3(1,2);
    
    gamma = gamma + triu(gamma,1)';
    [J] = symplecticForm(Nmodes);
    simon = min(real(eig(gamma + i*J)));

    if (simon<0)
% 	gamma =  physicalize(gamma, '1norm');
	gamma =  physicalize(gamma, '2norm');
    end
    simon = min(real(eig(gamma + i*J)));
    ppt = cmPeres(gamma, 'double', [2 2]);
    En = how_entangled(gamma, [2 2]);
    results(idx,1) = ppt; results(idx,2) = En; results(idx,3) = simon;
    start = last; last = last+SizeInterval; gamma = zeros(8,8);
    
end

%% Make some plots

% % Plot the results
% Nbins = 30;
% figure(1)
% hist(results(:,1), Nbins);
% title('PPT')
% 
% figure(2)
% hist(results(:,2), Nbins);
% title('Separability')
% 
% figure(3)
% hist(results(:,3), Nbins);
% title('Bona Fide')
% 

results

