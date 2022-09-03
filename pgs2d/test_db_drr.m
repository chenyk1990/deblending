% Benchmark marine (2D) field data example for deblending
%
% This example was used in
% Chen, Y., S. Fomel, and R. Abma, 2022, Joint deblending and source time
% inversion, Geophysics, doi: 10.1190/geo2022-0149.1.

close all; clc;clear;

%% Please change the directory path
% requiring the DRR package
% https://github.com/chenyk1990/MATdrr
addpath(genpath('~/MATdrr'));
addpath(genpath('../subroutines'));

%% please download data from https://drive.google.com/file/d/1ge0Mn_SB4LUsVgOBvATh0iISwGQahKh4/view?usp=sharing
load ~/Downloads/PGS/blended_data.mat
%% in this dataset
%there are two sources data3d(:,:,1:60) and data3d(:,:,61:120)
%each source contains 120 shots
%there are 60 receivers
%

% test.m
% A script for demonstrating 2D deblending in different cases
% Written by Yangkang Chen
%% here nt=1001;can be 1500, 1700, etc;
nt=1500;
dd=blended_data(1:1500,:,:);
dd=yc_transp(dd,23);
figure;yc_imagesc(dd(:,:,100),90);

figure;yc_imagesc(dd(:,:,1));
figure;yc_imagesc(squeeze(dd(:,2,:)));

mode=7;
nr=256; % number of receivers
[nt,nx,nr]=size(dd);
[parm.ntpad,parm.nxpad,parm.nr]=size(dd);
parm.nt=nt;parm.nx=nx;
parm.verb=1;
rt0=[1:nx]'*nt;
rand('state',201415);
rt1=floor([1:nx]'*nt*0.5+(rand(nx,1)-0.5)*nt);
parm.rt=rt1;
[bb]=dblendsr2d(dd,mode,parm);
figure;yc_imagesc([dd(:,:,1),bb(:,:,1)],90);
figure;yc_imagesc([squeeze(dd(:,1,:)),squeeze(bb(:,1,:))],90);

figure;plot(1:nx,rt0*0.004,'r*');hold on;plot(1:nx,rt1*0.004,'bo');
title('Shot shedule');legend('Traditional','Faster');
xlabel('Shot');ylabel('Shooting time (s)');

%% global SSA
dt=0.002;
lambda=0.5;
alpha=0.8; %control the stability
niter=10;
snr=[];
parm.nr=1;
dbb=zeros(size(bb));
for ir=1:nr
    m0=bb(:,:,ir);
    m=m0;
    for iter=1:niter
        m1=m+lambda*(bb(:,:,ir)-dblendsr2d(m,mode,parm));
%         m=fxmssa(m1,0,100,dt,5,0);
        m1(1:500,:) = yc_mf(m1(1:500,:),3,1,2);
        m1(501:1000,:) = yc_mf(m1(501:1000,:),6,1,2);
        m1(1001:end,:) = yc_mf(m1(1001:end,:),9,1,2);
        m=drr3d_win(m1,0,80,0.004,2,4,0,100,64,1,0.5,0.5,0.5);
        m = m + (m1-m) * alpha;  % alpha=0.75/maxfold
        snr(iter)=yc_snr(dd(:,:,ir),m);
        fprintf('ir=%d/%d, current iterations %d/%d, SNR=%g\n',ir,nr,iter,niter,snr(iter));
    end
    dbb(:,:,ir)=m;
    fprintf('ir=%d/%d, SNR=%g\n',ir,nr,snr(iter));
end

figure;yc_imagesc([dd(:,:,1),bb(:,:,1),dbb(:,:,1)],90);
figure;yc_imagesc([dd(:,:,64),bb(:,:,64),dbb(:,:,64)],90);
figure;yc_imagesc([dd(:,:,128),bb(:,:,128),dbb(:,:,128)],90);
figure;yc_imagesc([dd(:,:,200),bb(:,:,200),dbb(:,:,200)],90);

figure;yc_imagesc([squeeze(dd(:,1,:)),squeeze(bb(:,1,:)),squeeze(dbb(:,1,:))],90);
figure;yc_imagesc([squeeze(dd(:,100,:)),squeeze(bb(:,100,:)),squeeze(dbb(:,100,:))],90);
figure;yc_imagesc([squeeze(dd(:,128,:)),squeeze(bb(:,128,:)),squeeze(dbb(:,128,:))],90);
figure;yc_imagesc([squeeze(dd(:,64,:)),squeeze(bb(:,64,:)),squeeze(dbb(:,64,:))],90);
figure;yc_imagesc([squeeze(dd(:,200,:)),squeeze(bb(:,200,:)),squeeze(dbb(:,200,:))],90);

% save('pgs_drr.mat','dbb','rt0','rt1','-v7.3')
