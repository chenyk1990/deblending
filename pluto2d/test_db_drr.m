% Benchmark marine (2D) field data example for deblending
%
% This example was used in
% Chen, Y., M. Zhou, and R. Abma, 2022, Two practical ways for improving the deblending performance in marine towed-streamer acquisition.

close all; clc;clear;

%% Please change the directory path
% requiring the DRR package
% https://github.com/chenyk1990/MATdrr
% addpath(genpath('~/MATdrr'));
addpath(genpath('../subroutines'));

%% please download data from https://drive.google.com/file/d/1ge0Mn_SB4LUsVgOBvATh0iISwGQahKh4/view?usp=sharing
load pluto2d.mat

nt=1024;
dd=data(1:nt,:,:);
% dd=yc_transp(dd,23);
figure;yc_imagesc(dd(:,:,50),99);

figure;yc_imagesc(dd(:,:,1));
figure;yc_imagesc(squeeze(dd(:,2,:)));

mode=7;
[nt,nx,nr]=size(dd);
[parm.ntpad,parm.nxpad,parm.nr]=size(dd);
parm.nt=nt;parm.nx=nx;
parm.verb=1;
rt0=[1:nx]'*nt;
rand('state',201415);
rt1=floor([1:nx]'*nt*0.75+(rand(nx,1)-0.5)*nt);
% maxfold=2;ny=1;ntpad=nt;nxpad=nx;nypad=ny;rt=rt1;
% rt1 = db_gstimes(maxfold,nt,nx,ny,ntpad,nxpad,nypad);
parm.rt=rt1;
[bb]=dblendsr2d(dd,mode,parm);
figure;yc_imagesc([dd(:,:,1),bb(:,:,1)],99);
figure;yc_imagesc([squeeze(dd(:,1,:)),squeeze(bb(:,1,:))],99);

figure;plot(1:nx,rt0*0.004,'r*');hold on;plot(1:nx,rt1*0.004,'bo');
title('Shot shedule');legend('Traditional','Faster');
xlabel('Shot');ylabel('Shooting time (s)');

%% global SSA
dt=0.008;
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
%         m1(1:500,:) = yc_mf(m1(1:500,:),3,1,2);
%         m1(501:1000,:) = yc_mf(m1(501:1000,:),6,1,2);
%         m1(1001:end,:) = yc_mf(m1(1001:end,:),9,1,2);
%         m1 = yc_mf(m1,2,1,2);
        m=drr3d_win(m1,0,100,0.004,2,4,0,128,24,1,0.5,0.5,0.5);
        m = m + (m1-m) * alpha;  % alpha=0.75/maxfold
        snr(iter)=yc_snr(dd(:,:,ir),m);
        fprintf('ir=%d/%d, current iterations %d/%d, SNR=%g\n',ir,nr,iter,niter,snr(iter));
    end
    dbb(:,:,ir)=m;
    fprintf('ir=%d/%d, SNR=%g\n',ir,nr,snr(iter));
end

figure;yc_imagesc([dd(:,:,1),bb(:,:,1),dbb(:,:,1)],99);
figure;yc_imagesc([dd(:,:,60),bb(:,:,60),dbb(:,:,60)],99);

figure;yc_imagesc([squeeze(dd(:,1,:)),squeeze(bb(:,1,:)),squeeze(dbb(:,1,:))],99);
figure;yc_imagesc([squeeze(dd(:,100,:)),squeeze(bb(:,100,:)),squeeze(dbb(:,100,:))],99);
figure;yc_imagesc([squeeze(dd(:,128,:)),squeeze(bb(:,128,:)),squeeze(dbb(:,128,:))],99);
figure;yc_imagesc([squeeze(dd(:,64,:)),squeeze(bb(:,64,:)),squeeze(dbb(:,64,:))],99);

% save('pgs_drr.mat','dbb','rt0','rt1','-v7.3')

%% 3D FFT
lamda=0.7;maxfold=2;alpha=0.75/maxfold;%alpha=1;

nt0=1024;
d=dd;
d_b_com=zeros(size(d));
d_b=[];
nh=nr;ny=1;ntpad=nt;nxpad=nx;nypad=ny;rt=rt1;
% rt = db_gstimes(maxfold,nt,nx,ny,ntpad,nxpad,nypad);
for ih=1:nh
    d_bt = yc_Gamma(d(:,:,ih),rt,nt,nx,ny,ntpad,nxpad,nypad);
    d_b=[d_b,d_bt];
    d_b_com(:,:,ih)= yc_Gammai(d_bt,rt,nt,nx,1,ntpad,nxpad,1);
end

%starting model is zero
d_mod=zeros(size(d_b));
m_n1=zeros(size(d_b_com));
m_up=zeros(size(d_b_com));
vmx = prctile(reshape(abs(fftn(d_b_com(1:nt0,:,:))),nt0*nx*nh,1),99.99)

niter=50;
for iter=1:niter
    
    %%% compute the data misfit
    d_res = d_b - d_mod;  %d_a-W_aW_p^{-1}m
    for ih=1:nh
        m_up(:,:,ih) = db_Gammai(d_res(:,ih),rt,nt,nx,1,ntpad,nxpad,1);
    end 
    m_n = m_n1; % m_n and m_{n+1}
    m_n1 = m_n + m_up * lamda;
    mprime  = fftn(m_n1(1:nt0,:,:));
    vlim = ((niter-iter)/niter) * vmx ;
    for it=1:nt0
        for ix=1:nx
            for ih=1:nh
                if ( abs(mprime(it,ix,ih)) < vlim)
                    mprime(it,ix,ih) = 0.0 ;
                end
            end 
        end
    end 
%    mprime=(abs(mprime(1:nt0,:,:))>vlim).*mprime(1:nt0,nx,nh);
    
    m_n1 = [real( ifftn(mprime));zeros(nt-nt0,nx,nh)] ;
    
    %%% m_{n+1} = (1-alpha)*m_n + alpha*A^{-1}TA(m_n+\lamda*F^{pseudoinverse}(d_obs-Fm_n))
    m_temp = m_n1 - m_n;
    m_n1 = m_n + m_temp * alpha;  % alpha=0.75/maxfold
    
    %%% compute SNRs
    snr3(iter)=10*log10(sum(sum(sum(d(1:nt0,:,:).^2)))/sum(sum(sum((m_n1(1:nt0,:,:)-d(1:nt0,:,:)).^2))));
    %%% reforward modeling
    for ih=1:nh
        d_mod(:,ih) = yc_Gamma(m_n1(:,:,ih),rt,nt,nx,1,ntpad,nxpad,1);
    end 
    fprintf('iter = %d/%d \n',iter,niter);
    fprintf('Current SNR = %g \n\n',snr3(iter));
end

figure;yc_imagesc([squeeze(d(:,64,:)),squeeze(d_b_com(:,64,:)),squeeze(m_n1(:,64,:))],99);
figure;yc_imagesc([squeeze(d(:,120,:)),squeeze(d_b_com(:,120,:)),-squeeze(d(:,120,:))+squeeze(d_b_com(:,120,:)),squeeze(m_n1(:,120,:))],99);
figure;yc_imagesc([squeeze(d(:,2,:)),squeeze(d_b_com(:,2,:)),-squeeze(d(:,2,:))+squeeze(d_b_com(:,2,:)),squeeze(m_n1(:,2,:))],99);

figure;yc_imagesc([squeeze(d(:,:,30)),squeeze(d_b_com(:,:,30)),squeeze(m_n1(:,:,30))],99);

ir=10;figure;yc_imagesc([squeeze(d(:,:,ir)),squeeze(d_b_com(:,:,ir)),squeeze(m_n1(:,:,ir))],99);

