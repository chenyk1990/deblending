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
load yc_fieldsr.mat
%% in this dataset
%there are two sources data3d(:,:,1:60) and data3d(:,:,61:120)
%each source contains 120 shots
%there are 60 receivers
%

d1=data3d(:,:,1);
d2=data3d(:,:,61);

figure;
subplot(1,2,1);db_imagesc(d1);
subplot(1,2,2);db_imagesc(d2);

h1=1:120;
h2=1:120;

dt=0.004;
t=[0:1500-1]*dt;
nt=1500;
nx=120;
%% apply db_dithering
randn('state',202122);
shift1=floor(0.1*randn(1,size(d1,2))/dt);   % shift of data1 to data2
shift2=-shift1;                             % shift of data2 to data1

d1shift=db_dither(d1,shift1);
d2shift=db_dither(d2,shift2);

figure;
subplot(1,2,1);imagesc(h1,t,d1shift);
subplot(1,2,2);imagesc(h2,t,d2shift);

%% blend
d1b=d1+d2shift;
d2b=d2+d1shift;
figure;
subplot(1,2,1);imagesc(h1,t,d1b);
subplot(1,2,2);imagesc(h2,t,d2b);

%%
rand('state',2021222324);
[n1,n2]=size(d1);
D1=zeros(nt,nx);
D2=zeros(nt,nx);
% delta0=delta;
delta=shift2;

del=delta;      %ground truth
d1b=d1+db_dither(d2,del);
d2b=d2+db_dither(d1,-del);

figure;
subplot(1,2,1);db_imagesc(d1b);
subplot(1,2,2);db_imagesc(d2b);

%% blending for all receivers
data3db=zeros(size(data3d));
for ir=1:60
    data3db(:,:,ir)=data3d(:,:,ir)+db_dither(data3d(:,:,ir+60),del);
    data3db(:,:,ir+60)=data3d(:,:,ir+60)+db_dither(data3d(:,:,ir),-del);
end

%% first receiver
mask1=ones(size(d1));
mask1=db_mutterv(mask1,1,42,4.1);

mask2=ones(size(d2));
mask2=db_mutterv(mask2,61,42,4.1);

bd=d1b;
%% first receiver
d1b=data3db(:,:,1);
d2b=data3db(:,:,61);
d1=data3d(:,:,1);
d2=data3d(:,:,61);
bd=d1b;

%% correct source time
delta0=delta;

%% for all receivers
data3ddb=zeros(size(data3db));
niter=10;
for ir=1:60
    mask1=ones(size(d1));
    mask1=db_mutterv(mask1,ir,42,4.1);
    
    mask2=ones(size(d2));
    mask2=db_mutterv(mask2,60+ir,42,4.1);
    d1b=data3db(:,:,ir);
    d2b=data3db(:,:,60+ir);
    d1=data3d(:,:,ir);
    d2=data3d(:,:,60+ir);
    
    D1=zeros(nt,nx);
    D2=zeros(nt,nx);
    for iter=1:niter
        D1T=db_dither(D1,-delta0);
        D2T=db_dither(D2,delta0);
        D1u = D1 + 0.5*(d1b-(D1+D2T));    % updated model
        D2u = D2 + 0.5*(d2b-(D1T+D2));    % updated model
        %         D1= seislet_denoise_2d(D1u,dip,'ps',3,0.1,2,3); %This is to apply seislet-based sparse inversion
        %         D2= seislet_denoise_2d(D2u,dip2,'ps',3,0.1,2,3);%This is to apply seislet-based sparse inversion
        D1u=D1u.*mask1; %This is to apply a mask operator (muting) to remove noise before the first arrivals
        D2u=D2u.*mask2; %This is to apply a mask operator (muting) to remove noise before the first arrivals
        D1=drr3d_win(D1u,0,80,0.004,2,4,0,100,20,1,0.5,0.5,0.5);%This is to apply damped rank-reduction filter as a shaping operator to remove blending noise
        D2=drr3d_win(D2u,0,80,0.004,2,4,0,100,20,1,0.5,0.5,0.5);%This is to apply damped rank-reduction filter as a shaping operator to remove blending noise
        D1=D1.*mask1;
        D2=D2.*mask2;
        
        D1=D1+(d1b-D1).*(1-mask2).*mask1;
        D2=D2+(d2b-D2).*(1-mask1).*mask2;
        
        snr2(iter)=10*log10(sum(sum(d1.*d1))/sum(sum((d1-D1).*(d1-D1))));
        snr22(iter)=10*log10(sum(sum(d2.*d2))/sum(sum((d2-D2).*(d2-D2))));
 
    end
    data3ddb(:,:,ir)=D1;
    data3ddb(:,:,ir+60)=D2;
    fprintf('ir=%d,snr2=%g,snr22=%g\n',ir,snr2(niter),snr22(niter));
end

for ir=10:10
    figure(1)
    subplot(1,2,1);db_imagesc([data3d(:,:,ir),data3db(:,:,ir),data3ddb(:,:,ir),data3db(:,:,ir)-data3ddb(:,:,ir)]);title('Common receiver gather (source 1)');
    subplot(1,2,2);db_imagesc([data3d(:,:,ir+60),data3db(:,:,ir+60),data3ddb(:,:,ir+60),data3db(:,:,ir+60)-data3ddb(:,:,ir+60)]);title('Common receiver gather (source 2)');
end

for is=10:10
    figure(2)
    subplot(1,2,1);db_imagesc([squeeze(data3d(:,is,1:60)),squeeze(data3db(:,is,1:60)),squeeze(data3ddb(:,is,1:60)),squeeze(data3db(:,is,1:60))-squeeze(data3ddb(:,is,1:60))]);title('Common shot gather (source 1)');
    subplot(1,2,2);db_imagesc([squeeze(data3d(:,is,61:120)),squeeze(data3db(:,is,61:120)),squeeze(data3ddb(:,is,61:120)),squeeze(data3db(:,is,61:120))-squeeze(data3ddb(:,is,61:120))]);title('Common shot gather (source 2)');
end

save db_drr.mat data3d data3db data3ddb delta0




%% plot 3D comparison
figure('units','normalized','Position',[0.2 0.4 0.4, 1],'color','w');
subplot(2,2,1);db_mada3d(data3d(:,:,1:60));title('Unblended','Fontsize',12,'fontweight','normal');caxis([-0.1,0.1]);text(-30,0,'a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
xlabel('Shot #','Fontsize',12,'fontweight','normal','Rotation',-10,'horizontalalignment','center','verticalalignment','middle');
ylabel('Receiver #','Fontsize',12,'fontweight','normal','Rotation',20,'horizontalalignment','center','verticalalignment','middle');
zlabel('Time (s)','Fontsize',12,'fontweight','normal');
set(gca,'Linewidth',1.5,'Fontsize',12,'Fontweight','normal');
zlim([0,5]);

% subplot(3,4,2);db_mada3d(data3d(:,:,61:120));title('Source 2','Fontsize',12,'fontweight','normal');caxis([-0.1,0.1]);
subplot(2,2,2);db_mada3d(data3db(:,:,1:60));title('Blended','Fontsize',12,'fontweight','normal');caxis([-0.1,0.1]);text(-30,0,'b)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
xlabel('Shot #','Fontsize',12,'fontweight','normal','Rotation',-10,'horizontalalignment','center','verticalalignment','middle');
ylabel('Receiver #','Fontsize',12,'fontweight','normal','Rotation',20,'horizontalalignment','center','verticalalignment','middle');
zlabel('Time (s)','Fontsize',12,'fontweight','normal');
set(gca,'Linewidth',1.5,'Fontsize',12,'Fontweight','normal');
% ax = gca;
% ax.XAxisLocation='origin';
zlim([0,5]);

subplot(2,2,3);db_mada3d(data3ddb(:,:,1:60));title('Deblended (DRR)','Fontsize',12,'fontweight','normal');caxis([-0.1,0.1]);text(-30,0,'c)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
xlabel('Shot #','Fontsize',12,'fontweight','normal','Rotation',-10,'horizontalalignment','center','verticalalignment','middle');
ylabel('Receiver #','Fontsize',12,'fontweight','normal','Rotation',20,'horizontalalignment','center','verticalalignment','middle');
zlabel('Time (s)','Fontsize',12,'fontweight','normal');
set(gca,'Linewidth',1.5,'Fontsize',12,'Fontweight','normal');
zlim([0,5]);

% subplot(3,4,2);db_mada3d(data3d(:,:,61:120));title('Source 2','Fontsize',12,'fontweight','normal');caxis([-0.1,0.1]);
subplot(2,2,4);db_mada3d(data3db(:,:,1:60)-data3ddb(:,:,1:60));title('Noise (DRR)','Fontsize',12,'fontweight','normal');caxis([-0.1,0.1]);text(-30,0,'d)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
xlabel('Shot #','Fontsize',12,'fontweight','normal','Rotation',-10,'horizontalalignment','center','verticalalignment','middle');
ylabel('Receiver #','Fontsize',12,'fontweight','normal','Rotation',20,'horizontalalignment','center','verticalalignment','middle');
zlabel('Time (s)','Fontsize',12,'fontweight','normal');
set(gca,'Linewidth',1.5,'Fontsize',12,'Fontweight','normal');
% ax = gca;
% ax.XAxisLocation='origin';
zlim([0,5]);

% Create textarrow
annotation(gcf,'textarrow',[0.782986111111111 0.739583333333334],...
[0.840478087649402 0.809933598937583],'String',{'Blending noise'},'linewidth',2,'fontweight','bold','color','r');

% Create arrow
annotation(gcf,'arrow',[0.784722222222222 0.795138888888889],...
[0.837822045152722 0.817277556440903],'linewidth',2,'color','r');
print(gcf,'-dpng','-r300','test_db_drr.png');


