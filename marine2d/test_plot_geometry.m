%5 km x 5 km

gx_real=linspace(2.5,5,30);
gy_real=2.5*ones(size(gx_real));

sx_real=linspace(0,5,100);nn=length(sx_real)/2;
sx2_real=linspace(2.5,7.5,100);nn2=length(sx2_real)/2;
sy_real=[0*ones(1,nn),0*ones(1,nn)];
sy2_real=[5*ones(1,nn2),5*ones(1,nn2)];
figure;
title('Source and receiver geometry' ,'FontSize', 18, 'fontweight', 'b')
hold on
% line(gx_real-11500,gy_real,'color','b','linewidth',3,'linestyle','none','marker','v');
plot(gx_real,gy_real,'color','b','linewidth',3,'linestyle','none','marker','v');
plot(gx_real(1),gy_real(1),'color','m','linewidth',3,'linestyle','none','marker','v');

sx=sx_real(1:4:end);sy=sy_real(1:4:end);
sx2=sx2_real(1:4:end);sy2=sy2_real(1:4:end);
line(sx,sy,'color','r','linewidth',4,'linestyle','none','marker','*')
line(sx2,sy2,'color','g','linewidth',4,'linestyle','none','marker','*')

%% add second y and x axis
[hAx, hP1,hP2 ]=plotyy(zeros(1,10),zeros(1,10),...
                       zeros(1,10),zeros(1,10));
axes(hAx(1));
set(gca,'XTick',zeros(1,0),'YColor','k','YTick',zeros(1,0))
ylim([min(sy_real) max(sy2_real)])
xlim([min(sx_real) max(sx2_real)])
xlabel('Inline ','FontSize', 16, 'fontweight', 'b')
ylabel('Crossline','FontSize', 16, 'fontweight', 'b')

axes(hAx(2));
set(gca,'XTick',zeros(1,0),'YColor','k','YTick',zeros(1,0))
ylim([min(sy_real) max(sy2_real)])
xlim([min(sx_real) max(sx2_real)])
xlabel('Inline ','FontSize', 16, 'fontweight', 'b')

print(gcf,'-depsc','-r300','fig_geom.eps');
%% adding arrows
h1=annotation('arrow',[0.15,0.2],[0.16 0.16],'Linewidth',2,'color','r') %[begx,endx],[begy,endy]
h2=annotation('arrow',[0.4,0.45],[0.86 0.86],'Linewidth',2,'color','g') %[begx,endx],[begy,endy]

print (gcf, '-dpng', '-r300','marine2d_geometry.png');
