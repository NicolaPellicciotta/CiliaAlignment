cd '/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/29.3.19/beads/40rpm2/FL'
load('Res3.mat');

%% for fifure with the flow velocity
XX=Res.xx,YY=Res.yy;nu=Res.nu;nv=Res.nv
ppm=0.24
figure()
quiver(XX*ppm,YY*ppm,Res.uu1,Res.vv1,0.8,'k');
xlabel('X [um]');
ylabel(' Y[um]');

axis equal
ylabel('correlation function','Interpreter','latex');
xlabel('distance $r$ [$\mu$m]','Interpreter','latex');
x0=0;y0=0;width=400;height=800;
ylim([-3000,-1500])
xlim([100,800])
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',19);

%% for fifure for orientation fluctuation
ppm=0.24;
XX=Res.xx,YY=Res.yy;nu=Res.nu;nv=Res.nv;uu1=Res.uu1;vv1=Res.vv1
du=uu1-nanmean(uu1);dv=vv1-nanmean(vv1);
figure()
quiver(XX*ppm,YY*ppm,du,dv,1,'k');
xlabel('X [um]');
ylabel(' Y[um]');

axis equal
ylabel('correlation function','Interpreter','latex');
xlabel('distance $r$ [$\mu$m]','Interpreter','latex');
x0=0;y0=0;width=400;height=800;
ylim([-3000,-1500])
xlim([100,800])
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',19);
