%%%%% in this script we make figure 
% 1)for comparing few correlation functions
% 2) all the correlation lenght vs viscosity
% 3) viscosity vs concentration of MC  
%ciao bello!!!ciaooo ! <====3 
% e' un pene? RSPV by Roberta
%if it's yours it's more like <=3 by Sho...
%%% SHO!??? SHO is there ? 


subdir={...   
 'v0_0/FL','v0.5_5/FL','v1_4/FL','v1_5/FL','v2_5/FL'}
 
... 'v0_1/FL','v0.5_4/FL'...
... 'v0.5_5/FL','v1_4/FL','v1_5/FL'...
... 'v1.5_5/FL','v1.5_4/FL','v2_5/FL','v2_4/FL'...
...}; 
subdir = strcat('/media/np451/Seagate Backup Plus Drive1/DATA/viscosity_ependymal/PIV-10.8/',subdir)

visc=[0,0.5,1,1.5,2]
%     0,1,0.5...
%     0.5,1,1 ...
%     1.5,1.5,2,2 ...
% ]


%%%% orientation
cc=1;
figure()
ppm=0.14;
ms=30;lw=0.5;
for jj=1:2:numel(subdir)
   allfolders= regexp(Res.insert_name,'/','split'); 
   visc_string=allfolders{end-1};
    cleg{cc}=strcat(num2str(visc(jj)),' %MC');cc=cc+1;
   cd(subdir{jj}); load('Res3.mat');
    errorbar(Res.idr*ppm,Res.cc,Res.ecc,'.','MarkerSize',ms,'LineWidth',lw);hold on;
end
xlim([0,300]);
legend(cleg,'Interpreter','latex');x0=0;y0=0;width=400;height=300;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',13);
ylim([-0.1,1.1]);
%set(gca,'XScale', 'log', 'YScale', 'log');
saveas(gcf,'/media/np451/Seagate Backup Plus Drive1/DATA/viscosity_ependymal/DDM_tot/Figures_sigmoids/orientation_corr_plot_all.pdf')
%%
%%%% frequency

cc=1;
figure()
ppm=0.14;
ms=30;lw=0.5;
for jj=1:2:numel(subdir)
   allfolders= regexp(Res.insert_name,'/','split'); 
   visc_string=allfolders{end-1};
    cleg{cc}=strcat(num2str(visc(jj)),' %MC');cc=cc+1;
   cd(subdir{jj}); load('Res3.mat');
    errorbar(Res.freq.idr*ppm,Res.freq.cc,Res.freq.ecc,'.','MarkerSize',ms,'LineWidth',lw);hold on;
end
xlim([0,80]);
legend(cleg,'Interpreter','latex');x0=0;y0=0;width=400;height=300;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',13);
ylim([-0.1,1.1]);
%set(gca,'XScale', 'log', 'YScale', 'log');
saveas(gcf,'/media/np451/Seagate Backup Plus Drive1/DATA/viscosity_ependymal/DDM_tot/Figures_sigmoids/frequency_corr_plot_all.pdf')

%% compare correlation fucntion with days for viscosity data

%%%%% for plot polarisation of fov
%%% directories with results
clear all
subdir={...   
 'v0_0/FL','v0_1/FL','v0.5_4/FL'...
 'v0.5_5/FL','v1_4/FL','v1_5/FL'...
 'v1.5_5/FL','v1.5_4/FL','v2_5/FL','v2_4/FL'...
}; 
subdir = strcat('/media/np451/Seagate Backup Plus Drive1/DATA/viscosity_ependymal/PIV-10.8/',subdir)

visc=[
    0,0,0.5...
    0.5,1,1 ...
    1.5,1.5,2,2 ...
]

viscosity=[];
l=[];
lf=[];
nn=[];good=[];n=[];
for jj=1:numel(subdir)
   cd(subdir{jj}); load('Res3.mat');
   
     good=cat(1,good,1);
   ft = fittype('a*exp(-x/b)+c','independent','x');
    fo = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[0,1,-1],...
    'Upper',[1,Inf,1],...
    'StartPoint',[0.9, Res.idr(2)*10, 0]);
 
    fit_bins= Res.idr; fit_cc= Res.cc;max_cc=40;%floor(numel(cc(idth==0))/2);
    fit_out = fit(fit_bins(1:max_cc),fit_cc(1:max_cc),ft,fo);
   fit_bins= Res.freq.idr; fit_cc= Res.freq.cc;max_cc=10;%floor(numel(cc(idth==0))/2);
    fit_out2 = fit(fit_bins(1:max_cc),fit_cc(1:max_cc),ft,fo);

    
     l=cat(1,l,fit_out.b);
     lf=cat(1,lf,fit_out2.b);    
     
     viscosity=cat(1,viscosity,visc(jj));
     
    
end 
good=logical(good)
L=l(good);
Lf=lf(good);
Viscosity=viscosity(good);

  [G,idvisc]=findgroups(Viscosity);

%     plot frquecny correlation function compare
ppm=0.14;
close all
ms=8; lw=1;
Lc=accumarray(G,L,[],@median)*ppm; 
eLc=accumarray(G,L,[],@std)./sqrt(accumarray(G,L,[],@numel))*ppm ;
  
Lcf=accumarray(G,Lf,[],@median)*ppm; 
eLcf=accumarray(G,Lf,[],@std)./sqrt(accumarray(G,Lf,[],@numel))*ppm ;

%NN=accumarray(G,N,[],@median); 
%eNN=accumarray(G,N,[],@std)./sqrt(accumarray(G,N,[],@numel)) ;
%%
figure() 
errorbar(idvisc,Lcf,eLcf,'r^','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','r');hold on;
%errorbar(idvisc,Lcf,eLcf,'ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','r');
%xlabel('Days','Interpreter','latex') ;
%ylabel('Freq Correlation lenght [$\mu$m]','Interpreter','latex');
 
%xlim([4,30]);%title('polarisation correlation functions');
%legend({'control','0.8 dyne/cm$^2$'},'Interpreter','latex')
x0=0;y0=0;width=300;height=200;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',13);
%saveas(gcf,'/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/figures/corr_frequency_function_time.pdf')
saveas(gcf,'/media/np451/Seagate Backup Plus Drive1/DATA/viscosity_ependymal/DDM_tot/Figures_sigmoids/frequency_corr_tot.pdf')
saveas(gcf,'/u/homes/np451/Dropbox/alignment with the flow/figure_viscosity/frequency_corr_tot.pdf')

%% for orientation

figure() 
errorbar(idvisc,Lc,eLc,'bo','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','b');hold on;
%errorbar(idvisc,Lc,eLc,'ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','r');
%xlabel('Days','Interpreter','latex') ;
%ylabel('Freq Correlation lenght [$\mu$m]','Interpreter','latex');
 
%xlim([4,31]);%title('polarisation correlation functions');
%legend({'control','0.8 dyne/cm$^2$'},'Interpreter','latex')
x0=0;y0=0;width=300;height=200;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',13);
%saveas(gcf,'/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/figures/corr_frequency_function_time.pdf')
saveas(gcf,'/media/np451/Seagate Backup Plus Drive1/DATA/viscosity_ependymal/DDM_tot/Figures_sigmoids/orientation_corr_tot.pdf')
saveas(gcf,'/u/homes/np451/Dropbox/alignment with the flow/figure_viscosity/orientation_corr_tot.pdf')



%% compare correlation fucntion with days for viscosity data

%%%%% for plot polarisation of fov
%%% directories with results
clear all
subdir={...   
 'v0_0/FL','v0_1/FL','v0.5_4/FL'...
 'v0.5_5/FL','v1_4/FL','v1_5/FL'...
 'v1.5_5/FL','v1.5_4/FL','v2_5/FL','v2_4/FL'...
}; 
subdir = strcat('/media/np451/Seagate Backup Plus Drive1/DATA/viscosity_ependymal/PIV-10.8/',subdir)

visc=[
    0,0,0.5...
    0.5,1,1 ...
    1.5,1.5,2,2 ...
]

viscosity=[];
l=[];
lf=[];
nn=[];good=[];n=[];
for jj=1:numel(subdir)
   cd(subdir{jj}); load('Res3.mat');
   for fov=1:numel(Res.fov);
   

    good=cat(1,good,1)
     if isempty(Res.fov(fov).nu) | numel(Res.fov(fov).nu)<10
        l=cat(1,l,nan);
        lf=cat(1,lf,nan);
        good(end)=0;
    else   
        l=cat(1,l,Res.fit_out.b);
        lf=cat(1,lf,Res.fov(fov).freq.fit_out.b);
     end
     
     
     viscosity=cat(1,viscosity,visc(jj));
     
   d_lim=180;
   fXX=Res.fov(fov).x(Res.fov(fov).ind);
   fYY=Res.fov(fov).y(Res.fov(fov).ind);
   n= cat(1,n,numel(fXX));
% 
%    
%    time_temp=ones([numel(fXX),1])*time_temp;
%    time=cat(1,time,time_temp);

   near= find_nearcilia(fXX(:),fYY(:),d_lim);
   nn= cat(1,nn,near(:));
    
   end
end
good=logical(good);
L=l(good);
Lf=lf(good);
N=n(good);
  N=N.*(32^2)/(1200*1920);
  Viscosity=viscosity(good);

  [G,idvisc]=findgroups(Viscosity);

%     plot frquecny correlation function compare
ppm=0.14;
close all
ms=8; lw=1;
Lc=accumarray(G,L,[],@median)*ppm; 
eLc=accumarray(G,L,[],@std)./sqrt(accumarray(G,L,[],@numel))*ppm ;
  
Lcf=accumarray(G,Lf,[],@median)*ppm; 
eLcf=accumarray(G,Lf,[],@std)./sqrt(accumarray(G,Lf,[],@numel))*ppm ;

NN=accumarray(G,N,[],@median); 
eNN=accumarray(G,N,[],@std)./sqrt(accumarray(G,N,[],@numel)) ;
%%
figure() 
errorbar(idvisc,Lcf,eLcf,'k^','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
%errorbar(idvisc,Lcf,eLcf,'ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','r');
%xlabel('Days','Interpreter','latex') ;
%ylabel('Freq Correlation lenght [$\mu$m]','Interpreter','latex');
 
%xlim([4,30]);%title('polarisation correlation functions');
%legend({'control','0.8 dyne/cm$^2$'},'Interpreter','latex')
x0=0;y0=0;width=300;height=200;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',13);
%saveas(gcf,'/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/figures/corr_frequency_function_time.pdf')
saveas(gcf,'/media/np451/Seagate Backup Plus Drive1/DATA/viscosity_ependymal/DDM_tot/Figures_sigmoids/frequency_corr_tot.pdf')

%% for orientation

figure() 
errorbar(idvisc,Lc,eLc,'k^','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
%errorbar(idvisc,Lc,eLc,'ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','r');
%xlabel('Days','Interpreter','latex') ;
%ylabel('Freq Correlation lenght [$\mu$m]','Interpreter','latex');
 
%xlim([4,31]);%title('polarisation correlation functions');
%legend({'control','0.8 dyne/cm$^2$'},'Interpreter','latex')
x0=0;y0=0;width=300;height=200;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',13);
%saveas(gcf,'/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/figures/corr_frequency_function_time.pdf')
saveas(gcf,'/media/np451/Seagate Backup Plus Drive1/DATA/viscosity_ependymal/DDM_tot/Figures_sigmoids/orientation_corr_tot.pdf')

%% for density

figure() 
errorbar(idvisc,NN,eNN,'k^','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
%errorbar(idvisc,Lc,eLc,'ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','r');
%xlabel('Days','Interpreter','latex') ;
%ylabel('Freq Correlation lenght [$\mu$m]','Interpreter','latex');
 
%xlim([4,31]);%title('polarisation correlation functions');
%legend({'control','0.8 dyne/cm$^2$'},'Interpreter','latex')
x0=0;y0=0;width=300;height=200;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',13);
%saveas(gcf,'/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/figures/corr_frequency_function_time.pdf')
saveas(gcf,'/u/homes/np451/Dropbox/alignment with the flow/figure_viscosity/density_viscosity.pdf')

saveas(gcf,'/media/np451/Seagate Backup Plus Drive1/DATA/viscosity_ependymal/DDM_tot/Figures_sigmoids/density_viscosity.pdf')
%%
%%%metile cellulose vs viscosity
clear all
cd '/media/np451/Seagate Backup Plus Drive1/DATA/viscosity_ependymal'
load('MU.mat');
figure() ; ms=13;lw=0.5;


%errorbar(idvisc,Lc,eLc,'ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','r');
%xlabel('Days','Interpreter','latex') ;
%ylabel('Freq Correlation lenght [$\mu$m]','Interpreter','latex');
p=polyfit(MT(2:end),log(MU(2:end)-MU(1)),1);
muat2= MU(1)+exp(polyval(p,2));
plot(MT,MU,'mo','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','m');hold on;
hold on;
plot(2,muat2,'ko','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;

plot(0:0.1:2.1 ,MU(1)+exp(polyval(p,0:0.1:2.1)),'k-','LineWidth',2 )


% figure()
% plot(log(MT(2:end)),log(MU(2:end)),'o');hold on;
% plot(log(MT),(polyval(p,log(MT))),'-' )


xlim([0,2.1]);
%title('polarisation correlation functions');
legend({'data','$y=e^{\alpha x}$'},'Interpreter','latex')
x0=0;y0=0;width=350;height=300;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',13);
%set(gca,'XScale', 'log', 'YScale', 'log');
%saveas(gcf,'/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/figures/corr_frequency_function_time.pdf')
saveas(gcf,'/u/homes/np451/Dropbox/alignment with the flow/figure_viscosity/MC_viscosity.pdf')


%%
%%%% with old data extrapolated from graph
MT=[0,0.5,1,1.5,2];
MU=[0.0010,0.0018,0.0029,0.0038,0.0062]
p=polyfit(MT(2:end),log(MU(2:end)-MU(1)),1);
muat2= MU(1)+exp(polyval(p,2));
plot(MT,MU,'mo','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','m');hold on;
hold on;
%plot(2,muat2,'ko','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;

plot(0:0.1:2.1 ,MU(1)+exp(polyval(p,0:0.1:2.1)),'k-','LineWidth',2 )


% figure()
% plot(log(MT(2:end)),log(MU(2:end)),'o');hold on;
% plot(log(MT),(polyval(p,log(MT))),'-' )


xlim([0,2.1]);
%title('polarisation correlation functions');
legend({'data','$y=e^{\alpha x}$'},'Interpreter','latex')
x0=0;y0=0;width=350;height=300;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',13);
%set(gca,'XScale', 'log', 'YScale', 'log');
%saveas(gcf,'/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/figures/corr_frequency_function_time.pdf')
saveas(gcf,'/u/homes/np451/Dropbox/alignment with the flow/figure_viscosity/MC_viscosity_olddata.pdf')
%% figure for flow pattern
cd '/media/np451/Seagate Backup Plus Drive1/DATA/viscosity_ependymal/PIV-10.8/v2_5/FL'
load('Res3.mat');
figure();
quiver(Res.xx*0.14,Res.yy*0.14,Res.nu,Res.nv,1);
xlabel('X [um]');
ylabel(' Y[um]');
title('Flow pattern with validation');
axis image
x0=0;y0=0;width=400;height=400;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',13);
%set(gca,'XScale', 'log', 'YScale', 'log');
%saveas(gcf,'/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/figures/corr_frequency_function_time.pdf')
saveas(gcf,'/u/homes/np451/Dropbox/alignment with the flow/figure_viscosity/flow_v2_5.pdf')
%
