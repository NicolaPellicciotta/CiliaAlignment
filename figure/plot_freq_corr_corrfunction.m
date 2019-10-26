%%% correlation fucntion plots
cd '/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/29.3.19/beads/40rpm2/FL'
load('Res3.mat');
% frequency 
Ncf=400;ppm=0.24;
Nfov=numel(Res.fov);
ccft=zeros([Ncf,Nfov]);
for ii=1:Nfov
ccft(:,ii)=Res.fov(ii).freq.cc(1:Ncf);
end
ccf=nanmean(ccft,2);
idr=Res.fov(1).freq.idr*ppm;
    ft = fittype('a*exp(-x/b)+c','independent','x');
    fo = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[0,idr(2),-1],...
    'Upper',[1,Inf,1],...
    'StartPoint',[0.9, idr(2), 0]);
 
    fit_bins= idr; fit_cc=ccf;max_cc=30;
    fit_out = fit(fit_bins(1:max_cc),fit_cc(1:max_cc),ft,fo);

close all
figure();
ppm=0.24;
mso=10;lwo=2;
msf=4;lwf=0.5;fov=10;
%plot(Res.idr*ppm,Res.cc,'ro','MarkerSize',mso,'LineWidth',lwo,'MarkerFaceColor','r');hold on;
plot(idr(1:Ncf),ccf,'ro','MarkerSize',msf,'LineWidth',lwf,'MarkerFaceColor','r');
hold on;
%plot(fit_out)
xlim([0,40]);



cd '/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/29.3.19/beads/control/FL'
load('Res3.mat');
% frequency 
Ncf=400;ppm=0.24;
Nfov=numel(Res.fov);
ccft=zeros([Ncf,Nfov]);
for ii=1:Nfov
ccft(:,ii)=Res.fov(ii).freq.cc(1:Ncf);
end
ccf=nanmean(ccft,2);
idr=Res.fov(1).freq.idr*ppm;

    ft = fittype('exp(-x/b)','independent','x');
    fo = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[idr(2)],...
    'Upper',[Inf],...
    'StartPoint',[ idr(2)]);
 
    fit_bins= idr; fit_cc=ccf;max_cc=100;
    fit_out = fit(fit_bins(1:max_cc),fit_cc(1:max_cc),'exp1');

ppm=0.24;
plot(idr(1:Ncf),ccf,'k^','MarkerSize',msf,'LineWidth',lwf,'MarkerFaceColor','k');
hold on;
plot(idr(1:Ncf),fit_out(1:Ncf),'b-','LineWidth',1.5)
xlim([0,40]); ylim([-0.2,1.1])
%ylabel('correlation function','Interpreter','latex');
%xlabel('distance $r$ [$\mu$m]','Interpreter','latex');
x0=0;y0=0;width=250;height=200;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',13);
legend({'0.8 dyne/cm$^2$','control','$y=e^{-r/ \xi_{f}}$'},'Interpreter','latex');

saveas(gcf,'/u/homes/np451/Dropbox/alignment with the flow/figure3/freq_corr_single.pdf')


%set(gca,'XScale', 'log', 'YScale', 'log');
%%%%%%%%% figure orientation correlation single 
%%
close all
cd '/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/29.3.19/beads/40rpm2/FL'
load('Res3.mat');

ppm=0.24;
mso=4;lwo=0.5;
msf=5;lwf=1;fov=10;
figure();
plot(Res.idr*ppm,Res.cc,'ro','MarkerSize',mso,'LineWidth',lwo,'MarkerFaceColor','r');hold on;
%plot(Res.fov(fov).freq.idr(1:Ncf)*ppm,ccf,'bo','MarkerSize',msf,'LineWidth',lwf,'MarkerFaceColor','b')
xlim([0.5,300]);ylim([-0.1,1.1]);
%ylabel('correlation function','Interpreter','latex');
%xlabel('distance $r$ [$\mu$m]','Interpreter','latex');
fit_bins1= Res.idr*ppm; fit_cc1=Res.cc;max_cc=50;
fit_out_1 = fit(fit_bins1(1:max_cc),fit_cc1(1:max_cc),'exp1');

%set(gca,'XScale', 'log', 'YScale', 'log');

cd '/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/29.3.19/beads/control/FL'
load('Res3.mat');
fit_bins= Res.idr*ppm; fit_cc=Res.cc;max_cc=50;
fit_out_2 = fit(fit_bins(1:max_cc),fit_cc(1:max_cc),'exp1');

plot(Res.idr*ppm,Res.cc,'k^','MarkerSize',mso,'LineWidth',lwo,'MarkerFaceColor','k');hold on;
plot(Res.idr*ppm,fit_out_2(Res.idr*ppm),'b-','LineWidth',1);
plot(fit_bins1,fit_out_1(fit_bins1),'b-','LineWidth',1)
x0=0;y0=0;width=300;height=200;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',13);


legend({'0.8 dyne/cm$^2$','control','$y=e^{-r/ \xi_{o}}$'},'Interpreter','latex');
%set(gca,'XScale', 'log', 'YScale', 'log');

saveas(gcf,'/u/homes/np451/Dropbox/alignment with the flow/figure3/orientation_corr_single.pdf')

saveas(gcf,'/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/figures/correlations_orinetation.pdf')


%% compare correlation fucntion with days

%%%%% for plot polarisation of fov
%%% directories with results
clear all
subdirMarch={...   
   '22.3.19/beads/40rpm/FL','22.3.19/beads/40rpm2/FL','22.3.19/beads/control/FL'...
    '29.3.19/beads/40rpm/FL','29.3.19/beads/40rpm2/FL','29.3.19/beads/control/FL'...
    '5.4.19/beads/40rpm/FL','5.4.19/beads/40rpm2/FL'...
...    '5.4.19/beads/control/FL'...
    '12.4.19/beads/40rpm/FL','12.4.19/beads/40rpm2/FL'...,'12.4.19/beads/control/FL'... 
}; 
subdirMarch = strcat('/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/',subdirMarch);

subdirJuly={...   
...  '12.7.18/5A/FL', '12.7.18/5B/FL', '12.7.18/5C/FL'...
...   '16.7.18/5A/FL', '16.7.18/5B/FL', '16.7.18/5C/FL'...
   '16.7.18/0/FL'...
...  '12.7.18/15A/FL', '12.7.18/15B/FL', '12.7.18/15C/FL'...
   '29.7.18/0/FL'...
...   , '29.7.18/5A/FL', '29.7.18/5B/FL'...
  }; 
subdirJuly = strcat('/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/July2018/PIV/',subdirJuly);
 
subdirJune={...   
...   '13.6.18/30A/FL', '13.6.18/30B/FL'...
   '13.6.18/0/FL'...
...   '13.6.18/100A/FL', '13.6.18/1000/FL'...
   '18.6.18/0/FL'...
...   , '18.6.18/100/FL', '18.6.18/200/FL'...

   }; 
subdirJune = strcat('/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/June2018/PIV/',subdirJune);

subdir=[subdirMarch,subdirJuly,subdirJune]

pol=[];
good=[];
time=[];
flow=[];
l=[];
lf=[];
n=[];
for jj=1:numel(subdir)
   cd(subdir{jj}); load('Res3.mat');
   for fov=1:numel(Res.fov);
   
   % determine the experiment and date    
   allfolders= regexp(Res.insert_name,'/','split'); 
   flow_string=allfolders{end-1};
   if strcmp(allfolders{end-2},'beads'); time_string=allfolders{end-3};
   else time_string=allfolders{end-2};
   end

   
   if strcmp(flow_string,'40rpm') | strcmp(flow_string,'40rpm2'); Res.flow= 40;
   elseif strcmp(flow_string,'control') | strcmp(flow_string,'0'); Res.flow= 0;
   else Res.flow=nan;
   end
   flow=cat(1,flow,Res.flow);
   
   if     strcmp(time_string,'22.3.19') | strcmp(time_string,'12.7.18') | strcmp(time_string,'13.6.18'); 
       time_temp= 5;
   elseif strcmp(time_string,'29.3.19') | strcmp(time_string,'16.7.18');
        time_temp= 12;
   elseif strcmp(time_string,'5.4.19') | strcmp(time_string,'18.6.18') ; 
       time_temp= 19;
   elseif strcmp(time_string,'12.4.19')| strcmp(time_string,'29.7.18'); 
       time_temp= 26;
   end
   
   time=cat(1,time,time_temp); Res.time=time_temp;
   good=cat(1,good,Res.fov(fov).good);

     n=cat(1,n, numel(Res.fov(fov).nv));
   

   %%%% the experiment with 40rpm2 have teh flow direction reverted
   if ~strcmp(flow_string,'40rpm2'); pol= cat(1,pol,(Res.fov(fov).Polx));
   elseif strcmp(flow_string,'40rpm2'); pol= cat(1,pol,-(Res.fov(fov).Polx));
   end
     if isempty(Res.fov(fov).nu) | numel(Res.fov(fov).nu)<20
        l=cat(1,l,nan);
        lf=cat(1,lf,nan);
        good(end)=0;
    else   
        l=cat(1,l,Res.fov(fov).fit_out.b);
        lf=cat(1,lf,Res.fov(fov).freq.fit_out.b);
    end
    
   end
end
    good=logical(good);
    Pol=pol(good);
    Time=time(good);
    Flow= flow(good);
    L=l(good);
    Lf=lf(good);
    N=n(good).*(32^2)/(1200*1920);

  [G,idtime,idflow]=findgroups(Time,Flow);
  POL=accumarray(G,Pol,[],@median); 
  ePOL=accumarray(G,Pol,[],@std)./sqrt(accumarray(G,Pol,[],@numel)) ;
  NN=accumarray(G,N,[],@median); 
  eNN=accumarray(G,N,[],@std)./sqrt(accumarray(G,N,[],@numel)) ;
  ind_c= (idflow==0);
  ind_exp = ~ind_c; 

    %%     plot frquecny correlation function compare
ppm=0.24;
close all
ms=8; lw=1;
Lc=accumarray(G,L,[],@median)*ppm; 
eLc=accumarray(G,L,[],@std)./sqrt(accumarray(G,L,[],@numel))*ppm ;
  
Lcf=accumarray(G,Lf,[],@median)*ppm; 
eLcf=accumarray(G,Lf,[],@std)./sqrt(accumarray(G,Lf,[],@numel))*ppm ;

figure() 
errorbar(idtime(ind_c)+3,Lcf(ind_c),eLcf(ind_c),'k^','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
errorbar(idtime(ind_exp)+3,Lcf(ind_exp),eLcf(ind_exp),'ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','r');
%xlabel('Days','Interpreter','latex') ;
%ylabel('Freq Correlation lenght [$\mu$m]','Interpreter','latex');
 
xlim([4,30]);%title('polarisation correlation functions');
legend({'control','0.8 dyne/cm$^2$'},'Interpreter','latex')
x0=0;y0=0;width=300;height=200;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',13);
saveas(gcf,'/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/figures/corr_frequency_function_time.pdf')
saveas(gcf,'/u/homes/np451/Dropbox/alignment with the flow/figure3/freq_corr_time.pdf')

%% for orientation

figure() 
errorbar(idtime(ind_c)+3,Lc(ind_c),eLc(ind_c),'k^','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
errorbar(idtime(ind_exp)+3,Lc(ind_exp),eLc(ind_exp),'ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','r');
%xlabel('Days','Interpreter','latex') ;
%ylabel('Freq Correlation lenght [$\mu$m]','Interpreter','latex');
 
xlim([4,31]);%title('polarisation correlation functions');
legend({'control','0.8 dyne/cm$^2$'},'Interpreter','latex')
x0=0;y0=0;width=300;height=200;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',13);
saveas(gcf,'/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/figures/corr_frequency_function_time.pdf')
saveas(gcf,'/u/homes/np451/Dropbox/alignment with the flow/figure3/orientation_corr_time.pdf')





