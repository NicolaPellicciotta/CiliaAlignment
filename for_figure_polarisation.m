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

    %% plot 
  [G,idtime,idflow]=findgroups(Time,Flow);
  POL=accumarray(G,Pol,[],@median); 
  ePOL=accumarray(G,Pol,[],@std)./sqrt(accumarray(G,Pol,[],@numel)) ;
  NN=accumarray(G,N,[],@median); 
  eNN=accumarray(G,N,[],@std)./sqrt(accumarray(G,N,[],@numel)) ;
  
 
   
  ind_c= (idflow==0);
  ind_exp = ~ind_c; 
  ms=8;lw=1;
 errorbar(idtime(ind_c)+3,POL(ind_c),ePOL(ind_c),'k^','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
 errorbar(idtime(ind_exp)+3,POL(ind_exp),ePOL(ind_exp),'ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','r');
 % xlabel('Days','Interpreter','latex') ;ylabel('Polarisation','Interpreter','latex');
  xlim([4,30])
 x0=0;y0=0;width=350;height=300;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',13);
legend({'control','0.8 dyne/cm$^2$'},'Interpreter','latex')
%saveas(gcf,'/u/homes/np451/Dropbox/alignment with the flow/figure2/polarisation_time.pdf')
%%
figure()
 errorbar(idtime(ind_c),NN(ind_c),eNN(ind_c),'k^','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
 errorbar(idtime(ind_exp),NN(ind_exp),eNN(ind_exp),'ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','r');
%  xlabel('Days','Interpreter','latex') ;ylabel('density','Interpreter','latex');
  xlim([4,30])
 x0=0;y0=0;width=300;height=250;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',15);
saveas(gcf,'/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/figures/celldensity_time.pdf')


  
%   figure(2);
%   boxplot(Pol(ind1),GTime(ind1),'Colors','r');hold on;
%   boxplot(Pol(ind2),GTime(ind2),'Colors','k');
%     
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
x0=0;y0=0;width=300;height=300;
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
x0=0;y0=0;width=300;height=300;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',13);
saveas(gcf,'/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/figures/corr_frequency_function_time.pdf')
saveas(gcf,'/u/homes/np451/Dropbox/alignment with the flow/figure3/orientation_corr_time.pdf')




%%  plot for polarisation vs shear stress

clear all;
close all;

subdirJuly={...   
   '16.7.18/0/FL'...
   '29.7.18/0/FL'...
  }; 
subdirJuly = strcat('/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/July2018/PIV/',subdirJuly);
 
subdirJune={...   
   '13.6.18/0/FL'...
   '18.6.18/0/FL'...
   }; 
subdirJune = strcat('/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/June2018/PIV/',subdirJune);

subdirMarch={...   
'29.3.19/beads/40rpm/FL','29.3.19/beads/40rpm2/FL','29.3.19/beads/control/FL'...
'12.4.19/beads/40rpm/FL','12.4.19/beads/40rpm2/FL'...,'12.4.19/beads/control/FL'... 
}; 

subdirMarch = strcat('/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/',subdirMarch);

subdir=[subdirMarch,subdirJuly,subdirJune]

holes={...
    [1:2000,4900:10000],'none','none'...% [1:500,500:10000],'none','none'...
    'none',[3865:4075,6415:7195],...
    'none','none','none','none'...
    };

shear_ind=[];
pol=[];
good=[];
flow=[];
npo=[]
shear_um=1700;
for jj=1:numel(subdir)
   cd(subdir{jj}); load('Res3.mat');
   for fov=1:numel(Res.fov);

     
     
      allfolders= regexp(Res.insert_name,'/','split'); 
      flow_string=allfolders{end-1};
      if strcmp(allfolders{end-2},'beads'); time_string=allfolders{end-3};
      else time_string=allfolders{end-2};
      end

      if strcmp(time_string,'18.6.18')  
          shear_ind=cat(1,shear_ind,floor(1*(shear_um-1)/shear_um)+1);
      elseif strcmp(time_string,'13.6.18')
          shear_ind=cat(1,shear_ind,floor(2*(shear_um-1)/shear_um)+1);
      elseif strcmp(time_string,'16.7.18')
          shear_ind=cat(1,shear_ind,floor(3*(shear_um-1)/shear_um)+1);
       elseif strcmp(time_string,'29.7.18')
          shear_ind=cat(1,shear_ind,floor(4*(shear_um-1)/shear_um)+1);
      else
          shear_ind=cat(1,shear_ind,floor(Res.fov(fov).posx/shear_um)+1);
      end

      
      if strcmp(flow_string,'40rpm') | strcmp(flow_string,'40rpm2'); Res.flow= 40;
      elseif flow_string=='control' | flow_string=='0'; Res.flow= 0;
          
      end
      flow=cat(1,flow,Res.flow);
      nn= numel(Res.fov(fov).nu)./numel(Res.fov(fov).x);
      npo=cat(1,npo,nn);
   
      if  isa(holes{jj},'char'); if holes{jj}=='none';good=cat(1,good,1);end
      elseif any(Res.fov(fov).posx ==holes{jj}); good=cat(1,good,0); 
      else good=cat(1,good,1); 
      end
   
    %%%% the experiment with 40rpm2 have teh flow direction reverted
      if ~strcmp(flow_string,'40rpm2'); pol= cat(1,pol,(Res.fov(fov).Polx));
      elseif strcmp(flow_string,'40rpm2'); pol= cat(1,pol,-(Res.fov(fov).Polx));
      end
 
   end
end
    good=logical(good);
    Pol=pol(good);
    Npo=npo(good);
%    Time=time(good);
    Flow= flow(good);
    Shear_ind=shear_ind(good);
    %%%
    [G,idshear,idflow]=findgroups(Shear_ind,Flow);
    POL=accumarray(G,Pol,[],@median); 
    ePOL=accumarray(G,Pol,[],@std)./sqrt(accumarray(G,Pol,[],@numel)) ;
    NPO=accumarray(G,Npo,[],@median); 
    eNPO=accumarray(G,Npo,[],@std)./sqrt(accumarray(G,Npo,[],@numel)) ;

    
  ind_c= (idflow==0);
  ind_exp = ~ind_c; 
  ms=10;lw=2;
  shear_c= shear_calculate('tapered',40,idshear(ind_c)*shear_um);
  shear_exp= shear_calculate('tapered',40,idshear(ind_exp)*shear_um);

  
  %%%% polarisation vs shear stress
  
 errorbar(shear_c,POL(ind_c),ePOL(ind_c),'k^','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
 errorbar(shear_exp,POL(ind_exp),ePOL(ind_exp),'ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','r');
  xlabel('Shear stress [dyne/c$m^2$]','Interpreter','latex') ;ylabel('Polarisation','Interpreter','latex');
  %xlim([0,5])
 x0=0;y0=0;width=600;height=500;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',15);
saveas(gcf,'/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/figures/polarisation_shear.pdf')

    %%%%% numbero fo cells vs shear stress
figure()
    errorbar(shear_c,NPO(ind_c),eNPO(ind_c),'k^','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
 errorbar(shear_exp,NPO(ind_exp),eNPO(ind_exp),'ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','r');
 xlabel('Shear stress [dyne/c$m^2$]','Interpreter','latex') ;ylabel('Cell density','Interpreter','latex');
 %xlim([0,5])
 x0=0;y0=0;width=600;height=500;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',15);
saveas(gcf,'/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/figures/density_shear.pdf')


%% plot for seeing if polarisation correlate with synchronisation 


%%%%% for plot polarisation of fov
%%% directories with results
clear all
subdir={...
...            '22.3.19/beads/40rpm/FL','22.3.19/beads/40rpm2/FL','22.3.19/beads/control/FL'...
...              '29.3.19/beads/40rpm/FL'...
            '29.3.19/beads/40rpm2/FL'...
            '29.3.19/beads/control/FL'...
...             '5.4.19/beads/40rpm/FL','5.4.19/beads/40rpm2/FL','5.4.19/beads/control/FL'...
           '12.4.19/beads/40rpm/FL'...
          '12.4.19/beads/40rpm2/FL'...,
        }; 
%%%%%% hole for seeing the correlation lenght
holes={...
...   'none','none','none'...
...   [1:2000,4900:10000]...
    'none'...
   'none'...
...   'none','none','none'...
    'none'...
  [3865:4075,6415:7195]...
    };

subdir = strcat('/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/',subdir);
    
pol=[];good=[];time=[];flow=[];
l=[];
lf=[];df=[];
medf=[];nf=[];m=[];mv=[];
for jj=1:numel(subdir)
   cd(subdir{jj}); load('Res3.mat');
   for fov=1:numel(Res.fov);
   
   % determine the experiment and date    
   allfolders= regexp(Res.insert_name,'/','split'); 
   flow_string=allfolders{end-1};
   time_string=allfolders{end-3};
   
   if strcmp(flow_string,'40rpm') | strcmp(flow_string,'40rpm2'); Res.flow= 40;
   elseif flow_string=='control'; Res.flow= 0;
   end
   flow=cat(1,flow,Res.flow);
   m=cat(1,m,Res.fov(fov).Mm);
   mv=cat(1,mv,mean(Res.fov(fov).v1m( Res.fov(fov).ind)));
   
   if     strcmp(time_string,'22.3.19'); time_temp= 5;
   elseif strcmp(time_string,'29.3.19'); time_temp= 12;
   elseif strcmp(time_string,'5.4.19'); time_temp= 19;
   elseif strcmp(time_string,'12.4.19'); time_temp= 26;
   end
   time=cat(1,time,time_temp); Res.time=time_temp;
%   good=cat(1,good,Res.fov(fov).good);





       if  isa(holes{jj},'char'); if holes{jj}=='none';good=cat(1,good,1);end
      elseif any(Res.fov(fov).posx ==holes{jj}); good=cat(1,good,0); 
      else good=cat(1,good,1); 
      end
   


   %%%% the experiment with 40rpm2 have teh flow direction reverted
   if ~strcmp(flow_string,'40rpm2'); pol= cat(1,pol,(Res.fov(fov).Polx));
   elseif strcmp(flow_string,'40rpm2'); pol= cat(1,pol,(Res.fov(fov).Polx));
   end
%     if ~strcmp(flow_string,'40rpm2'); pol= cat(1,pol,(Res.fov(fov).Pol));
%     elseif strcmp(flow_string,'40rpm2'); pol= cat(1,pol,(Res.fov(fov).Pol));
%     end
    if isempty(Res.fov(fov).nu) | numel(Res.fov(fov).nu)<20
        l=cat(1,l,nan);
        lf=cat(1,lf,nan);
        good(end)=0;
    else   
        l=cat(1,l,Res.fov(fov).fit_out.b);
        lf=cat(1,lf,Res.fov(fov).freq.fit_out.b);
    end
    
   F=Res.fov(fov).freq.F4;
   df =cat(1,df,(nanstd(F(:))/nanmean(F(:))));
   medf= cat(1,medf,nanmedian(F(:)));
   nf= cat(1,nf,(numel(F(~isnan(F(:))))));
   
   end
     
end
close all
    good=logical(good);
    Pol=pol(good);
    Time=time(good);
    Flow= flow(good);
    L=l(good);
    Lf=lf(good);
    Df=df(good);
    Medf=medf(good);
    Nf=nf(good); 
    Mv=mv(good);
    M=m(good);
%%%% plot 
  [G,idtime,idflow]=findgroups(Time,Flow);
  POL=accumarray(G,Pol,[],@median); 
  ePOL=accumarray(G,Pol,[],@std)./sqrt(accumarray(G,Pol,[],@numel)) ;

  Lc=accumarray(G,L,[],@median); 
  eLc=accumarray(G,L,[],@std)./sqrt(accumarray(G,L,[],@numel)) ;
  
   Lcf=accumarray(G,Lf,[],@median); 
   eLcf=accumarray(G,Lf,[],@std)./sqrt(accumarray(G,Lf,[],@numel)) ;

   DF=accumarray(G,Df,[],@median); 
  eDF=accumarray(G,Df,[],@std)./sqrt(accumarray(G,Df,[],@numel)) ;

  
   ms=10;lw=2;
% 
ind_c= (idflow==0);
 ind_exp = ~ind_c; 
  ms=10;lw=2;
 errorbar(idtime(ind_c),Lc(ind_c),eLc(ind_c),'k^','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
 errorbar(idtime(ind_exp),Lc(ind_exp),eLc(ind_exp),'ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','r');
  xlabel('Days','Interpreter','latex') ;ylabel('Polarisation','Interpreter','latex');
 xlim([4,30])
 x0=0;y0=0;width=600;height=500;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',15);
saveas(gcf,'/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/figures/polarisation_time.pdf')

%%  plot polarisation and flow velocity

figure(5) 

x=Pol(Flow==0);y=Mv(Flow==0);n_bins=5;
plot(x,y,'k^','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
[histw,ehistw,m_bins]=hist_nico2(x,y,5,[]);
errorbar(m_bins,histw,ehistw,'-ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
%plot(L(Flow==0),Df(Flow==0),'ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
 xlabel('$Polarisation y$','Interpreter','latex') ;ylabel('Flow velocity [px/frame]','Interpreter','latex');
 title('Flow velocity [px/frame]');
 x0=0;y0=0;width=600;height=500;
set(gcf,'position',[x0,y0,width,height]);
set(gca,'FontSize',15);
%set(gca,'XScale', 'log', 'YScale', 'log');
saveas(gcf,'/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/figures/Pol_flowvelocity.pdf')


%%

figure(5) 

x=Nf(Flow==0);y=M(Flow==0);n_bins=5;
plot(x,y,'k^','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
[histw,ehistw,m_bins]=hist_nico2(x,y,5,[]);
errorbar(m_bins,histw,ehistw,'-ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
%plot(L(Flow==0),Df(Flow==0),'ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
 xlabel('$Polarisation y$','Interpreter','latex') ;ylabel('Flow velocity [px/frame]','Interpreter','latex');
 title('Flow velocity [px/frame]');
 x0=0;y0=0;width=600;height=500;
set(gcf,'position',[x0,y0,width,height]);
set(gca,'FontSize',15);
%set(gca,'XScale', 'log', 'YScale', 'log');
saveas(gcf,'/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/figures/Pol_flowvelocity.pdf')


%%
    
figure(5) 

x=Pol(Flow==40);y=Nf(Flow==40);n_bins=5;
plot(x,y,'k^','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
[histw,ehistw,m_bins]=hist_nico2(x,y,5,[]);
errorbar(m_bins,histw,ehistw,'-ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
%plot(L(Flow==0),Df(Flow==0),'ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
 xlabel('$\zeta_{or}$','Interpreter','latex') ;ylabel('$\zeta_{f}$','Interpreter','latex');
 title('correlation orientation lenght vs frequency length');
 x0=0;y0=0;width=600;height=500;
set(gcf,'position',[x0,y0,width,height]);
set(gca,'FontSize',15);
%set(gca,'XScale', 'log', 'YScale', 'log');
saveas(gcf,'/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/figures/corr_orientation_vs_freq.pdf')

figure(4)

x= Pol(Flow==40);y=Df(Flow==40);n_bins=5;
x= Pol;y=Df;n_bins=5;
plot(x,y,'k^','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
[histw,ehistw,m_bins]=hist_nico2(x,y,5,[]);
errorbar(m_bins,histw,ehistw,'-ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
xlabel('$\zeta_{or}$','Interpreter','latex') ;ylabel('$\sigma_{f}$[Hz]','Interpreter','latex');
 title('Polarisation vs frequency detuning');
 x0=0;y0=0;width=600;height=500;
set(gcf,'position',[x0,y0,width,height]);
set(gca,'FontSize',15);
saveas(gcf,'/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/figures/polarisation_vs_freq_detuning.pdf')

figure();
plot(Pol(Flow==40),Medf(Flow==40),'k^','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;

