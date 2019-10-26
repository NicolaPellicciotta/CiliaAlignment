%%%%% for plot polarisation of fov
%%% directories with results
clear all
subdir={...   
   '22.3.19/beads/40rpm/FL','22.3.19/beads/40rpm2/FL','22.3.19/beads/control/FL'...
    '29.3.19/beads/40rpm/FL','29.3.19/beads/40rpm2/FL','29.3.19/beads/control/FL'...
    '5.4.19/beads/40rpm/FL','5.4.19/beads/40rpm2/FL','5.4.19/beads/control/FL'...
    '12.4.19/beads/40rpm/FL','12.4.19/beads/40rpm2/FL'...,'12.4.19/beads/control/FL'... 
}; 

subdir = strcat('/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/',subdir);
    
pol=[];
good=[];
time=[];
flow=[];
l_x=[];l_y=[];
lf_x=[];lf_y=[];
for jj=1:numel(subdir)
   cd(subdir{jj}); load('Res.mat');
   for fov=1:numel(Res.fov);
   
   % determine the experiment and date    
   allfolders= regexp(Res.insert_name,'/','split'); 
   flow_string=allfolders{end-1};
   time_string=allfolders{end-3};
   
   if strcmp(flow_string,'40rpm') | strcmp(flow_string,'40rpm2'); Res.flow= 40;
   elseif flow_string=='control'; Res.flow= 0;
   end
   flow=cat(1,flow,Res.flow);
   
   if     strcmp(time_string,'22.3.19'); time_temp= 5;
   elseif strcmp(time_string,'29.3.19'); time_temp= 12;
   elseif strcmp(time_string,'5.4.19'); time_temp= 19;
   elseif strcmp(time_string,'12.4.19'); time_temp= 26;
   end
   time=cat(1,time,time_temp); Res.time=time_temp;
   good=cat(1,good,Res.fov(fov).good);

   %%%% the experiment with 40rpm2 have teh flow direction reverted
   if ~strcmp(flow_string,'40rpm2'); pol= cat(1,pol,(Res.fov(fov).Polx));
   elseif strcmp(flow_string,'40rpm2'); pol= cat(1,pol,-(Res.fov(fov).Polx));
   end
   l_x=cat(1,l_x,Res.fov(fov).fit_out_x.b);
   l_y=cat(1,l_y,Res.fov(fov).fit_out_y.b);
   lf_x=cat(1,lf_x,Res.fov(fov).freq.fit_out_x.b);
   lf_y=cat(1,lf_y,Res.fov(fov).freq.fit_out_y.b);
   
   end
     
end
    good=logical(good);
    Pol=pol(good);
    Time=time(good);
    Flow= flow(good);
    L_x=l_x(good);
    L_y=l_y(good);
    Lf_x=lf_x(good);
    Lf_y=lf_y(good);


    %% plot 
  [G,idtime,idflow]=findgroups(Time,Flow);
  POL=accumarray(G,Pol,[],@median); 
  ePOL=accumarray(G,Pol,[],@std)./sqrt(accumarray(G,Pol,[],@numel)) ;
 
   
  ind_c= (idflow==0);
  ind_exp = ~ind_c; 
  ms=10;lw=2;
 errorbar(idtime(ind_c),POL(ind_c),ePOL(ind_c),'k^','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
 errorbar(idtime(ind_exp),POL(ind_exp),ePOL(ind_exp),'ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','r');
  xlabel('Days','Interpreter','latex') ;ylabel('Polarisation','Interpreter','latex');
  xlim([4,30])
 x0=0;y0=0;width=600;height=500;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',15);
saveas(gcf,'/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/figures/polarisation_time.pdf')

  
%   figure(2);
%   boxplot(Pol(ind1),GTime(ind1),'Colors','r');hold on;
%   boxplot(Pol(ind2),GTime(ind2),'Colors','k');
%     
%     
Lc_x=accumarray(G,L_x,[],@median); 
  eLc_x=accumarray(G,L_x,[],@std)./sqrt(accumarray(G,L_x,[],@numel)) ;
  
  Lc_y=accumarray(G,L_y,[],@median); 
  eLc_y=accumarray(G,L_y,[],@std)./sqrt(accumarray(G,L_y,[],@numel)) ;
  
  Lcf_x=accumarray(G,Lf_x,[],@median); 
  eLcf_x=accumarray(G,Lf_x,[],@std)./sqrt(accumarray(G,Lf_x,[],@numel)) ;
  
  Lcf_y=accumarray(G,Lf_y,[],@median); 
  eLcf_y=accumarray(G,Lf_y,[],@std)./sqrt(accumarray(G,Lf_y,[],@numel)) ;
  

figure() 
errorbar(idtime(ind_c),Lc_y(ind_c),eLc_y(ind_c),'k^','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
 errorbar(idtime(ind_exp),Lc_y(ind_exp),eLc_y(ind_exp),'ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','r');
 xlabel('Days','Interpreter','latex') ;ylabel('Polarisation','Interpreter','latex');
 xlim([4,30]);title('polarisation correlation functions');
 x0=0;y0=0;width=600;height=500;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',15);
saveas(gcf,'/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/figures/corr_function_time.pdf')


%%  plot for polarisation vs shear stress

clear all;
close all
subdir={...   
'29.3.19/beads/40rpm/FL','29.3.19/beads/40rpm2/FL','29.3.19/beads/control/FL'...
'12.4.19/beads/40rpm/FL','12.4.19/beads/40rpm2/FL'...,'12.4.19/beads/control/FL'... 
}; 

holes={...
    [1:2000,4900:10000],'none','none'...% [1:500,500:10000],'none','none'...
    'none',[3865:4075,6415:7195]...
    };

shear_ind=[];
pol=[];
good=[];
flow=[];
lambda=[];
shear_um=1700;
subdir = strcat('/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/',subdir);
for jj=1:numel(subdir)
   cd(subdir{jj}); load('Res.mat');
   for fov=1:numel(Res.fov);
     shear_ind=cat(1,shear_ind,floor(Res.fov(fov).posx/shear_um)+1);

     
     
      allfolders= regexp(Res.insert_name,'/','split'); 
      flow_string=allfolders{end-1};
      time_string=allfolders{end-3};
   
      if strcmp(flow_string,'40rpm') | strcmp(flow_string,'40rpm2'); Res.flow= 40;
      elseif flow_string=='control'; Res.flow= 0;
      end
      flow=cat(1,flow,Res.flow);
   
      if  isa(holes{jj},'char'); if holes{jj}=='none';good=cat(1,good,1);end
      elseif any(Res.fov(fov).posx ==holes{jj}); good=cat(1,good,0); 
      else good=cat(1,good,1); 
      end
   
    %%%% the experiment with 40rpm2 have teh flow direction reverted
      if ~strcmp(flow_string,'40rpm2'); pol= cat(1,pol,(Res.fov(fov).Polx));
      elseif strcmp(flow_string,'40rpm2'); pol= cat(1,pol,-(Res.fov(fov).Polx));
      end
      lambda=cat(1,lambda,Res.fov(fov).fit_out.b);
            

   end
end
    good=logical(good);
    Pol=pol(good);
%    Time=time(good);
    Flow= flow(good);
    Lambda=lambda(good);
    Shear_ind=shear_ind(good);
    %%%
    [G,idshear,idflow]=findgroups(Shear_ind,Flow);
    POL=accumarray(G,Pol,[],@median); 
    ePOL=accumarray(G,Pol,[],@std)./sqrt(accumarray(G,Pol,[],@numel)) ;
    LAMBDA=accumarray(G,Lambda,[],@median); 
    eLAMBDA=accumarray(G,Lambda,[],@std)./sqrt(accumarray(G,Lambda,[],@numel)) ;
  
  ind_c= (idflow==0);
  ind_exp = ~ind_c; 
  ms=10;lw=2;
  shear_c= shear_calculate('tapered',40,idshear(ind_c)*shear_um);
  shear_exp= shear_calculate('tapered',40,idshear(ind_exp)*shear_um);

 errorbar(shear_c,POL(ind_c),ePOL(ind_c),'k^','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
 errorbar(shear_exp,POL(ind_exp),ePOL(ind_exp),'ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','r');
  xlabel('Shear stress [dyne/c$m^2$]','Interpreter','latex') ;ylabel('Polarisation','Interpreter','latex');
  %xlim([0,5])
 x0=0;y0=0;width=600;height=500;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',15);
saveas(gcf,'/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/figures/polarisation_shear.pdf')

%% plot for seeing if polarisation correlate with synchronisation 


%%%%% for plot polarisation of fov
%%% directories with results
clear all
subdir={...   
   '22.3.19/beads/40rpm/FL','22.3.19/beads/40rpm2/FL','22.3.19/beads/control/FL'...
    '29.3.19/beads/40rpm/FL'...
'29.3.19/beads/40rpm2/FL','29.3.19/beads/control/FL'...
    '5.4.19/beads/40rpm/FL','5.4.19/beads/40rpm2/FL','5.4.19/beads/control/FL'...
    '12.4.19/beads/40rpm/FL','12.4.19/beads/40rpm2/FL'...,'12.4.19/beads/control/FL'... 
}; 
%%%%%% hole for seeing the correlation lenght
holes={...
    'none','none','none'...
    [1:2000,4900:10000],'none','none'...
    'none','none','none'...
    'none',[3865:4075,6415:7195]...
    };

subdir = strcat('/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/',subdir);
    
pol=[];good=[];time=[];flow=[];
l_x=[];l_y=[];
lf_x=[];lf_y=[];df=[];
medf=[];
for jj=1:numel(subdir)
   cd(subdir{jj}); load('Res.mat');
   for fov=1:numel(Res.fov);
   
   % determine the experiment and date    
   allfolders= regexp(Res.insert_name,'/','split'); 
   flow_string=allfolders{end-1};
   time_string=allfolders{end-3};
   
   if strcmp(flow_string,'40rpm') | strcmp(flow_string,'40rpm2'); Res.flow= 40;
   elseif flow_string=='control'; Res.flow= 0;
   end
   flow=cat(1,flow,Res.flow);
   
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
   elseif strcmp(flow_string,'40rpm2'); pol= cat(1,pol,-(Res.fov(fov).Polx));
   end
%     if ~strcmp(flow_string,'40rpm2'); pol= cat(1,pol,(Res.fov(fov).Pol));
%     elseif strcmp(flow_string,'40rpm2'); pol= cat(1,pol,(Res.fov(fov).Pol));
%     end


   l_x=cat(1,l_x,Res.fov(fov).fit_out_x.b);
   l_y=cat(1,l_y,Res.fov(fov).fit_out_y.b);
   lf_x=cat(1,lf_x,Res.fov(fov).freq.fit_out_x.b);
   lf_y=cat(1,lf_y,Res.fov(fov).freq.fit_out_y.b);
   % medf= nanmean(Res.fov(fov).freq.F(:));
   % emedf= nanmedian(abs(Res.fov(fov).freq.F(:)-medf));
   F=Res.fov(fov).freq.F;
   df =cat(1,df,nanstd(F(:))/nanmean(F(:)));
   medf= cat(1,medf,nanmedian(F(:)));
  % df =cat(1,df,nanmean(Res.fov(fov).freq.F(:)));
  % df =cat(1,df,emedf);
   
   end
     
end
close all
    good=logical(good);
    Pol=pol(good);
    Time=time(good);
    Flow= flow(good);
    L_x=l_x(good);
    L_y=l_y(good);
    Lf_x=lf_x(good);
    Lf_y=lf_y(good);
    Df=df(good);
    Medf=medf(good);

%%%% plot 
  [G,idtime,idflow]=findgroups(Time,Flow);
  POL=accumarray(G,Pol,[],@median); 
  ePOL=accumarray(G,Pol,[],@std)./sqrt(accumarray(G,Pol,[],@numel)) ;
   Lc_x=accumarray(G,L_x,[],@median); 
  eLc_x=accumarray(G,L_x,[],@std)./sqrt(accumarray(G,L_x,[],@numel)) ;
   Lc_y=accumarray(G,L_y,[],@median); 
  eLc_y=accumarray(G,L_y,[],@std)./sqrt(accumarray(G,L_y,[],@numel)) ;
  Lcf_x=accumarray(G,Lf_x,[],@median); 
  eLcf_x=accumarray(G,Lf_x,[],@std)./sqrt(accumarray(G,Lf_x,[],@numel)) ;
   Lcf_y=accumarray(G,Lf_y,[],@median); 
  eLcf_y=accumarray(G,Lf_y,[],@std)./sqrt(accumarray(G,Lf_y,[],@numel)) ;
   DF=accumarray(G,Df,[],@median); 
  eDF=accumarray(G,Df,[],@std)./sqrt(accumarray(G,Df,[],@numel)) ;

  
  ms=10;lw=2;

ind_c= (idflow==0);
 ind_exp = ~ind_c; 
  ms=10;lw=2;
 errorbar(idtime(ind_c),Lcf_y(ind_c),DF(ind_c),'k^','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
 errorbar(idtime(ind_exp),Lcf_y(ind_exp),DF(ind_exp),'ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','r');
  xlabel('Days','Interpreter','latex') ;ylabel('Polarisation','Interpreter','latex');
 % xlim([4,30])
 x0=0;y0=0;width=600;height=500;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',15);
saveas(gcf,'/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/figures/polarisation_time.pdf')

    
figure(5) 
L=mean([L_y,L_x],2);
Lf=mean([Lf_y,Lf_x],2);
x= Lf_x(Flow==40);y=Lf_y(Flow==40);n_bins=5;
plot(x,y,'k^','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
[histw,ehistw,m_bins]=hist_nico2(x,y,[],[20:20:100]);
errorbar(m_bins,histw,ehistw,'-ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
%plot(L(Flow==0),Df(Flow==0),'ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
 xlabel('$\zeta_{or}$','Interpreter','latex') ;ylabel('$\zeta_{f}$','Interpreter','latex');
 title('correlation orientation lenght vs frequency length');
 x0=0;y0=0;width=600;height=500;
set(gcf,'position',[x0,y0,width,height]);
set(gca,'FontSize',15);
saveas(gcf,'/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/figures/corr_orientation_vs_freq.pdf')

figure(4)

x= Pol(Flow==40);y=Df(Flow==40);n_bins=5;
plot(x,y,'k^','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
[histw,ehistw,vinterval]=hist_nico(x,y,[0.6:0.1:1]);
errorbar(vinterval,histw,ehistw,'-ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
plot(x,y,'k^','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
%plot(Pol(Flow==0),Df(Flow==0),'ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
 xlabel('$\zeta_{or}$','Interpreter','latex') ;ylabel('$\sigma_{f}$[Hz]','Interpreter','latex');
 title('Polarisation vs frequency detuning');
 x0=0;y0=0;width=600;height=500;
set(gcf,'position',[x0,y0,width,height]);
set(gca,'FontSize',15);
saveas(gcf,'/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/figures/polarisation_vs_freq_detuning.pdf')

figure();
plot(Df(Flow==40),Medf(Flow==40),'k^','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;

%% general staff, recalculate correlation functions

plot()






