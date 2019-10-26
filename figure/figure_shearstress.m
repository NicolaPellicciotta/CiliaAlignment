%%  plot for polarisation vs shear stress

clear all;
close all;

subdirJuly={...   
   '16.7.18/0/FL'...
   '29.7.18/0/FL'...
   '16.7.18/5A/FL', '16.7.18/5B/FL', '16.7.18/5C/FL'... %%% straight
  }; 
subdirJuly = strcat('/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/July2018/PIV/',subdirJuly);
 
subdirJune={...   
   '13.6.18/0/FL'...
   '18.6.18/0/FL'...
...   '18.6.18/100/FL',
   '18.6.18/200/FL'...  %% straight
   }; 
subdirJune = strcat('/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/June2018/PIV/',subdirJune);

subdirMarch={...   
'29.3.19/beads/40rpm/FL','29.3.19/beads/40rpm2/FL','29.3.19/beads/control/FL'...
'12.4.19/beads/40rpm/FL','12.4.19/beads/40rpm2/FL'...,'12.4.19/beads/control/FL'... 
'5.4.19/beads/40rpm2/FL','5.4.19/beads/40rpm/FL'
}; 

subdirMarch = strcat('/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/',subdirMarch);

subdir=[subdirMarch,subdirJuly,subdirJune]

holes={...
    [1:2000,4900:10000],'none','none'...% [1:500,500:10000],'none','none'...
    'none',[3865:4075,6415:7195],...
    'none','none','none','none',...
 'none','none','none'...
 'none','none','none','none'
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

     pol_sign=1;   %%% sign that accounts for the direction of the flow respective to camera 
      allfolders= regexp(Res.insert_name,'/','split'); 
      flow_string=allfolders{end-1};
      if strcmp(allfolders{end-2},'beads'); time_string=allfolders{end-3};
      else time_string=allfolders{end-2};
      end

      if strcmp(flow_string,'40rpm') | strcmp(flow_string,'40rpm2'); Res.flow= 40;
      elseif strcmp(flow_string,'control') | strcmp(flow_string,'0'); Res.flow= 0;
          
      end
      
      if strcmp(time_string,'18.6.18')  & strcmp(flow_string,'0')
          shear_temp= shear_calculate('tapered',40,shear_um*(floor(1*(shear_um-1)/shear_um)+1));
      elseif strcmp(time_string,'13.6.18')  & strcmp(flow_string,'0')
          shear_temp= shear_calculate('tapered',40,shear_um*(floor(2*(shear_um-1)/shear_um)+1));
      elseif strcmp(time_string,'16.7.18')  & strcmp(flow_string,'0')
         shear_temp= shear_calculate('tapered',40,shear_um*(floor(3*(shear_um-1)/shear_um)+1));
      elseif strcmp(time_string,'29.7.18')  & strcmp(flow_string,'0')
        shear_temp= shear_calculate('tapered',40,shear_um*(floor(4*(shear_um-1)/shear_um)+1));
      elseif strcmp(flow_string,'100') | strcmp(flow_string,'200') 
         shear_temp= shear_calculate('straight',10);  Res.flow= 10; pol_sign=-1;
      elseif strcmp(flow_string,'5A') | strcmp(flow_string,'5B') | strcmp(flow_string,'5C')
         shear_temp= shear_calculate('straight',5);  Res.flow= 5;
      elseif strcmp(time_string,'5.4.19') 
         shear_temp= shear_calculate('straight',40);  Res.flow= 15;
      else
         shear_temp= shear_calculate('tapered',40,shear_um*(floor(Res.fov(fov).posx/shear_um)));
      end

       shear_ind=cat(1,shear_ind,shear_temp);
       
       

      flow=cat(1,flow,Res.flow);
      nn= numel(Res.fov(fov).nu)./numel(Res.fov(fov).x);
      npo=cat(1,npo,nn);
   
      if  isa(holes{jj},'char'); if holes{jj}=='none';good=cat(1,good,1);end
      elseif any(Res.fov(fov).posx ==holes{jj}); good=cat(1,good,0); 
      else good=cat(1,good,1); 
      end
   
    %%%% the experiment with 40rpm2 have teh flow direction reverted
      if ~strcmp(flow_string,'40rpm2'); 
      elseif strcmp(flow_string,'40rpm2');pol_sign=-1;;
      end
      
      pol= cat(1,pol,pol_sign*(Res.fov(fov).Polx));
 
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
  ind_tap = (idflow==40); 
  ind_str = (idflow~=40) & (idflow~=0); 
  ms=10;lw=2;
  
  shear_tap= idshear(ind_tap);
shear_c= idshear(ind_c);
shear_str= idshear(ind_str);

%%
  %%%% polarisation vs shear stress
 ms=8 ; lw=1
% errorbar(shear_c,POL(ind_c),ePOL(ind_c),'k^','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
errorbar(shear_str,POL(ind_str),ePOL(ind_str),'ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','r');hold on;
errorbar(shear_tap,POL(ind_tap),ePOL(ind_tap),'bd','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','b');hold on;
% xlabel('Shear stress [dyne/c$m^2$]','Interpreter','latex') ;ylabel('Polarisation','Interpreter','latex');
 ylim([0.5,1]); %xlim([0,5])
 x0=0;y0=0;width=300;height=300;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',13);
legend({'straight','tapered',},'Interpreter','latex')
saveas(gcf,'/u/homes/np451/Dropbox/alignment with the flow/figure2/polarisation_shear_lintap.pdf')
saveas(gcf,'/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/figures/polarisation_shear_lintap.pdf')
%%
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

%%



