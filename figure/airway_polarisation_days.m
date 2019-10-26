%%  plot for polarisation vs shear stress

clear all;
close all;


subdir={...   
    'half_diff/26.3.19/10rpm/FL'...
...    'half_diff/26.3.19/10rpm2/FL'...
   'half_diff/26.3.19/control/FL'...
...    'half_diff/27.3.19/10rpm/FL','half_diff/27.3.19/10rpm2/FL','half_diff/27.3.19/control/FL'...
...   'end_diff/10.4.19/beads/40rpm/FL'...
   'end_diff/10.4.19/beads/40rpm2/FL'...
   'end_diff/10.4.19/beads/control/FL'...
}; subdir=strcat('/media/np451/Seagate Backup Plus Drive1/DATA/airway_paper/fluid/PIV/',subdir)


holes={...
 'none','none','none'...    
 'none','none','none'... 
  'none','none','none'...   
  };

shear_ind=[];
pol=[];
good=[];
flow=[];
npo=[]
shear_um=1700;
cc=1;
for jj=1:numel(subdir)
   cd(subdir{jj}); load('Res3.mat');
   for fov=1:numel(Res.fov);

        pol_sign=1;   %%% sign that accounts for the direction of the flow respective to camera 
      allfolders= regexp(Res.insert_name,'/','split'); 
      flow_string=allfolders{end-1};
%      if strcmp(allfolders{end-2},'beads'); time_string=allfolders{end-3};
%      else time_string=allfoldersend{end-2};
%      end

      if strcmp(flow_string,'40rpm') | strcmp(flow_string,'40rpm2'); Res.flow= 40;
                shear_temp= shear_calculate('straight',Res.flow);
      elseif strcmp(flow_string,'10rpm') | strcmp(flow_string,'10rpm2'); Res.flow= 10;
                shear_temp= shear_calculate('step',Res.flow,Res.fov(fov).posx);
      elseif strcmp(flow_string,'control') | strcmp(flow_string,'0'); Res.flow= 0;
                shear_temp= shear_calculate('step',Res.flow,Res.fov(fov).posx);
      end
 

        shear_ind= cat(1,shear_ind,shear_temp);
      
      flow=cat(1,flow,Res.flow);

      nn= numel(Res.fov(fov).nu)./numel(Res.fov(fov).x);
      npo=cat(1,npo,nn);
   
      if  isa(holes{jj},'char'); if holes{jj}=='none';good=cat(1,good,1);end
      elseif any(Res.fov(fov).posx ==holes{jj}); good=cat(1,good,0); 
      else good=cat(1,good,1); 
      end
   
      if isnan(shear_ind(end)); good(end)=0;end
    %%%% the experiment with 40rpm2 have then flow direction reverted

        if strcmp(flow_string,'40rpm2') | strcmp(flow_string,'10rpm2');pol_sign=-1;
        end
        disp(pol_sign)
        if isempty(Res.fov(fov).Polx);
            pol= cat(1,pol,nan);good(end)=0;
        else pol= cat(1,pol,pol_sign*Res.fov(fov).Polx);
        end
   end
end
    good=logical(good);
    Pol=pol(good);
    Npo=npo(good);
%    Time=time(good);
    Flow= flow(good);
    Shear_ind=shear_ind(good);
    Flow(end-12:end)=-2 %%%% it is to distinguish the first from the second control
    %%%
    [G,idshear,idflow]=findgroups(Shear_ind,Flow);
    POL=accumarray(G,Pol,[],@median); 
    ePOL=accumarray(G,Pol,[],@std)./sqrt(accumarray(G,Pol,[],@numel)) ;
    NPO=accumarray(G,Npo,[],@median); 
    eNPO=accumarray(G,Npo,[],@std)./sqrt(accumarray(G,Npo,[],@numel)) ;

    
  ind_c_half= (idflow==0);
ind_c_end= (idflow==-2);

  ind_end = (idflow==40); 
  ind_half = (idflow==10); 
  ms=10;lw=2;
  
%  shear_end= idshear(ind_end);
%shear_c= idshear(ind_c);
%shear_half= idshear(ind_half);

%%
  %%%% polarisation vs shear stress
 ms=8 ; lw=1
 errorbar(15,mean(POL(ind_c_half)),mean(ePOL(ind_c_half)),'k^','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
errorbar(15,mean(POL(ind_half)),mean(ePOL(ind_half)),'ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','r');hold on;
errorbar(35,mean(POL(ind_c_end)),mean(ePOL(ind_c_end)),'k>','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
errorbar(35,mean(POL(ind_end)),mean(ePOL(ind_end))*1.1,'bd','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','b');hold on;

% xlabel('Shear stress [dyne/c$m^2$]','Interpreter','latex') ;ylabel('Polarisation','Interpreter','latex');
 %ylim([0.5,1]); %xlim([0,5])
 x0=0;y0=0;width=300;height=300;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',13);
legend({'control','0.2 dyne/cm$^2$','control','0.2 dyne/cm$^2$'},'Interpreter','latex')
%legend({'control','14 DIV'},'Interpreter','latex')
saveas(gcf,'/home/np451/Dropbox/phd_thesis/Chapter5/Figs/polarisation_days.pdf')


