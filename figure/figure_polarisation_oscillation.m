%%%%% for plot polarisation of fov
%%% directories with results
clear all

subdirMarch={...   
   '29.3.19/beads/control/FL'...
 }; 
subdirMarch = strcat('/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/',subdirMarch);

subdirOsc={...  
    'sample1/FL'...
}; 
subdirOsc = strcat('/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/oscillation_longtime/March2019/27.3.19/beads/PIV/',subdirOsc);

subdir=[subdirMarch,subdirOsc];

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

   
     if strcmp(flow_string,'control') | strcmp(flow_string,'0'); Res.flow= 0;
    else Res.flow=40;
   end
   flow=cat(1,flow,Res.flow);
   
   
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
    Flow= flow(good);
    L=l(good);
    Lf=lf(good);
    N=n(good).*(32^2)/(1200*1920);

    %% plot 
  [G,idflow]=findgroups(Flow);
  POL=accumarray(G,Pol,[],@median); 
  ePOL=accumarray(G,Pol,[],@std)./sqrt(accumarray(G,Pol,[],@numel)) ;
  NN=accumarray(G,N,[],@median); 
  eNN=accumarray(G,N,[],@std)./sqrt(accumarray(G,N,[],@numel)) ;
  
 
   
  ind_c= (idflow==0);
  ind_exp = ~ind_c; 
  ms=8;lw=1;
  eta= 0.8*10^-3 %% Pas
  flow_z= 2.5*10^-3
  z=7*10^-6;
  shear = eta * (flow_z)/z;
 errorbar(15,POL(ind_c),ePOL(ind_c),'k^','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
 errorbar(15,POL(ind_exp),ePOL(ind_exp),'ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','r');
 % xlabel('Days','Interpreter','latex') ;ylabel('Polarisation','Interpreter','latex');
  xlim([4,30])
  ylim([-0.3,0.3]);
 x0=0;y0=0;width=350;height=300;
 
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',13);
xlabel('Days in vitro DIV','Interpreter','latex');
ylabel('Polarisation $\Phi$','Interpreter','latex');
legend({'control',strcat(num2str(shear,1), ' dyne/cm$^2$ at 10Hz')},'Interpreter','latex')
saveas(gcf,'/u/homes/np451/Dropbox/alignment with the flow/figure_sup/oscillation_label.pdf')
