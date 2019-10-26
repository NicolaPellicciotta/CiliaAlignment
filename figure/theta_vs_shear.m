%%  plot for variation of direction vs shear
subdirMarch={...   
...'29.3.19/beads/40rpm/FL',
...'29.3.19/beads/40rpm2/FL'...,'29.3.19/beads/control/FL'...
'12.4.19/beads/40rpm/FL'...,'12.4.19/beads/40rpm2/FL'...,'12.4.19/beads/control/FL'... 
}; 

subdirMarch = strcat('/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/',subdirMarch);

subdir=[subdirMarch]

holes={...
'none'
    };

shear_ind=[];shear_vec=[];
pol=[];
good=[];
flow=[];
npo=[];
nu=[];nv=[];
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


       shear_temp= shear_calculate('tapered',40,shear_um*(floor(Res.fov(fov).posx/shear_um)+1));

       shear_ind=cat(1,shear_ind,shear_temp);
       nu=cat(1,nu,Res.fov(fov).nu(:));
       nv=cat(1,nv,Res.fov(fov).nv(:));
       dummy=shear_temp*ones(size(Res.fov(fov).nu));
       shear_vec= cat(1,shear_vec, dummy(:));
       
      if strcmp(flow_string,'40rpm') | strcmp(flow_string,'40rpm2'); Res.flow= 40;
      elseif strcmp(flow_string,'control') | strcmp(flow_string,'0'); Res.flow= 0;
          
      end
      
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

%%
% 
%     good=logical(good);
%     Pol=pol(good);
%     Npo=npo(good);
% %    Time=time(good);
%     Flow= flow(good);
%     Shear_ind=shear_ind(good);
%     %%%

    [G,idshear]=findgroups(shear_vec);
    du=nu-nanmean(nu),dv=nv-nanmean(nv);cc=1;
    theta=angle(complex(nu,nv))-angle(complex(nanmean(nu(:)),nanmean(nv(:)))); 
    colors={'v','red','blue','green','y'};
for ii=flip(2:5);
    x= nu(G==ii); y=dv(G==ii);shear=shear_vec(G==ii);the=theta(G==ii);
    [N,edges]=histcounts(the,[-pi:0.1:pi],'Normalization','probability');
    edge=(edges(1:end-1)+edges(2:end))/2
    %plot(edge,smooth(N),'-','LineWidth',2);hold on;
    patch(edge,smooth(N),colors{ii},'FaceAlpha',0.6);hold on;
    cleg{cc}=strcat(num2str(shear(1),2) ,' dyne/cm$^2$');cc=cc+1;
end
xlabel('dv');
ylabel('probability');
%xlim([-1,0]);
legend(cleg,'Interpreter','latex');
 x0=0;y0=0;width=600;height=500;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',18);
saveas(gcf,'/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/figures/probailitytheta_shear_12.4.19.40rpm.pdf')
