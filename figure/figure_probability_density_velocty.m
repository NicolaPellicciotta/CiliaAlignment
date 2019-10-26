%%%%% for plot polarisation of fov
%%% directories with results
clear all
subdirMarch={...   
  ...'22.3.19/beads/40rpm/FL',
  '22.3.19/beads/40rpm2/FL'...
  '22.3.19/beads/control/FL'...
...    '29.3.19/beads/40rpm/FL',
'29.3.19/beads/40rpm2/FL','29.3.19/beads/control/FL'...
  '5.4.19/beads/40rpm/FL'...
'5.4.19/beads/40rpm2/FL'...
...    '5.4.19/beads/control/FL'...
    '12.4.19/beads/40rpm/FL'...
    '12.4.19/beads/40rpm2/FL'...,'12.4.19/beads/control/FL'... 
}; 
subdirMarch = strcat('/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/',subdirMarch);

subdirJuly={...   
...  '12.7.18/5A/FL', '12.7.18/5B/FL', '12.7.18/5C/FL'...
...   '16.7.18/5A/FL', '16.7.18/5B/FL', '16.7.18/5C/FL'...
...   '16.7.18/0/FL'...
...  '12.7.18/15A/FL', '12.7.18/15B/FL', '12.7.18/15C/FL'...
 ... '29.7.18/0/FL'...
...   , '29.7.18/5A/FL', '29.7.18/5B/FL'...
  }; 
subdirJuly = strcat('/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/July2018/PIV/',subdirJuly);
 
subdirJune={...   
...   '13.6.18/30A/FL', '13.6.18/30B/FL'...
...  '13.6.18/0/FL'...
...   '13.6.18/100A/FL', '13.6.18/1000/FL'...
...  '18.6.18/0/FL'...
...   , '18.6.18/100/FL', '18.6.18/200/FL'...

   }; 
subdirJune = strcat('/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/June2018/PIV/',subdirJune);

subdir=[subdirMarch,subdirJuly,subdirJune]

DIV={5,5,12,12,19,19,26,26}...,5,12,5,12}
...DIV={19,19,26,5,12};


time=[];nn=[];den=[];mm=[];
medf=[];nf=[];m=[];mv=[];f=[];mm=[];;ddf=[];
for jj=1:numel(subdir)
   cd(subdir{jj}); load('Res3.mat');
   for fov=1:numel(Res.fov);
   
   % determine the experiment and date    
   allfolders= regexp(Res.insert_name,'/','split'); 
   flow_string=allfolders{end-1};
   time_string=allfolders{end-3};
   
   
%   good=cat(1,good,Res.fov(fov).good);

   d_lim=180;
   fXX=Res.fov(fov).x(Res.fov(fov).ind);
   fYY=Res.fov(fov).y(Res.fov(fov).ind);
    time_temp=ones([numel(fXX),1])*DIV{jj};
   time=cat(1,time,time_temp); 
   f_temp= Res.fov(fov).freq.F32(Res.fov(fov).ind);
   f=cat(1,f,f_temp(:));
   
   near= find_nearcilia(fXX(:),fYY(:),d_lim);
   nn= cat(1,nn,near(:));
    ind=Res.fov(fov).ind;
   m_temp=sqrt(Res.fov(fov).u1m(ind).^2+ Res.fov(fov).v1m(ind).^2);

   if fov==1
        cd(Res.fov(fov).file_mat(1:end- numel(Res.fov(fov).filename)))
        fl_var=load(Res.fov(fov).filename); fps=fl_var.mo.FrameRate;
   end
   m_temp=m_temp*fps;
%  m_temp=sqrt(Res.fov(fov).u1m(ind).^2)
   mm=cat(1,mm, m_temp(:));
   
   end
     
end
close all
 
 rho=nn/(pi*(d_lim/32)^2); 
 time=time+3;
    
%% density ditribution with DIV
divs=unique(time);colors={'red','green','blue','yellow'};cc=1;
figure()
for jj=1:numel(divs)
x=rho(time==divs(jj)); mrho(jj)=median(x);
[N,edges]=histcounts(x,[0:0.01:0.5],'Normalization','probability');
    edge=(edges(1:end-1)+edges(2:end))/2
    %plot(edge,smooth(N),'-','LineWidth',2);hold on;
    patch(edge,smooth(N),colors{jj},'FaceAlpha',0.6);hold on;
    cleg{cc}=strcat(num2str(divs(jj),2) ,' DIV');cc=cc+1;
end


for jj=1:numel(divs)
    vy=0:0.05:0.1;vx=mrho(jj).*ones([numel(vy),1]);hold on;
    plot(vx,vy,colors{jj},'LineWidth',2)
end
legend(cleg,'Interpreter','latex');

% xlabel('$\rho$','Interpreter','latex') ;ylabel('Probability','Interpreter','latex');
% title('correlation orientation lenght vs frequency length');
 x0=0;y0=0;width=350;height=300;
 xlim([0,0.4]);
set(gcf,'position',[x0,y0,width,height]);
set(gca,'FontSize',13);
saveas(gcf,'/u/homes/np451/Dropbox/alignment with the flow/figure4/density_distribution_time.pdf')


%%
figure();
boxplot(rho,time);
 x0=0;y0=0;width=350;height=300;
ylabel('$\rho$','Interpreter','latex')
set(gcf,'position',[x0,y0,width,height]);
set(gca,'FontSize',13);
saveas(gcf,'/u/homes/np451/Dropbox/alignment with the flow/figure4/median_density_time.pdf')

%% density ditribution with DIV
divs=unique(time);colors={'red','green','blue','yellow'};cc=1;
figure()
for jj=1:numel(divs)
x=f(time==divs(jj)); mf(jj)=median(x);
[N,edges]=histcounts(x,[14:1:32],'Normalization','probability');
    edge=(edges(1:end-1)+edges(2:end))/2
    %plot(edge,smooth(N),'-','LineWidth',2);hold on;
   area(edge,smooth(N),'FaceColor',colors{jj},'FaceAlpha',0.6);hold on;
    cleg{cc}=strcat(num2str(divs(jj),2) ,' DIV');cc=cc+1;
end


for jj=1:numel(divs)
    vy=0:0.05:0.15;vx=mf(jj).*ones([numel(vy),1]);hold on;
    plot(vx,vy,colors{jj},'LineWidth',2)
end
legend(cleg,'Interpreter','latex');

 xlabel('CBF [Hz]','Interpreter','latex') ;ylabel('Probability','Interpreter','latex');
% title('correlation orientation lenght vs frequency length');
 x0=0;y0=0;width=600;height=500;
 xlim([14,32]);
set(gcf,'position',[x0,y0,width,height]);
set(gca,'FontSize',18);


figure();
boxplot(f,time);
 x0=0;y0=0;width=600;height=500;
ylabel('$CBF [Hz]$','Interpreter','latex')
set(gcf,'position',[x0,y0,width,height]);
set(gca,'FontSize',18);



%% velocity ditribution with DIV
divs=unique(time);colors={'red','green','blue','yellow'};cc=1;
figure()
for jj=1:numel(divs)
x=mm(time==divs(jj)); mrho(jj)=median(x);
[N,edges]=histcounts(x,[0:5:120],'Normalization','probability');
    edge=(edges(1:end-1)+edges(2:end))/2
   % plot(edge,smooth(N),strcat(colors{jj},'-'),'LineWidth',2);hold on;
    area(edge,smooth(N),'FaceColor',colors{jj},'FaceAlpha',0.6);hold on;
   % patch(edge,smooth(N),colors{jj},'FaceAlpha',0.6);hold on;
    cleg{cc}=strcat(num2str(divs(jj),2) ,' DIV');cc=cc+1;
end



for jj=1:numel(divs)
    vy=0:0.05:0.15;vx=mrho(jj).*ones([numel(vy),1]);hold on;
    plot(vx,vy,colors{jj},'LineWidth',2)
end
legend(cleg,'Interpreter','latex');

 xlabel('flow magnitude $\mu$m/s','Interpreter','latex') ;ylabel('Probability','Interpreter','latex');
% title('correlation orientation lenght vs frequency length');
 x0=0;y0=0;width=600;height=500;
 xlim([0,120]);
set(gcf,'position',[x0,y0,width,height]);
set(gca,'FontSize',18);

