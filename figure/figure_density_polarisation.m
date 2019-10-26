%% general staff, calculate staff in the insert

%%%%% for plot polarisation of fov
%%% directories with results
clear all
subdir={...
...          '22.3.19/beads/40rpm/FL'...
       '22.3.19/beads/40rpm2/FL'...
 ...        '22.3.19/beads/control/FL'...
...              '29.3.19/beads/40rpm/FL'...
           '29.3.19/beads/40rpm2/FL'...
 ...           '29.3.19/beads/control/FL'...
  ...               '5.4.19/beads/40rpm/FL'...
                '5.4.19/beads/40rpm2/FL'...
...                '5.4.19/beads/control/FL'...
          '12.4.19/beads/40rpm/FL'...
         '12.4.19/beads/40rpm2/FL'...,
        }; 
%%%%%% hole for seeing the correlation lenght
holes={...
   'none'...
...  'none'...
    'none'...
...   [1:2000,4900:10000]...
   'none'...
...   'none'...
...   'none'...
    'none'...
 'none'...
   'none'...
 [3865:4075,6415:7195]...
    };

subdir = strcat('/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/',subdir);
    
pol=[];good=[];time=[];flow=[];
l=[];
lf=[];df=[];nn=[];pp=[];
medf=[];nf=[];m=[];mv=[];f=[];ddf=[];uv=[];
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
 
   if     strcmp(time_string,'22.3.19'); time_temp= 5;
   elseif strcmp(time_string,'29.3.19'); time_temp= 12;
   elseif strcmp(time_string,'5.4.19'); time_temp= 19;
   elseif strcmp(time_string,'12.4.19'); time_temp= 26;
   end

       if  isa(holes{jj},'char'); if holes{jj}=='none';good=cat(1,good,1);end
      elseif any(Res.fov(fov).posx ==holes{jj}); good=cat(1,good,0); 
      else good=cat(1,good,1); 
      end
   

   %%%% the experiment with 40rpm2 have teh flow direction reverted
 
   if ~strcmp(flow_string,'40rpm2');
        flow_sign=1;
   
   elseif strcmp(flow_string,'40rpm2'); 
        flow_sign=-1;
   end
   
     pol= cat(1,pol,flow_sign*(Res.fov(fov).Polx));
     pp= cat(1,pp,flow_sign*Res.fov(fov).nv(:)); 
   
   
   if isempty(Res.fov(fov).nu) | numel(Res.fov(fov).nu)<20
        l=cat(1,l,nan);
        lf=cat(1,lf,nan);
        good(end)=0;
    else   
        l=cat(1,l,Res.fov(fov).fit_out.b);
        lf=cat(1,lf,Res.fov(fov).freq.fit_out.b);
    end
   ind= Res.fov(fov).ind;
   F=Res.fov(fov).freq.F32(ind);
   f=cat(1,f,F(:));

%   m_temp=sqrt(Res.fov(fov).u1m(ind).^2+ Res.fov(fov).v1m(ind).^2);
   %%%%% find frame rate
if fov==1
    cd(Res.fov(fov).file_mat(1:end- numel(Res.fov(fov).filename)))
   fl_var=load(Res.fov(fov).filename); fps=fl_var.mo.FrameRate;
end
   mv_temp=sqrt(Res.fov(fov).v1m(ind).^2)
   mv_temp=mv_temp*fps;
   mv=cat(1,mv, mv_temp(:));
      
   uv_temp=Res.fov(fov).v1m(ind)
   uv_temp=flow_sign*uv_temp*fps;
   uv=cat(1,uv, uv_temp(:));
   
   m_temp=sqrt(Res.fov(fov).v1m(ind).^2+Res.fov(fov).u1m(ind).^2)
   m_temp=m_temp*fps;
   m=cat(1,m, m_temp(:));

   df =cat(1,df,(nanstd(F(:))/nanmean(F(:))));
   medf= cat(1,medf,nanmedian(F(:)));
   nf= cat(1,nf,(numel(F(~isnan(F(:))))));
    d_lim=180;
   fXX=Res.fov(fov).x(Res.fov(fov).ind);
   fYY=Res.fov(fov).y(Res.fov(fov).ind);
   F32= Res.fov(fov).freq.F32(Res.fov(fov).ind);
   
   time_temp=ones([numel(fXX),1])*time_temp;
   time=cat(1,time,time_temp);

   near= find_nearcilia(fXX(:),fYY(:),d_lim);
   [near2,df] = find_nearcilia_df(fXX(:),fYY(:),F32,d_lim)
   nn= cat(1,nn,near(:));
   ddf=cat(1,ddf,df(:));
%   ciao=nan(size(Res.fov(fov).x));
%   ciao(Res.fov(fov).ind)=m_temp(:);
%    
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

 rho=nn/(pi*(d_lim/32)^2);   
    
%% density and polarisation    
 %%%% now it is a density of cell
figure();ms=3;lw=0.5;
x=rho;y=pp;n_bins=5;
plot(x,y,'k^','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
[histw,ehistw,m_bins]=hist_nico2(x,y,[],[0:0.05:0.4]);
errorbar(m_bins,histw,ehistw,'-ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
%plot(L(Flow==0),Df(Flow==0),'ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
 xlabel('$\rho $','Interpreter','latex') ;ylabel('$\Phi$','Interpreter','latex');
 title('local density vs polarisation');
 x0=0;y0=0;width=600;height=500;
set(gcf,'position',[x0,y0,width,height]);
set(gca,'FontSize',15);
%% density and frequency

rho=nn/(pi*(d_lim/32)^2); %%%% now it is a density of cell
figure();ms=10;lw=2;
x=rho;y=f;n_bins=5;
%plot(x,y,'k^','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
[histw,ehistw,m_bins]=hist_nico2(x,y,[],[0:0.05:0.4]);
errorbar(m_bins,histw,ehistw,'-bo','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','b');hold on;
%plot(L(Flow==0),Df(Flow==0),'ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
 xlabel('$\rho $','Interpreter','latex') ;ylabel('$CBF$','Interpreter','latex');
% title('local density vs polarisation');
 x0=0;y0=0;width=600;height=500;
set(gcf,'position',[x0,y0,width,height]);
set(gca,'FontSize',15);
%% polarisation vs density

figure();ms=10;lw=2;
x=pp;y=mm;n_bins=5;
%plot(x,y,'ko','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
%[histw,ehistw,m_bins]=hist_nico2(x,y,[],[1:7:40]);
[histw,ehistw,m_bins]=hist_nico2(x,y,[],[0.5:0.1:1])
errorbar(m_bins,histw,ehistw,'-ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','r');hold on;
%plot(L(Flow==0),Df(Flow==0),'ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
 xlabel('$\Phi$','Interpreter','latex') ;ylabel('flow $\mu$m/s','Interpreter','latex');
 %title('correlation orientation lenght vs frequency length');
 x0=0;y0=0;width=600;height=500;
set(gcf,'position',[x0,y0,width,height]);
set(gca,'FontSize',15);

%% density and flow magnitude
figure();ms=10;lw=2;
x=rho(time==5);y=m(time==5);n_bins=5;
%plot(x,y,'ko','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
%[histw,ehistw,m_bins]=hist_nico2(x,y,[],[1:7:40]);
[histw,ehistw,m_bins]=hist_nico2(x,y,[],[0:0.05:0.4])
errorbar(m_bins,histw,ehistw,'-ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','r');hold on;
%plot(L(Flow==0),Df(Flow==0),'ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
 xlabel('$\rho$','Interpreter','latex') ;ylabel('flow $\mu$m/s','Interpreter','latex');
% title('correlation orientation lenght vs frequency length');
 x0=0;y0=0;width=600;height=500;
set(gcf,'position',[x0,y0,width,height]);
set(gca,'FontSize',18);

x=rho(time==26);y=m(time==26);n_bins=5;
%plot(x,y,'ko','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
%[histw,ehistw,m_bins]=hist_nico2(x,y,[],[1:7:40]);
[histw,ehistw,m_bins]=hist_nico2(x,y,[],[0:0.05:0.4])
errorbar(m_bins,histw,ehistw,'-bd','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','b');hold on;
%plot(L(Flow==0),Df(Flow==0),'ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
 xlabel('$\rho$','Interpreter','latex') ;ylabel('flow $\mu$m/s','Interpreter','latex');
% title('correlation orientation lenght vs frequency length');
 x0=0;y0=0;width=600;height=500;
set(gcf,'position',[x0,y0,width,height]);
set(gca,'FontSize',18);
legend({'DIV 8','DIV 29'},'Interpreter','latex')

%% density and flow magnitude over time
figure();ms=10;lw=2;

divs=unique(time);colors={'red','green','blue','yellow'};cc=1;
figure()
for jj=1:numel(divs)
x=rho(time==divs(jj));y=m(time==divs(jj));n_bins=5;
[histw,ehistw,m_bins]=hist_nico2(x,y,[],[0:0.05:0.4]);
errorbar(m_bins,histw,ehistw,strcat(colors{jj},'-o'),'MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','r');hold on;
cleg{cc}=strcat(num2str(divs(jj),2) ,' DIV');cc=cc+1;
end

xlabel('$\rho$','Interpreter','latex') ;ylabel('flow $\mu$m/s','Interpreter','latex');
 x0=0;y0=0;width=600;height=500;
set(gcf,'position',[x0,y0,width,height]);
set(gca,'FontSize',18);

legend(cleg,'Interpreter','latex')


%%
figure();ms=3;lw=0.5;
x=rho;y=ddf;n_bins=5;
plot(x,y,'ko','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
[histw,ehistw,m_bins]=hist_nico2(x,y,[],[0:0.05:0.4]);
errorbar(m_bins,histw,ehistw,'-ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
%plot(L(Flow==0),Df(Flow==0),'ro','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','k');hold on;
 xlabel('$\rho$','Interpreter','latex') ;ylabel('$\Delta f/\overline{f}$','Interpreter','latex');
 %title('correlation orientation lenght vs frequency length');
 x0=0;y0=0;width=600;height=500;
set(gcf,'position',[x0,y0,width,height]);
set(gca,'FontSize',15);
