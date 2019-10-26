%% figure all flow 
clear all
%%%%%%%%%% Load data to find mean polarisation over multiple field of view

%%%%%% this script is suitable with PIV_script_maker
%%%% it read in order of date
ppm= 0.14*2;  %%% 20X
%path_dir = '/home/np451/Desktop/ependymal/June2018/July/16.7.18/PIV/';
path_dir = '/media/np451/Seagate Backup Plus Drive1/DATA/airway_paper/fluid'
cd(path_dir); 

subdir={...   
...    'half_diff/26.3.19/10rpm/FL','half_diff/26.3.19/10rpm2/FL','half_diff/26.3.19/control/FL'...
...    'half_diff/27.3.19/10rpm/FL','half_diff/27.3.19/10rpm2/FL','half_diff/27.3.19/control/FL'...
    'end_diff/10.4.19/beads/40rpm/FL','end_diff/10.4.19/beads/40rpm2/FL','end_diff/10.4.19/beads/control/FL'...
}; subdir=strcat('PIV/',subdir)


posrange={...
    'all','all','all'...
   'all','all','all'...
    'all','all','all'...
    };
allfov={...
    0,0,0,...
    0,0,0,...
    0,0,0,...
    };

BFdir={...   
...    'half_diff/26.3.19/10rpm/BF','half_diff/26.3.19/10rpm2/BF','half_diff/26.3.19/control/BF'...
...    'half_diff/27.3.19/10rpm/BF','half_diff/27.3.19/10rpm2/BF','half_diff/27.3.19/control/BF'...
    'end_diff/10.4.19/beads/40rpm/BF','end_diff/10.4.19/beads/40rpm2/BF','end_diff/10.4.19/beads/control/BF'...
}; %BFdir=strcat(BFdir)


%%
for jj=1:numel(subdir)
    
start=0;
cd(path_dir);cd(subdir{jj}); d=dir('*FL*.mat');

clear xx;clear yy;clear uu;clear vv;clear uu1;clear vv1;
clear XX;clear YY;clear UU;clear VV;clear UU1;clear VV1;

clear Res;clear FF;
%Pol=nan([size(d,1),1]);  %%% general array where to put all the polarisations
%Mn=nan([size(d,1),1]);  %%% general array where to put all the polarisations
Res.insert_name= strcat(path_dir,'/',subdir{jj});
Res.posrange=posrange{jj};

%%
for ii=1:size(d,1)  %%%% load all coordinates and velocities in 4 arrays (xx,yy,uu,vv)
    %% start loop for over the fov in the directory
    filename= d(ii).name;

    cd(path_dir);cd((subdir{jj}));
    fl_var= load(filename);
    u= fl_var.u;v= fl_var.v;u1= fl_var.u1;v1= fl_var.v1; x= fl_var.x;y= fl_var.y;
    um=nanmean(u,3); vm=nanmean(v,3); u1m=nanmean(u1,3); v1m=nanmean(v1,3);
    fl_var.mo;

    indxL=strfind(filename,'_X')+2;
    posx=[];cc=0; while ~isempty(str2num(filename(indxL+cc))) |filename(indxL+cc)=='-' ; 
    posx=strcat(posx,filename(indxL+cc));cc=cc+1;end;posx=str2num(posx);
    indyL=strfind(filename,'Y')+1;
    posy=[];cc=0; while  ~isempty(str2num(filename(indyL+cc)))| filename(indyL+cc)=='-'; 
    posy=strcat(posy,filename(indyL+cc));cc=cc+1;end;cc=0;posy=str2num(posy);
    if posx <8500 
    if  isa(Res.posrange,'char'); if Res.posrange=='all';Res.fov(ii).good=1;end
    elseif posx>=Res.posrange(1) & posx<=Res.posrange(2); Res.fov(ii).good=1;
    else  Res.fov(ii).good=0;   
    end
    Res.fov(ii).posx=posx;Res.fov(ii).posy=posy; Res.fov(ii).filename=filename;
    Res.fov(ii).file_mat= strcat(path_dir,'/',subdir{jj},'/',filename);

    %% load the std from BF videos and calculate mask and frequency map
  
    
%     if isfile('Res.mat'); old_Res=load('Res.mat');
%         ind=old_Res.Res.fov(ii).ind;
%         F= old_Res.Res.fov(ii).freq.F;
%     else
       
    cd(path_dir);cd(BFdir{jj});
    d_BF=dir(strcat('*X',num2str(posx),'*Y',num2str(posy),'*.movie') );
    mo=moviereader(d_BF.name);
    f_lim=[5,25];prob_cilia=0.5;box_size=4;area_min=10;
    [F_temp,s,BW_temp] = find_cilia(mo,f_lim,prob_cilia,box_size);
    [F4,BW] = remove_debris(F_temp,area_min,f_lim);
    ss=imresize(BW,box_size,'nearest');    
    IF= imresize(F4,box_size,'nearest');
    
    if size(ss,2)<fl_var.mo.width; dummy=zeros([fl_var.mo.height,fl_var.mo.width]); 
        dummy_fs=zeros([fl_var.mo.height,fl_var.mo.width,size(fs,3)],'uint8');
        dummy_fs(:,end-size(ss,2)+1:end,:)=fs; fs=dummy_fs;
        dummy(:,end-size(ss,2)+1:end)= ss; ss=logical(dummy);
    end; 
    
       if allfov{jj}==1; ss=logical(ones(size(ss)));end;
    
    box_ss=[];F32=[];bs=16; n_good_px= 0.2;
    for b=1:numel(x);
        y_box=floor(y(b)-bs);
        x_box=floor(x(b)-bs);
        temp_ss =ss(floor(y(b)-bs):floor(y(b)+bs),floor(x(b)-bs):floor(x(b)+bs));
        temp_f32=IF(floor(y(b)-bs):floor(y(b)+bs),floor(x(b)-bs):floor(x(b)+bs));
        box_ss(b)=mean(temp_ss(:));
        if box_ss(b)>n_good_px
            F32(b)=nanmedian(temp_f32(:));
        else F32(b)=nan;
        end
    end
    ind= box_ss>n_good_px;
    
    cd(path_dir);cd((subdir{jj}));
    F32=reshape(F32,size(x));
    IF(isnan(IF))=0;
    
     Res.fov(ii).freq.F4=F4; Res.fov(ii).freq.BW =BW; Res.fov(ii).freq.IF=IF;
     Res.fov(ii).freq.F_temp=F_temp; Res.fov(ii).freq.BW_temp =BW_temp;
     Res.fov(ii).freq.area_min=area_min;Res.fov(ii).freq.f_lim=f_lim;
     
    %% find polarisation fov  
    M= sqrt(u1m(ind).^2 +v1m(ind).^2);Res.fov(ii).Mm=median(M);
    nu= u1m(ind)./M; nv= v1m(ind)./M;
    pol=sqrt(nanmean(nu)^2+ nanmean(nv)^2); polx= nanmean(nv); 
    Res.fov(ii).Pol= pol;   Res.fov(ii).x=x;   Res.fov(ii).y=y; Res.fov(ii).u1m=u1m; Res.fov(ii).v1m=v1m;
    Res.fov(ii).ind=ind;    Res.fov(ii).nu=nu; Res.fov(ii).nv=nv; Res.fov(ii).Polx= polx;

    %% make plot of the field of view and save it
    figure(5);
%    imshow(ss);hold on;
    imagesc(IF); colormap(gray);colorbar();
    hold on; axis image
    quiver(x(ind),y(ind),u1m(ind),v1m(ind),'r','LineWidth',2);
    xlabel('X [$\mu$m]','Interpreter','latex');
    ylabel('Y[$\mu$m]','Interpreter','latex');
    title(strcat('X',num2str(posx),' Y',num2str(posy),' good= ',num2str(Res.fov(ii).good),' Pol=',num2str(polx,2)),'Interpreter','latex');
    
    x0=0;y0=0;width=1000;height=500;
    set(gcf,'position',[x0,y0,width,height]);set(gca,'FontSize',15);
    saveas(gcf,strcat(filename,'_PIV_std.png'));
    close(5)


    %% find frequency correlation, fit with exponential and plot 
    box_size=4;
    [fx,fy]=meshgrid(1:box_size:(size(F4,2))*box_size,1: box_size: (size(F4,1))*box_size);
    bin_res= 32; 
%    [idr,idth,cc,ecc,n_cc] = corr_frequency_theta_func(fx(BW),fy(BW),F4(BW),bin_res);
    [idr,idth,cc,ecc,n_cc] = corr_frequency_theta_func(x(ind),y(ind),F32(ind),bin_res);

    Res.fov(ii).freq.idr =idr; Res.fov(ii).freq.idth =idth; Res.fov(ii).freq.cc=cc; Res.fov(ii).freq.ecc=ecc; Res.fov(ii).freq.n_cc=n_cc;
    Res.fov(ii).freq.F32=F32;

    if numel(idr)>6  & ~any(isnan(cc(1:20)))  
    
        
    ft = fittype('a*exp(-x/b)+c','independent','x');
    fo = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[0,1,-1],...
    'Upper',[1,Inf,1],...
    'StartPoint',[0.9, idr(2)*10, 0]);
 
    fit_bins= idr; fit_cc=cc;max_cc=20;%floor(numel(cc(idth==0))/2);
    fit_out = fit(fit_bins(1:max_cc),fit_cc(1:max_cc),ft,fo);
    Res.fov(ii).freq.fit_out=fit_out;

    figure(4)
    plot(fit_out,idr,cc,'k^');hold on;
    else
    figure(4)
    plot(idr,cc,'k^');hold on;Res.fov(ii).freq.fit_out=nan;
    end      
    %% find orientation correlation fov

    bin_res= 32; 
    [idr,idth,cc,ecc,n_cc] = corr_orientation_theta_func(x(ind),y(ind),nu,nv,bin_res);
    Res.fov(ii).idr =idr; Res.fov(ii).idth =idth; Res.fov(ii).cc=cc; Res.fov(ii).ecc=ecc; Res.fov(ii).n_cc=n_cc;
    
    %%%%% fit correlation fucntion with exponential to find tipic lenght 
    if numel(idr)>6
    
    ft = fittype('a*exp(-x/b)+c','independent','x');
    fo = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[0,idr(2),-1],...
    'Upper',[1,Inf,1],...
    'StartPoint',[0.9, idr(2)*10, 0]);
 
    fit_bins= idr; fit_cc=cc;max_cc=floor(numel(cc)/2);
    fit_out = fit(fit_bins(1:max_cc),fit_cc(1:max_cc),ft,fo);
    Res.fov(ii).fit_out=fit_out;
   
    figure(4)
     plot(fit_out,idr,cc,'ko');hold on;
    else
     figure(4)
     plot(idr,cc,'ko');hold on;  
     Res.fov(ii).fit_out=nan;
    end   
     legend({'freq_corr','orientation_corr'});
     xlabel('pixel [$\mu$m]','Interpreter','latex');
     ylabel(' correlation','Interpreter','latex');
     title('Orientation correlation function','Interpreter','latex');
     
    x0=0;y0=0;width=1000;height=500;
    set(gcf,'position',[x0,y0,width,height]);set(gca,'FontSize',15);
    saveas(gcf,strcat(filename,'_correlations.png'));
    close(4)

        
    %% total flow properties

    xx(start+1:start+numel(x(ind))) = x(ind) +posy/ppm;
    yy(start+1:start+numel(x(ind))) = y(ind) -posx/ppm;

    
    uu(start+1:start+numel(x(ind))) = um(ind);
    vv(start+1:start+numel(x(ind))) = vm(ind);
    uu1(start+1:start+numel(x(ind))) = u1m(ind);
    vv1(start+1:start+numel(x(ind))) = v1m(ind);
    FF(start+1:start+numel(x(ind))) = F32(ind);  %%%%% load all the frequencies

%%%%%%%copiato
       
    start= start+ numel(x(ind));
    end
end
Res.xx=xx;Res.yy=yy;Res.uu1=uu1;Res.vv1=vv1; Res.FF=FF;

XX= xx; YY= yy; UU= uu; VV= vv; UU1= uu1; VV1= vv1;

M= sqrt(UU1.^2 +VV1.^2);nu= UU1./M;nv= VV1./M;
Res.nu=nu;Res.nv=nv;Res.M=M;
Mm=median(M);
pol=sqrt(nanmean(nu)^2+ nanmean(nv)^2);
polx= nanmean(nu); 
UM= nanmean(nu);
VM= nanmean(nv);
Res.Pol=pol;
Res.Polx=polx;


figure(2);
subplot(2,2,1)
quiver(XX,YY,UU1,VV1);
xlabel('X [um]');
ylabel(' Y[um]');
title('Flow pattern with validation');
axis image

subplot(2,2,2)
quiver(XX,YY,nu,nv,.6);
xlabel('X [um]');
ylabel(' Y[um]');
title('Flow pattern validation and normalisation');
axis image


subplot(2,2,4)
rose(angle(complex(nu,nv)),18); 
title(strcat('Orientation distribution;  mean pol=',num2str(polx))); 

x0=0;y0=0;width=1000;height=1000;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',15);
    %%
bin_res=32;
%[bins,cc,ecc,n_cc] = corr_orientation_func(XX,YY,nu,nv,bin_res);
[idr,idth,cc,ecc,n_cc] = corr_orientation_theta_func(XX,YY,nu,nv,bin_res);
Res.idr =idr; Res.idth =idth; Res.cc=cc; Res.ecc=ecc; Res.n_cc=n_cc;
       
    ft = fittype('a*exp(-x/b)+c','independent','x');
    fo = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[0,1,-1],...
    'Upper',[1,Inf,1],...
    'StartPoint',[0.9, idr(2)*10, 0]);
 
    fit_bins= idr; fit_cc=cc;max_cc=20;%floor(numel(cc(idth==0))/2);
    fit_out = fit(fit_bins(1:max_cc),fit_cc(1:max_cc),ft,fo);
    Res.fit_out=fit_out;

    subplot(2,2,3)
    plot(fit_out,idr,cc,'ko');hold on;
    xlim([0,bin_res*40]);ylim([-1,1]);
    xlabel('pixel [$\mu$m]','Interpreter','latex');
    ylabel('orientation correlation','Interpreter','latex');
    title('Orientation correlation function total','Interpreter','latex');


    %% total frequency correlation plot and save into variable
    bin_res=32;
    [idr,idth,cc,ecc,n_cc] = corr_frequency_theta_func(XX,YY,FF,bin_res,pi/2);
    Res.freq.idr =idr; Res.freq.idth =idth; Res.freq.cc=cc; Res.freq.ecc=ecc; Res.freq.n_cc=n_cc;


    ft = fittype('a*exp(-x/b)+c','independent','x');
    fo = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[0,1,-1],...
    'Upper',[1,Inf,1],...
    'StartPoint',[0.9, idr(2)*10, 0]);
 
    fit_bins= idr; fit_cc=cc;max_cc=20;%floor(numel(cc(idth==0))/2);
    fit_out = fit(fit_bins(1:max_cc),fit_cc(1:max_cc),ft,fo);
    Res.freq.fit_out=fit_out;

 %   subplot(2,2,3)
    plot(fit_out,idr,cc,'^r');hold on;
    legend({'orientation','frequency'});
    xlabel('pixel [$\mu$m]','Interpreter','latex');
    ylabel('spatial correlation','Interpreter','latex');
    title('Frequency and orientation correlation function total','Interpreter','latex');
    xlim([0,bin_res*40]);ylim([-1,1]);
   
    x0=0;y0=0;width=1000;height=1000;
    set(gcf,'position',[x0,y0,width,height])
    set(gca,'FontSize',15);


    saveas(gcf,'flow_pattern.png');
    close all;
   

    save('Res3.mat','Res');

    cd(path_dir);
end
%% for creating binary maps to identify cells 






%% do all standard deviation and then cretae threshold std in ss; 


%long= {'0/BF','5A/BF','5B/BF'}%,'5C/BF'};
longBF={...
    '22.3.19/beads/40rpm/BF','22.3.19/beads/40rpm2/BF','22.3.19/beads/control/BF'...
    '29.3.19/beads/40rpm/BF','29.3.19/beads/40rpm2/BF','29.3.19/beads/control/BF'...
    '5.4.19/beads/40rpm/BF','5.4.19/beads/40rpm2/BF', '5.4.19/beads/control/BF'...
    '12.4.19/beads/40rpm/BF','12.4.19/beads/40rpm2/BF','12.4.19/beads/control/BF'...    
}

prob_cilia_insert={...
0.96,0.96,0.96...
0.95,0.95,0.95...
0.95,0.95,0.95...
0.95,0.95,0.97...    
}



for jj=1:numel(longBF)

%% load files
insert=longBF{jj};
data_dir= strcat('/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/',insert)
a_folder=strcat('/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/',insert);

mkdir(a_folder); 

cd(data_dir) 
direc = dir('*X*Y*.movie');
N_files= size(direc,1);

%%

for i=1:N_files
    cd(data_dir)
    exp_name = direc(i).name;
    exp_name=exp_name(1:end-6);

%    if exist(strcat(a_folder,'/',exp_name)) == 0

        mo=moviereader(direc(i).name);frameload=round(2*mo.FrameRate/20)
        fs=mo.read([30,30+frameload]); 
        s=std(double(fs),[],3);s=mat2gray(s);
        sm=medfilt2(s,[5,5]);sm=wiener2(sm,[5,5]);
        [pf,edg]=histcounts(sm(:),'Normalization','probability');
        cpf=cumsum(pf); prob_cilia= prob_cilia_insert{jj}
        edg=edg(2:end); T=edg(cpf>prob_cilia ); T=T(1);
        BW= imbinarize(sm,T);
        figure(1);
        subplot(1,2,1)
        title(strcat(exp_name,' standard deviation'));
        imagesc(s);
%        axis equal
        subplot(1,2,2)
        title(strcat(exp_name,' threshold=',num2str(T)));
        imagesc(BW);
        x0=0;
        y0=0;
        width=1200;
        height=400;
        set(gcf,'position',[x0,y0,width,height])
        fig=figure(1); 
        saveas(fig,strcat(exp_name,'_mask_std.png'));
        close(1); 

        clear fs; 
        save(strcat(exp_name,'.mat'));

end

end




