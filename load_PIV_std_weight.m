%% figure all flow 

%%%%%%%%%% Load data to find mean polarisation over multiple field of view

%%%%%% this script is suitable with PIV_script_maker
%%%% it read in order of date
ppm= 0.14*2;  %%% 20X
%path_dir = '/home/np451/Desktop/ependymal/June2018/July/16.7.18/PIV/';
path_dir = '/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019'
cd(path_dir); 
%subdir= {'0/FL','30A/FL','30B/FL','100A/FL','100B/FL','1000/FL'};
%BFdir={'0/BF','30A/BF','30B/BF','100A/BF','100B/BF','1000/BF'};
%subdir= {'0/FL','100/FL','200/FL'};
%BFdir={'0/BF','100/BF','200/BF'};


subdir={...
%    '5.4.19/beads/40rpm/FL','5.4.19/beads/40rpm2/FL','5.4.19/beads/control/FL'...
%    '12.4.19/beads/40rpm/FL','12.4.19/beads/40rpm2/FL','12.4.19/beads/control/FL'...    
%   '22.3.19/beads/40rpm/FL','22.3.19/beads/40rpm2/FL','22.3.19/beads/control/FL'...
    '29.3.19/beads/40rpm/FL','29.3.19/beads/40rpm2/FL','29.3.19/beads/control/FL'...
}; subdir=strcat('PIV/',subdir)

%subdir= {'0/FL','5A/FL','5B/FL','5C/FL'};
%BFdir={'0/BF','5A/BF','5B/BF','5C/BF'};
BFdir={...
%    '5.4.19/beads/40rpm/BF','5.4.19/beads/40rpm2/BF','5.4.19/beads/control/BF'...
%    '12.4.19/beads/40rpm/BF','12.4.19/beads/40rpm2/BF','12.4.19/beads/control/BF'...    
%   '22.3.19/beads/40rpm/BF','22.3.19/beads/40rpm2/BF','22.3.19/beads/control/BF'...
    '29.3.19/beads/40rpm/BF','29.3.19/beads/40rpm2/BF','29.3.19/beads/control/BF'...
}

%subdir={'5A/FL','5B/FL','5C/FL','15A/FL','15B/FL','15C/FL'};
%BFdir={'5A/BF','5B/BF','5C/BF','15A/BF','15B/BF','15C/BF'};

%BFdir = {'100B/BF'};

for jj=1:numel(subdir)
    
start=0;
cd(path_dir);cd(subdir{jj}); d=dir('*FL*.mat');

clear xx;clear yy;clear uu;clear vv;clear uu1;clear vv1;
clear XX;clear YY;clear UU;clear VV;clear UU1;clear VV1;


for ii=1:size(d,1)  %%%% load all coordinates and velocities in 4 arrays (xx,yy,uu,vv)
    filename= d(ii).name;

    cd(path_dir);cd((subdir{jj}));
    load(filename);
%    subdir= {'0/FL','100/FL','200/FL'};
%    BFdir={'0/BF','100/BF','200/BF'};
    cd(path_dir);cd(subdir{jj});
    d=dir(strcat('*FL*.mat'));
    
    indxL=strfind(filename,'_X')+2;
    posx=[];cc=0; while ~isempty(str2num(filename(indxL+cc))) |filename(indxL+cc)=='-' ; 
    posx=strcat(posx,filename(indxL+cc));cc=cc+1;end;posx=str2num(posx);
    indyL=strfind(filename,'Y')+1;
    posy=[];cc=0; while  ~isempty(str2num(filename(indyL+cc)))| filename(indyL+cc)=='-'; 
    posy=strcat(posy,filename(indyL+cc));cc=cc+1;end;cc=0;posy=str2num(posy);
% 
%     indxR=strfind(filename,'Y')-1;
%     indyR=strfind(filename,'16Jul')-2;
%     posx= str2num(filename(indxL :indxR));
%     posy= str2num(filename(indyL:end-23));
    
    %%%%%load the std from BF videos
    cd(path_dir);cd(BFdir{jj});
    d_BF=dir(strcat('*X',num2str(posx),'*Y',num2str(posy),'*.mat') );
    bf_var= load(d_BF.name);

    cd(path_dir);cd((subdir{jj}));
    
    %bf_var.ss(ss>0)=1;
    ss=bf_var.BW;
%    if size(ss,2)<mo.width; dummy=zeros([mo.height,mo.width]); 
%        dummy(:,end-size(ss,2)+1:end)= ss;end; ss=dummy;
%    if size(ss,2)>mo.width; dummy=zeros([mo.height,mo.width]); 
%        dummy=ss(:,end-size(ss,2)+1:end);end; ss=dummy;
    bs=16;
    box_ss=[];
    for b=1:numel(x);
        y_box=floor(y(b)-bs);
        x_box=floor(x(b)-bs);
        temp_ss =ss(floor(y(b)-bs):floor(y(b)+bs),floor(x(b)-bs):floor(x(b)+bs),:);
        box_ss(b)=mean(temp_ss(:));
    end
    
    n_good_px= 0.1
    ind= box_ss>n_good_px;
        
    xx(start+1:start+numel(x(ind))) = x(ind) +posy/ppm;
    yy(start+1:start+numel(x(ind))) = y(ind) -posx/ppm;
    um=nanmean(u,3);
    vm=nanmean(v,3);
    u1m=nanmean(u1,3);
    v1m=nanmean(v1,3);
    
    
    uu(start+1:start+numel(x(ind))) = um(ind);
    vv(start+1:start+numel(x(ind))) = vm(ind);
    uu1(start+1:start+numel(x(ind))) = u1m(ind);
    vv1(start+1:start+numel(x(ind))) = v1m(ind);

%%%%%%%copiato
       
    start= start+ numel(x(ind));
end


 XX= xx;
 YY= yy;
 UU= uu;
 VV= vv;
 UU1= uu1;
 VV1= vv1;

% 
% %%%%% SET new lattice with space dx and do a velocy average over a radius dx/2
% 
% dx=64;
% px_max=(floor(max(xx(:))/dx)+1)*dx;
% px_min=(floor(min(xx(:))/dx)+1)*dx;
% %%%% set the new axis limit
% py_max=(floor(max(yy(:))/dx)+1)*dx;
% py_min=(floor(min(yy(:))/dx)+1)*dx;
% 
% XX= repmat((px_min:dx:px_max),[numel((py_min:dx:py_max)),1]);
% YY= repmat((py_min:dx:py_max)',[1,numel((px_min:dx:px_max))]);
% UU= zeros(size(XX));
% VV= zeros(size(YY));
% UU1= zeros(size(XX));
% VV1= zeros(size(YY));
% 
% 
% min_d=dx;
% 
% 
% for i=1:size(XX,1)
%     for kk=1:size(XX,2)
%         dr= sqrt( (xx(:)-XX(i,kk)).^2 + (yy(:)-YY(i,kk)).^2 ) ;
% 
%         if numel(dr(:)<min_d)~=0
%             UU1(i,kk) = nanmean(uu1(dr(:)<min_d));
%             VV1(i,kk) = nanmean(vv1(dr(:)<min_d));
%             UU(i,kk) = nanmean(uu(dr(:)<min_d));
%             VV(i,kk) = nanmean(vv(dr(:)<min_d));
%         else
%             UU1(i,kk) = NaN;
%             VV1(i,kk) = NaN;
%             UU(i,kk) = NaN;
%             VV(i,kk) = NaN;
%         end
%     end
% end 

%
M= sqrt(UU1.^2 +VV1.^2);
nu= UU1./M;
nv= VV1./M;

figure(2);
subplot(2,2,1)
quiver(XX,YY,UU1,VV1);
xlabel('X [um]');
ylabel(' Y[um]');
title('Flow pattern with validation');
axis image
%saveas(gcf,'FlowPattern.pdf');
%close all;


% figure();
% quiver(XX,YY,UU,VV,1.4);
% xlabel('X [um]');
% ylabel(' Y[um]');
% title('Flow pattern without validation');
% axis image
% saveas(gcf,'FlowPattern_noval.pdf');
% close all;

subplot(2,2,2)
quiver(XX,YY,nu,nv,.6);
xlabel('X [um]');
ylabel(' Y[um]');
title('Flow pattern validation and normalisation');
axis image
%saveas(gcf,'FlowPattern_norm.pdf');
%close all;



%%% correlation 

fXX= XX;
fYY= YY;
fUU= nu-mean(nu(:));
fVV= nv-mean(nv(:));




dist_x= repmat(fXX(:), [1,numel(fXX(:))]);
dist_y= repmat(fYY(:), [1,numel(fYY(:))]);
DR= sqrt((dist_x - dist_x').^2 + (dist_y - dist_y').^2) ;

%CORR = repmat(a(:),[1,numel(a)]) - (repmat(a(:),[1,numel(a)]))'; 
CORRu  = repmat(fUU(:), [1,numel(fUU(:))]);
CORRv  = repmat(fVV(:), [1,numel(fVV(:))]);

CORR = (CORRu.*CORRu' + CORRv.*CORRv')./( sqrt(CORRu.^2 + CORRv.^2).*sqrt(CORRu'.^2 + CORRv'.^2)); 
bin_res = 96;
bins = 0:bin_res:max(DR(:));

clear cc;
clear n;
for i=1:(numel(bins)-1)

       CORR_bin= (CORR( (DR(:) > bins(i)) & (DR(:) <= bins(i+1))));
       cc(i)= mean(CORR_bin(~isnan(CORR_bin)));
       ecc(i) = mean((CORR_bin(~isnan(CORR_bin))- cc(i)).^2);
       n(i) = numel(~isnan(CORR_bin));
       ecc(i)=sqrt(ecc(i)/n(i));
end

subplot(2,2,3)
plot((1:numel(cc))*bin_res*ppm,cc,'o');
xlabel(' [px]');
ylabel(' correlation[um]');
title('Orientation correlation function');
%axis image
%saveas(gcf,'orientation_correlation.pdf');
%close all;
% 
% %%%%%% end correlation
% 
% %%%%% polarisation
% 
M= sqrt(UU1.^2 +VV1.^2);
nu= UU1./M;
nv= VV1./M;
Mm=median(M);
pol=sqrt(nanmean(nu)^2+ nanmean(nv)^2);
polx= nanmean(nu); 
UM= nanmean(nu);
VM= nanmean(nv);

%%%%% end polarisation 

subplot(2,2,4)
rose(angle(complex(nu,nv)),18); 
title(strcat('Orientation distribution;  mean pol=',num2str(polx))); 

x0=0;
y0=0;
width=1000;
height=1000;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',15);


saveas(gcf,'flow_pattern.png');
close all;
% 
% 
% 
 save('post.mat')


cd(path_dir);
end

%% this is the script to find correlation between the number of first neighbours and the polarisation 

ppm= 0.14*2;  %%% 20X


%path_dir = '/home/np451/Desktop/ependymal/June2018/13.6.18/PIV/';
path_dir = '/home/np451/Desktop/ependymal/June2018/July/12.7.18/PIV/';
cd(path_dir); 

%subdir= {'0/FL','30A/FL','30B/FL','100A/FL','100B/FL','1000/FL'};
%subdir= {'0/FL','30B/FL','100A/FL','1000/FL'};
subdir= {'0F/FL'}%,'5A/FL_backup','5B/FL','5C/FL'};
%subdir={'0F/FL','5A/FL_backup','5B/FL','5C/FL','15A/FL','15B/FL','15C/FL'};

for jj=1:numel(subdir)

cd(path_dir);cd(subdir{jj}); d=dir('*FL*.mat');
a=load('post.mat');
%subdir= {'0/FL','30A/FL','30B/FL','100A/FL','100B/FL','1000/FL'};
%subdir= {'0/FL','30B/FL','100A/FL','1000/FL'};
%subdir= {'0/FL','100/FL','200/FL'};
fXX= a.XX;
fYY= a.YY;
%fUU= a.nu*0;%UU1;
%fVV= a.nv;%VV1;
%%%%% for polarisation only along one axis
nu=a.nu*0;
nv=a.nv;
%%%% for mean velocity
nu=a.UU1;
nv=a.VV1;


d_lim= 500;  %%%% px? 

dist_x= repmat(fXX(:), [1,numel(fXX(:))]);
dist_y= repmat(fYY(:), [1,numel(fYY(:))]);
DR= sqrt((dist_x - dist_x').^2 + (dist_y - dist_y').^2) ;

N_near=sum(DR<d_lim,2);

%near_bin= min(N_near)+ (1:6)*(max(N_near)-min(N_near))/6;

%near_bin= mean(N_near)-std(N_near) + (1:8)*(std(N_near)/4);
%near_bin= (0:10)*(mean(N_near)+ std(N_near))/8;
near_bin=10:10:100%5:3:50;

for n=1:(numel(near_bin)-1)
    near_ind= N_near > near_bin(n) & N_near <= near_bin(n+1);
    near_pol(n)= sqrt(mean(nu(near_ind))^2 + mean(nv(near_ind))^2);
    NN_near(n)=sum(near_ind);
    
    near_e(n)=sqrt(2*mean(nu(near_ind))^2 *(std(nu(near_ind))/sqrt(NN_near(n)) )^2 ...
    + 2*mean(nv(near_ind))^2 *(std(nv(near_ind))/sqrt(NN_near(n)))^2);
        
end

bin = (near_bin(1:end-1)+near_bin(2:end))/2;



figure(1); errorbar(bin,near_pol,near_e ,'-o'); hold on;
title(subdir{jj});

figure(2); plot(bin,NN_near,'o'); hold on;
title(subdir{jj});


end
%leg={'0','30-A','30-B','100-A','100-B','200'};
%leg={'0','30','100','200'};
%leg={'0/FL','5A/FL','5B/FL','5C/FL'};
%leg={'5A/FL','5B/FL','5C/FL','15A/FL','15B/FL','15C/FL'};
leg=subdir;
figure(1);legend(leg);
xlabel('number of near ciliated cells');ylabel('polarisation');title('');
figure(2);legend(leg);
%% Divide the channel in square and calculate local properties 

%%%%% SET new lattice with space dx and do a velocy average over a radius dx/2
ppm= 0.14*2;%%% 20X
fps=10;


%path_dir = '/home/np451/Desktop/ependymal/June2018/13.6.18/PIV/';
path_dir = '/home/np451/Desktop/ependymal/June2018/July/16.7.18/PIV/';
cd(path_dir); 
subdir= {'0/FL','5B/FL','5C/FL'};
%subdir= {'5B/FL'};

for jj=1:numel(subdir)

cd(path_dir);cd(subdir{jj}); d=dir('*FL*.mat');
a=load('post.mat');

xx= a.XX;
yy= a.YY;

%%%%% load normalised velocity
nu=a.nu*0;
nv=a.nv;
%%%% load velocity
uu1=a.UU1;
vv1=a.VV1;

dx=250;
px_max=(floor(max(xx(:))/dx)+1)*dx;
px_min=(floor(min(xx(:))/dx)+1)*dx;
%%%% set the new axis limit
py_max=(floor(max(yy(:))/dx)+1)*dx;
py_min=(floor(min(yy(:))/dx)+1)*dx;

XX= repmat((px_min:dx:px_max),[numel((py_min:dx:py_max)),1]);
YY= repmat((py_min:dx:py_max)',[1,numel((px_min:dx:px_max))]);
UL= zeros(size(XX));    %%% velocity along 
VL= zeros(size(YY)); %%% velocity along x
CL= zeros(size(XX)); %%% concentration in each box
PL = zeros(size(XX)); %%% polarisation in each box

min_d=dx;

for i=1:size(XX,1)
    for kk=1:size(XX,2)
        dr= sqrt( (xx(:)-XX(i,kk)).^2 + (yy(:)-YY(i,kk)).^2 ) ;

        if numel(dr(:)<min_d)~=0
            UL(i,kk) = nanmean(uu1(dr(:)<min_d));
            VL(i,kk) = nanmean(vv1(dr(:)<min_d));
            CL(i,kk) = sum((dr(:)<min_d));
            PL(i,kk) = sqrt( mean( nu((dr(:)<min_d)))^2 + mean( nv((dr(:)<min_d)))^2 );
        else
            UL(i,kk) = NaN; VL(i,kk) = NaN; CL(i,kk) = NaN; PL(i,kk) = NaN;
        end
    end
end 

MUV= sqrt(UL.^2 + VL.^2);

ind=CL(:)>10;
%%%% polarisation vs concentration of cells
nbins=10;
figure(1);
[histw,ehistw,vinterval]=hist_nico(CL(ind),PL(ind),nbins);
errorbar(vinterval,histw,ehistw,'o','MarkerSize',10,'LineWidth',1.5);hold on;
title('concentration vs polarisation')
xlabel('concentration of cells');
ylabel('polarisation');
%set(gca, 'YScale', 'log')

%%%% Flow velocity vs concentration of cells

nbins=10;
figure(2);
[histw,ehistw,vinterval]=hist_nico(CL(ind),MUV(ind)*ppm*fps ,nbins);
errorbar(vinterval,histw,ehistw,'o','MarkerSize',10,'LineWidth',1.5);hold on;
title('concentration vs velocity')
xlabel('concentration of cells');
ylabel('velocity');


%%%% polarisation vs velocity of cells

nbins=10;
figure(3);

[histw,ehistw,vinterval]=hist_nico(PL(ind),MUV(ind)*ppm*fps ,nbins);
errorbar(vinterval,histw,ehistw,'o','MarkerSize',10,'LineWidth',1.5); hold on;
title('polarisation vs velocity');
xlabel('polarisation');
ylabel('velocity');

end

leg=subdir;
legend(leg)
%% this script find the ciliated cell density in the channel and polarisation at different distance from the channel

ppm= 0.14*2;  %%% 20X
%path_dir = '/home/np451/Desktop/ependymal/June2018/13.6.18/PIV/';
path_dir = '/home/np451/Desktop/ependymal/June2018/16.7.18/PIV/';
cd(path_dir); 
%subdir= {'0/FL','30A/FL','30B/FL','100A/FL','100B/FL','1000/FL'};
%subdir= {'0/FL','100/FL','200/FL'};
%subdir={'5A/FL','5B/FL','5C/FL','15A/FL','15B/FL','15C/FL'};
subdir={'0/FL','5A/FL','5B/FL','5C/FL'};

xres=11;
x_slices= (1:xres)*3800/xres;
for sub=1:numel(subdir)

cd(path_dir);cd(subdir{sub});
load('post.mat');
%subdir= {'0/FL','30A/FL','30B/FL','100A/FL','100B/FL','1000/FL'};
%subdir= {'0/FL','100/FL','200/FL'};
for g=1:(numel(x_slices)-1)
    ind_x= (XX>x_slices(g) & XX<=x_slices(g+1));
    Pcilia(g)= numel(XX(XX>x_slices(g) & XX<=x_slices(g+1)));
    pol_xslice(g)= sqrt(mean(UU1(ind_x))^2 + mean(VV1(ind_x))^2); 
end


R{sub}.Pcilia=Pcilia/(sum(Pcilia(:)));
R{sub}.x_bin = (x_slices(1:end-1) + x_slices(2:end))./2;
R{sub}.density=numel(XX);
R{sub}.pol= pol_xslice;
clear Pcilia
end



%%%%%% plot polarisation at different distances from the channel
%leg={'0','30-A','30-B','100-A','100-B','200'};
%flow=[0,30,30,100,100,200];
%leg={'5A','5B','5C','15A','15B','15C'};
leg={'0','5A','5B','5C'};
flow=[5,5,5,15,15,15].*40;


figure(3)
for b=1:numel(R)
 plot(R{b}.x_bin,R{b}.Pcilia,'-o'); hold on;

end
xlabel('X position [px]'); ylabel('Cilia Probaility ditribution');

legend(leg)


%%%%%% plot polarisation at different distances from the channel

figure(4)
for b=1:numel(R)
 plot(R{b}.x_bin,R{b}.pol,'-o','LineWidth',2); hold on;

end
xlabel('X position [px]'); ylabel('polarisation');

legend(leg)


%%%%% plot density of cells vs experiment
figure(5)
for b=1:numel(R)
 plot(flow(b),R{b}.density,'o','MarkerSize',5,'LineWidth',3); hold on;

end
xlabel('flow [ul/min]'); ylabel('N of pixels with high threshold-std');
legend(leg)



%% PLOTS compare correlation function and mean orientation
%path_dir = '/home/np451/Desktop/ependymal/June2018/13.6.18/PIV/';
path_dir = '/home/np451/Desktop/ependymal/June2018/29.7.18/PIV/';

cd(path_dir); 

%subdir={'0ulmin','10ulmin','100ulmin'};
%subdir= {'0/FL','30A/FL','30B/FL','100A/FL','100B/FL','1000/FL'};
%flow=[0,30,30,100,100,200];
%subdir= {'0/FL','100/FL','200/FL'};
%flow=[0,100,200];
subdir={'0/FL','5A/FL','5B/FL'}%,'5C/FL'};
%subdir={'5A/FL','5B/FL','5C/FL','15A/FL','15B/FL','15C/FL'};
%flow=[5,5,5,15,15,15].*40;
flow=[0,5,5,5].*40;
for h=1:numel(subdir)
cd(path_dir);cd(subdir{h});
load('post.mat');
%subdir= {'0/FL','30A/FL','30B/FL','100A/FL','100B/FL','1000/FL'};
%subdir= {'0/FL','100/FL','200/FL'};
figure(1); plot((1:numel(cc))*bin_res*ppm,cc,'o','MarkerSize',5,'LineWidth',3); hold on;
figure(2); plot(flow(h), polx ,'o','MarkerSize',5,'LineWidth',3); hold on;
R{h}.mean_pol=polx;
R{h}.corr=cc;
end

cd(path_dir);

figure(1); xlabel(' distance [um]'); ylabel(' correlation[um]');
xlim([0,1000]);
%leg={'0','30-A','30-B','100-A','100-B','200'};
%leg={'0','100','200'};
%leg={'5A','5B','5C','15A','15B','15C'};
leg={'0','5A','5B'}%,'5C'};

legend(leg);
%saveas(gcf,'correlation_compare.pdf');

figure(2);xlabel('flow [ul/min]');ylabel('polarisation');
%saveas(gcf,'polarisation_compare.pdf');
legend(leg);


np(1)=R{1}.mean_pol;
np(2) = (R{2}.mean_pol+R{3}.mean_pol)/2;
np(3)=  (R{4}.mean_pol+R{5}.mean_pol)/2;
np(4)=R{6}.mean_pol;
figure();
plot([0,30,100,200], np,'o','LineWidth',2);
xlabel('flow [ul/min]');ylabel('polarisation');





%% do all standard deviation and then cretae threshold std in ss; 


%long= {'0/BF','5A/BF','5B/BF'}%,'5C/BF'};
longBF={...
%    '5.4.19/beads/40rpm/BF','5.4.19/beads/40rpm2/BF',...
    '5.4.19/beads/control/BF'...
    '12.4.19/beads/40rpm/BF','12.4.19/beads/40rpm2/BF','12.4.19/beads/control/BF'...    
   '22.3.19/beads/40rpm/BF','22.3.19/beads/40rpm2/BF','22.3.19/beads/control/BF'...
    '29.3.19/beads/40rpm/BF','29.3.19/beads/40rpm2/BF','29.3.19/beads/control/BF'...
}


for insert=longBF;

%% load files
insert=insert{1};
%data_dir = strcat('/run/user/10704/gvfs/smb-share:server=sf3.bss.phy.private.cam.ac.uk,share=space/np451/ependymalJune/4.7.18/',insert);
%data_dir =strcat('/media/np451/Seagate Expansion Drive/ependymal/27.7.18/',insert);
%a_folder= strcat('/home/np451/Desktop/ependymal/June2018/27.7.18/PIV/',insert);

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

        mo=moviereader(direc(i).name);
        fs=mo.read([30,30+frameload]); frameload=round(2*mo.FrameRate/20)
        s=std(double(fs),[],3);s=mat2gray(s);
        sm=medfilt2(s,[5,5]);sm=wiener2(sm,[5,5]);
        [pf,edg]=histcounts(sm(:),'Normalization','probability');
        cpf=cumsum(pf); prob_cilia= 0.97
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

        
%        axis equal
%        fig=figure(1); 
%        saveas(fig,strcat(exp_name,'_mask_std.png'));
%        close(1); 
        
%         %%interactive search for good threeshold for the std
%         disp('Select Noise');
%         BN=roipoly(mat2gray(s));
%         v_noise= mean(s(BN));
%         %%%% select focused cilia
%         disp('Select very focused bead');
%         BC=roipoly(mat2gray(s));
%         v_good= mean(s(BC));
%         p_val= v_good -abs((v_good - v_noise)/3) % -(abs(v_good-v_halo)/3);
%        p_val=1.5; 
%        ss=s;ss(ss<p_val)=0;ss(ss>p_val)=1;
        
%         cd(a_folder);       
% 
%         figure(1);
%         subplot(2,1,1)
%         title(strcat(exp_name,' standard deviation'));
%         imagesc(s);
%         axis equal
%         subplot(2,1,2)
%         title(strcat(exp_name,' p_val=',num2str(p_val)));
%         imagesc(ss);
%         axis equal
%         fig=figure(1); 
%         saveas(fig,strcat(exp_name,'_mask_std.png'));
%         close(1); 
% 
% 
%         clear fs; 
%         save(strcat(exp_name,'.mat'));

        
        
    
%    end
end

end

%%

%%%load std and clean the useless variables tha mess up everything
%%% do all standard deviation and then cretae threshold std in ss; 
%long= {'100A/BF','1000/BF','0/BF','30A/BF',
long= {'5A/BF','5B/BF','5C/BF','15A/BF','15B/BF','15C/BF'};
for insert=long;

%% load files
insert=insert{1};
a_folder= strcat('/media/np451/Seagate Expansion Drive/ependymal/12.7.18/',insert)
cd(a_folder) 
direc = dir('*BF*.mat');
N_files= size(direc,1);

%%

for i=1:N_files
    cd(a_folder)
    exp_name = direc(i).name;
    exp_name=exp_name(1:end-4);
    a=load(strcat(exp_name,'.mat'));
    s=a.s;ss=a.ss;
    save(strcat('std_',exp_name,'.mat'),'ss','s');
  
end

end


%% single figure to check validity of PIV filtered with std

ppm= 0.14*2  %%% 20X


%path_dir = '/home/np451/Desktop/ependymal/June2018/13.6.18/PIV/';
path_dir = '/u/homes/np451/Desktop/ependymal/June2018/June/18.6.18/PIV/';
cd(path_dir); 

subdir= {'0/FL'}
BFdir= {'0/BF'};


for jj=1:numel(subdir)

    
start=0;
cd(path_dir);cd(subdir{jj}); d=dir('*FL*.mat');


clear xx;clear yy;clear uu;clear vv;clear uu1;clear vv1;
clear XX;clear YY;clear UU;clear VV;clear UU1;clear VV1;


    ii=1;  %%%% load all coordinates and velocities in 4 arrays (xx,yy,uu,vv)
    filename= d(ii).name;

    cd(path_dir);cd((subdir{jj}));
    load(filename);
    cd(path_dir);cd(subdir{jj});
    d=dir(strcat('*FL*.mat'));
    
    indxL=strfind(filename,'_X')+2;
    indyL=strfind(filename,'Y')+1;
    indxR=strfind(filename,'Y')-1;
    indyR=strfind(filename,'18Jun')-2;
    posx= str2num(filename(indxL :indxR));
    posy= str2num(filename(indyL:indyR));
    
    %%%%%load the std from BF videos
    cd(path_dir);cd(BFdir{jj});
    d_BF=dir(strcat('std*X',num2str(posx),'Y',num2str(posy),'*.mat') );
    load(d_BF.name);
    cd(path_dir);cd((subdir{jj}));
    
    ss(ss>0)=1;
    
    bs=16;
    for b=1:numel(x);
        y_box=floor(y(b)-bs);
        x_box=floor(x(b)-bs);
        temp_ss =ss(floor(y(b)-bs):floor(y(b)+bs),floor(x(b)-bs):floor(x(b)+bs),:);
        box_ss(b)=mean(temp_ss(:));
    end
    
    ind= box_ss>0.07;
        
    xx(start+1:start+numel(x(ind))) = x(ind) 
    yy(start+1:start+numel(x(ind))) = y(ind) 
    um=nanmean(u,3);
    vm=nanmean(v,3);
    u1m=nanmean(u1,3);
    v1m=nanmean(v1,3);
    
    
    uu(start+1:start+numel(x(ind))) = um(ind);
    vv(start+1:start+numel(x(ind))) = vm(ind);
    uu1(start+1:start+numel(x(ind))) = u1m(ind);
    vv1(start+1:start+numel(x(ind))) = v1m(ind);

%%%%%%%copiato
       
    start= start+ numel(x(ind));

figure();
imshow(ss);hold on;
quiver(xx,yy,uu1,vv1,'r','LineWidth',2);
xlabel('X [um]');
ylabel(' Y[um]');
title('Flow pattern with validation');
axis image


cd(path_dir);
end





