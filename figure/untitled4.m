%% figure all flow 
clear all
%%%%%%%%%% Load data to find mean polarisation over multiple field of view

%%%%%% this script is suitable with PIV_script_maker
%%%% it read in order of date
ppm= 0.14*2;  %%% 20X
%path_dir = '/home/np451/Desktop/ependymal/June2018/July/16.7.18/PIV/';
path_dir = '/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019'
cd(path_dir); 

subdir={...   
  '22.3.19/beads/40rpm/FL','22.3.19/beads/40rpm2/FL','22.3.19/beads/control/FL'...
    '29.3.19/beads/40rpm/FL','29.3.19/beads/40rpm2/FL','29.3.19/beads/control/FL'...
   '5.4.19/beads/40rpm/FL','5.4.19/beads/40rpm2/FL','5.4.19/beads/control/FL'...
    '12.4.19/beads/40rpm/FL','12.4.19/beads/40rpm2/FL'...,'12.4.19/beads/control/FL'... 
}; subdir=strcat('PIV/',subdir)


posrange={...
    'all','all','all'...
    [2112,4824],[0,4000],[0,4000]...
    'all','all','all'...
    [0,4000],[0,3000],...'none'
    }

BFdir={...   
   '22.3.19/beads/40rpm/BF','22.3.19/beads/40rpm2/BF','22.3.19/beads/control/BF'...
   '29.3.19/beads/40rpm/BF','29.3.19/beads/40rpm2/BF','29.3.19/beads/control/BF'...
   '5.4.19/beads/40rpm/BF','5.4.19/beads/40rpm2/BF','5.4.19/beads/control/BF'...
   '12.4.19/beads/40rpm/BF','12.4.19/beads/40rpm2/BF',...'12.4.19/beads/control/BF'... 
};% BFdir=strcat('PIV/',BFdir)

for jj=5:numel(subdir)
    
start=0;
cd(path_dir);cd(subdir{jj}); d=dir('*FL*.mat');


%%
for ii=1:size(d,1)  %%%% load all coordinates and velocities in 4 arrays (xx,yy,uu,vv)
    %% start loop for over the fov in the directory
    filename= d(ii).name;

    cd(path_dir);cd((subdir{jj}));
    fl_var= load(filename);
    fl_var.mo;
    
    
    cd(fl_var.data_dir);
    d2=dir(strcat(filename(1:end-5),'*.movie'));
    mo=moviereader(d2(1).name);
    fs=mo.read();
    s=std(double(fs),[],3);
    cd(path_dir);cd((subdir{jj}));
    
    figure();
    subplot(1,2,1);
    imagesc(imadjust(mat2gray(s)));title(filename)
    axis equal
    title('')

    
    
    indxL=strfind(filename,'_X')+2;
    posx=[];cc=0; while ~isempty(str2num(filename(indxL+cc))) |filename(indxL+cc)=='-' ; 
    posx=strcat(posx,filename(indxL+cc));cc=cc+1;end;posx=str2num(posx);
    indyL=strfind(filename,'Y')+1;
    posy=[];cc=0; while  ~isempty(str2num(filename(indyL+cc)))| filename(indyL+cc)=='-'; 
    posy=strcat(posy,filename(indyL+cc));cc=cc+1;end;cc=0;posy=str2num(posy);


    
    %% load the std from BF videos and calculate mask and frequency map
  

    cd(path_dir);cd(BFdir{jj});
    d_BF=dir(strcat('*X',num2str(posx),'*Y',num2str(posy),'*.movie') );
    mo=moviereader(d_BF.name);
    fs=mo.read([30,70]);
    subplot(2,2,2);
    s=std(double(fs),[],3);
    cd(path_dir);cd((subdir{jj}));
    
    imagesc(imadjust(mat2gray(s)));title(filename);
    axis equal
    
    x0=0;y0=0;width=1000;height=500;

    set(gcf,'position',[x0,y0,width,height])
    set(gca,'FontSize',15);
    saveas(gcf,strcat(filename,'std_beads_and_bf.pdf'));
    close all
    
end
end