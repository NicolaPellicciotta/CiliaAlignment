% script for figure 1 of the paper 

subdir={...   
  '22.3.19/beads/40rpm/FL','22.3.19/beads/40rpm2/FL','22.3.19/beads/control/FL'...
    '29.3.19/beads/40rpm/FL','29.3.19/beads/40rpm2/FL','29.3.19/beads/control/FL'...
   '5.4.19/beads/40rpm/FL','5.4.19/beads/40rpm2/FL','5.4.19/beads/control/FL'...
    '12.4.19/beads/40rpm/FL','12.4.19/beads/40rpm2/FL'...,'12.4.19/beads/control/FL'... 
};
subdir = strcat('/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/',subdir);
 

BFdir={...   
   '22.3.19/beads/40rpm/BF','22.3.19/beads/40rpm2/BF','22.3.19/beads/control/BF'...
   '29.3.19/beads/40rpm/BF','29.3.19/beads/40rpm2/BF','29.3.19/beads/control/BF'...
   '5.4.19/beads/40rpm/BF','5.4.19/beads/40rpm2/BF','5.4.19/beads/control/BF'...
   '12.4.19/beads/40rpm/BF','12.4.19/beads/40rpm2/BF',...'12.4.19/beads/control/BF'... 
}
BFdir = strcat('/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/',BFdir);
 


for jj=1:numel(subdir)
cd(subdir{jj}); d=dir('*FL*.mat');

for ii=1:size(d,1)
    filename= d(ii).name;
    cd((subdir{jj}));
    fl_var= load(filename);
    cd(fl_var.data_dir);
    d2=dir(strcat(filename(1:end-5),'*.movie'));
    mo=moviereader(d2(1).name);
    fs=mo.read();
    s=std(double(fs),[],3);
    cd(subdir{jj})
    figure();imagesc(imadjust(mat2gray(s)));title(filename)
    axis equal
    saveas(gcf,strcat(filename,'std_beads.pdf'));
    close all
end
end

