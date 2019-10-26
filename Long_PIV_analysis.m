
%% Standard PIV Settings
s = cell(10,2); % To make it more readable, let's create a "settings table"
%Parameter                       %Setting           %Options
s{1,1}= 'Int. area 1';           s{1,2}=256;         % window size of first pass
s{2,1}= 'Step size 1';           s{2,2}=64;         % step of first pass
s{3,1}= 'Subpix. finder';        s{3,2}=1;          % 1 = 3point Gauss, 2 = 2D Gauss
s{4,1}= 'Mask';                  s{4,2}=[];         % If needed, generate via: imagesc(image); [temp,Mask{1,1},Mask{1,2}]=roipoly;
s{5,1}= 'ROI';                   s{5,2}=[];         % Region of interest: [x,y,width,height] in pixels, may be left empty
s{6,1}= 'Nr. of passes';         s{6,2}=3;          % 1-4 nr. of passes
s{7,1}= 'Int. area 2';           s{7,2}=128;         % second pass window size
s{8,1}= 'Int. area 3';           s{8,2}=64;         % third pass window size
s{9,1}= 'Int. area 4';           s{9,2}=32;         % fourth pass window size
s{10,1}='Window deformation';    s{10,2}='*linear'; % '*spline' is more accurate, but slower

%%% Standard image preprocessing settings
p = cell(8,1);
%Parameter                       %Setting           %Options
p{1,1}= 'ROI';                   p{1,2}=s{5,2};     % same as in PIV settings
p{2,1}= 'CLAHE';                 p{2,2}=0;          % 1 = enable CLAHE (contrast enhancement), 0 = disable
p{3,1}= 'CLAHE size';            p{3,2}=50;         % CLAHE window size
p{4,1}= 'Highpass';              p{4,2}=0;          % 1 = enable highpass, 0 = disable
p{5,1}= 'Highpass size';         p{5,2}=15;         % highpass size
p{6,1}= 'Clipping';              p{6,2}=0;          % 1 = enable clipping, 0 = disable
p{7,1}= 'Wiener';                p{7,2}=0;          % 1 = enable Wiener2 adaptive denaoise filter, 0 = disable
p{8,1}= 'Wiener size';           p{8,2}=3;          % Wiener2 window size

%%
%long= {'3_FL','5_FL','6_FL'};%,'3_BF','5_BF','6_BF'};
%long= {'0/FL','5B/FL'}%,'5B/FL'};
%long={'5A/FL','5B/FL','5C/FL','15A/FL','15B/FL','15C/FL'}
long={...
    ...'11.10.19/control/FL/',
 %   '11.10.19/5rpm/FL/','11.10.19/5rpm_2/FL/'...
...    '18.10.19/8rpm/FL','18.10.19/8rpm_2/FL','18.10.19/control/FL'...
    '22.10.19/40rpm/FL','22.10.19/40rpm_2/FL','22.10.19/control/FL'...
 %    '5.4.19/beads/40rpm/FL','5.4.19/beads/40rpm2/FL','5.4.19/beads/control/FL'...
%    '12.4.19/beads/40rpm/FL','12.4.19/beads/40rpm2/FL','12.4.19/beads/control/FL'...    
%  '22.3.19/beads/40rpm/FL','22.3.19/beads/40rpm2/FL','22.3.19/beads/control/FL'...
%    '29.3.19/beads/40rpm/FL','29.3.19/beads/40rpm2/FL','29.3.19/beads/control/FL'...
}

   

for insert=long;

%% load files
insert=insert{1};
%data_dir = strcat('/run/user/10704/gvfs/smb-share:server=sf3.bss.phy.private.cam.ac.uk,share=space/np451/ependymalJune/4.7.18/',insert);
data_dir= strcat('/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/Octorber2019/',insert)
a_folder=strcat('/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/Octorber2019/PIV/',insert);
%a_folder=strcat('/home/np451/Desktop/ependymal/4.7.18/PIV/');
mkdir(a_folder); 

cd(data_dir) 
suffix='.movie';
direc = dir('*FL*.movie');
%direc = dir('*BF*.movie');
N_files= size(direc,1);

%%

for i=1:N_files
    cd(data_dir)
    exp_name = direc(i).name;
    exp_name=exp_name(1:end-6);

%    if exist(strcat(a_folder,'/',exp_name)) == 0
        
%        mkdir(strcat(a_folder,'/',exp_name));
         
        mo=moviereader(direc(i).name);
        movie_path=strcat(data_dir,mo.Filename)
        frame_stack=mo.read();
        frame_stack=(frame_stack)-movmin(frame_stack,40,3);
        max_int= double(max(frame_stack(:)));
        if max_int>256        
        frame_stack= uint8((double(frame_stack)/max_int)*255);
        else frame_stack= uint8(frame_stack);
        end
%        frame_stack=imadjustn(frame_stack);
        for kk=1:size(frame_stack,3); 
            frame_stack(:,:,kk)= wiener2(frame_stack(:,:,kk),[5,5]);
        end
%        frame_stack= uint8(double(frame_stack)/2^(8));  %%%% converting images from 16 to 8 bit
%        frame_stack= uint8((frame_stack));  %%%% for images at 8 bit

%        cd(strcat(a_folder,'/',exp_name));


        %%  PIV

        [X,Y,U,V] = PIV_GetData(frame_stack,s,p);
        [x,y,u,v] =  PIV_ChangeFormat(X,Y,U,V);

%        save(strcat(exp_name,'PIV_result.mat'),'X','Y','U','V','x','y','u','v','p','s');
        close('all');
        
        cd(a_folder);

        %% PostProcessing 

        ulim= [-10,10];
        vlim= [-10,10];

        [U1,V1]= PIV_Validation(X,Y,U,V,ulim,vlim);
        [x,y,u1,v1] =  PIV_ChangeFormat(X,Y,U1,V1);

%        save('PIV_post.mat','U1','V1','ulim','vlim','u1','v1');
        %%%%%% save a figure of the first frame
        ss= std(double(frame_stack),[],3);
        figure(1);
        subplot(2,2,1)
        title(strcat(exp_name,'frame 1'));
        imagesc(ss);
      
        %%%%%% save the data before processing
        subplot(2,2,2)
        title(strcat(exp_name,' PreProcessing'));
        quiver(x,-y,nanmean(u,3),-nanmean(v,3),3);
  
        %%%% figure 
        subplot(2,2,3)
        title(strcat(exp_name,' PostProcessing'));
        quiver(x,-y,mean(u1,3),-mean(v1,3),3);
        fig=figure(1);
        saveas(fig,strcat(exp_name,'_PostPro.png'));
        close(1); 
  
        clear frame_stack; 
        save(strcat(exp_name,'.mat'));
%    end
end

end


