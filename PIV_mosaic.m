




%%%%% PIV Mosaic
ppm= 0.0014 ; %%%% mm  (using 4X objective 1 pixel is 1.4 um) 

pos=    {'CC','BC','B2','TC','T2','CB','CT','BT','TB','TT','BB'};
y_shift=floor([  0 ,-1.5,-3  ,1.5 , 3  , 0  , 0  ,-1.5, 1.5, 1.5,-1.5]/ppm); 
x_shift=-floor([  0 , 0  , 0  , 0  , 0  ,-1.5, 1.5, 1.5,-1.5, 1.5,-1.5]/ppm); 


path_fold = '/Users/nicola/Desktop/6Apr/CONTROL/PIV_64_64/'

cd(path_fold) 
suffix='*.movie';
direc = dir([path_fold,'4X_*']);
N_files= size(direc,1);

figure(1)
title('Flow pattern  PIV with validation')
xlabel('[mm]');
ylabel('[mm]');
xlim([-2,4.5])
ylim([-4.5,2])

hold on;

figure(2)
title('Flow pattern  PIV without validation')
xlabel('pixel');
ylabel('pixel');
xlim([-2,4.5])
ylim([-4.5,2])

hold on;





for i=1:N_files

    indx=0;    
    exp_name= direc(i).name;
    for j=1:size(pos,2);
        if (direc(i).isdir ==1) & (sum(pos{j} == exp_name(4:5)) == 2); 
        indx=j;
        end
    end
    disp(indx)
    if indx~=0   %%% indx corrispond to the
    cd(strcat(path_fold, exp_name));
    
    load('PIV_result.mat');
    load('PIV_post.mat');
    xn= (x-x_shift(indx))*ppm;
    yn= (y-y_shift(indx))*ppm;
    
    figure(1)
    quiver(xn,-yn,mean(u1,3),mean(-v1,3),2,'b');
    hold on;
    
    figure(2)
    quiver(xn,-yn,mean(u,3),mean(-v,3),2,'b');
    hold on;
    end
    
        
end




