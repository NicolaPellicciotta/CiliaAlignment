%%%%%% this script is suitable with PIV_script_maker
%%%% it read in order of date

%% setting area to explore
Lx=6000;    %%%% in um
Ly=6000;

obj_n= 10;  %%%% magnification
ppm= 5.84/obj_n;
obj= strcat(num2str(obj_n),'X_');

% setting frmae size
fsx= round(1200*ppm);  %%% all in micron
fsy= round(1200*ppm);

posx= (0:round((Lx/fsx)))* fsx;
posy= (0:round((Ly/fsx)))* fsy;

posx = (posx-round(Lx/2));   %%%% position in pixel unit
posy = (posy-round(Ly/2));


x_shift=[];y_shift=[];
for x=posx; for y=posy; x_shift=[x_shift,x]; y_shift=[y_shift,y]; end; end;


path_fold = '/home/np451/Desktop/Cilia_data/22.6/3_BF/'

cd(path_fold)

start=0;
for i=1:numel(x_shift)  %%%% load all coordinates and velocities in 4 arrays (xx,yy,uu,vv)
    

    cd(path_fold)
%%%1    file= dir(['*X',num2str(x_shift(i)),'Y',num2str(y_shift(i)),'*.movie']);
    exp_name= dir(['*X',num2str(x_shift(i)),'Y',num2str(y_shift(i)),'*']);

    if isempty(exp_name)==0;
    
        exp_name=exp_name.name;
    %     indx=0;    
    %     exp_name= direc(i).name;
    %     for j=1:size(pos,2);
    %         if (direc(i).isdir ==1) & (sum(pos{j} == exp_name(4:5)) == 2); 
    %         indx=j;
    %         end
    %     end
    %     disp(indx)
    %     if indx~=0   %%% indx corrispond to the

    %%%1    exp_name= file.name(1:end-6);
        cd(strcat(path_fold, exp_name));
        disp(i);


        load('PIV_result.mat');
        load('PIV_post.mat');
        xn= x+y_shift(i)/ppm;     %%%%% weird but it work seems rotation of 90 degree anticlockwise? (x,y)-> (y,-x)
        yn= y-x_shift(i)/ppm;
        u1m = nanmean(u1,3);
        v1m = nanmean(v1,3);
        um = nanmean(u,3);
        vm = nanmean(v,3);




        xx(start+1:start+numel(xn)) = xn;
        yy(start+1:start+numel(yn)) = yn;
        uu1(start+1:start+numel(u1m)) = u1m;
        vv1(start+1:start+numel(v1m)) = v1m;
        uu(start+1:start+numel(um)) = um;
        vv(start+1:start+numel(vm)) = vm;

        start= start+ numel(x(:));
    end
    
end
    
cd(path_fold);

%% SET new lattice with space dx and do a velocy average over a radius dx/2

dx=128;
px_max=(floor(max(xx(:))/dx)+1)*dx;  %%%% set the new axis limit
py_max=(floor(max(yy(:))/dx)+1)*dx;

XX= repmat((-px_max:dx:px_max),[numel((-py_max:dx:py_max)),1]);
YY= repmat((-py_max:dx:py_max)',[1,numel((-px_max:dx:px_max))]);
UU= zeros(size(XX));
VV= zeros(size(YY));
UU1= zeros(size(XX));
VV1= zeros(size(YY));



min_d=dx;

uu(isnan(uu))=0;
vv(isnan(vv))=0;
uu1(isnan(uu))=0;
vv1(isnan(vv))=0;



for i=1:size(XX,1)
    for j=1:size(XX,2)
        dr= sqrt( (xx(:)-XX(i,j)).^2 + (yy(:)-YY(i,j)).^2 ) ;
        UU1(i,j) = mean(uu1(dr(:)<min_d));
        VV1(i,j) = mean(vv1(dr(:)<min_d));
        UU(i,j) = mean(uu(dr(:)<min_d));
        VV(i,j) = mean(vv(dr(:)<min_d));
    
    
    end
end 


%% figures

close('all')
Arrow_l=100;

figure(1);
quiver(XX*ppm,-YY*ppm,UU,-VV,Arrow_l);
xlabel('X [um]');
ylabel(' Y[um]');
title('Flow pattern without validation.pdf');
axis image

%%% draw circle
r=3250;
xc= 0%;+fsx/2; 
yc= 0;%-+fsx/2; %%%% insert

%%% INLET (0,-2951)   ->(weird coordinates)  (-2951,0) 
r_inlet=500, x_inlet= xc- 2500  ; y_inlet=yc+0  %%%% inlet
r_outlet=500, x_outlet= xc +2500  ; y_outlet=yc +0  %%%% inlet
viscircles([xc,yc],r);
viscircles([x_inlet,y_inlet],r_inlet,'Color','k');
viscircles([x_outlet,y_outlet],r_outlet,'Color','y');

fig=figure(1);
saveas(fig,'Flow pattern without validation.pdf');



figure(2);
quiver(XX*ppm,-YY*ppm,UU1,-VV1,Arrow_l);
xlabel('[mm]');
ylabel('[mm]');
title('Flow pattern with validation');
axis image;
%%% draw circle

viscircles([xc,yc],r);
viscircles([x_inlet,y_inlet],r_inlet,'Color','k');
viscircles([x_outlet,y_outlet],r_outlet,'Color','y');

fig=figure(2);
saveas(fig,'Flow pattern with validation.pdf');


figure(3);
quiver(xx*ppm,-yy*ppm,uu,-vv,Arrow_l);
xlabel('X [um]');
ylabel(' Y[um]');
title('Flow pattern preAverage without validation');
axis image;
%%% draw circle

viscircles([xc,yc],r);
viscircles([x_inlet,y_inlet],r_inlet,'Color','k');
viscircles([x_outlet,y_outlet],r_outlet,'Color','y');

fig=figure(3);
saveas(fig,'Flow pattern preAverage without validation.pdf');

