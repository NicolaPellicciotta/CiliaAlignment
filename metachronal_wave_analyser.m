path_fold = '/home/np451/Desktop/Cilia_data/17.6/1_BF/'

cd(path_fold);
d=dir('*X505Y-1598*.movie');
filename=d.name;
mo=moviereader(filename);
frame_stack=mo.read();
Iqtau=DDM_core(frame_stack,200);
A= mean(Iqtau,2);
plot(A)

%% DDM asymmetryc
cilia.height = 1200;
cilia.width = 1200;
cilia.max_tau = floor(size(frame_stack,3)/2);

fs = single(frame_stack(1:cilia.height, 1:cilia.width,:));

Iqtau2 = zeros(cilia.height,cilia.width,cilia.max_tau,'single');

for tau = 1:cilia.max_tau
    
    Iqtau2(:,:,tau) =  ...
        mean(abs( 1/(cilia.height*cilia.width).* fft2(fs(:,:,1:end-tau) - fs(:,:,tau+1:end))).^2,3); % not shifted
    
end %for

Iqtau2(1,1,:) = 0; %to avoid very high centre number (mode 0)

Iq = fftshift(mean(Iqtau2,3)); % average over tau and shiftsave





%% find maxima

MIq = max(Iq(:));

[my,mx] = find(Iq == MIq);
cc = cilia.width/2+1;

lambda = cilia.height ./ hypot(my-cc,mx-cc)  % in pixels ofc;
prop_angle_deg = atan2d(my-cc,mx-cc)        % comes out negative because of reverse ydir in imagesc;

hf = figure(20);
clf
hf.Units = 'centimeters';
hf.Position([3 4]) = [8 6];
ha = axes('Position',[0.13 0.15 0.83 0.87] );
hi = imagesc([-10:10],[-10:10],...
    Iq(cc+[-10:10],cc+[-10:10])./10);
axis image
hc = colorbar;
hxl = xlabel('mode_x','FontSize',12);
hxl.Units = 'Normalized';
hxl.Position(2) = -0.1;
hyl = ylabel('mode_y','FontSize',12);
hyl.Units = 'Normalized';
hyl.Position(1) = -0.09;
hapos = ha.Position;
hc.Position(1) = 0.73;
drawnow
ha.Position = hapos;
hc.Label.String = '<I(q,\tau)>_\tau, [a.u.]';
hc.Label.FontSize = 12;
hc.Label.Units = 'normalized';
hc.Label.Position(1) = 1.6;
hold on
plot(0,0,'+w','MarkerSize',8)
plot(mx-cc,my-cc,'+r','MarkerSize',8)
% print2svg(hf,'Iq.svg')