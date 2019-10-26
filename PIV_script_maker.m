fid = fopen('BF_0.5x0.5_40X.xml','wt');
fprintf(fid,'<temika>'); fprintf(fid,'\n');

%% open temika
fprintf(fid,'\t<microscope>'); fprintf(fid,'\n');
fprintf(fid,'\t\t<eclipsetie>'); fprintf(fid,'\n');
fprintf(fid,'\t\t\t<light_path>L100</light_path>'); fprintf(fid,'\n');
fprintf(fid,'\t\t</eclipsetie>'); fprintf(fid,'\n');
fprintf(fid,'\t</microscope>'); fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\t\t<timestamp>0</timestamp>'); fprintf(fid,'\n');
fprintf(fid,'\n');fprintf(fid,'\n');fprintf(fid,'\n');

%% setting area to explore
%%%%%% it work with field of you 1200 px

Lx=500;
Ly=500;

obj_n= 40;  %%%% magnification
ppm= 5.84/obj_n;
obj= strcat(num2str(obj_n),'X_');

% setting frmae size
fsx= round(1200*ppm);
fsy= round(1200*ppm);

posx= (0:round((Lx/fsx)))* fsx;
posy= (0:round((Ly/fsx)))* fsy;

posx = posx-round(Lx/2);
posy = posy-round(Ly/2);


pathdir='/home/np451/data/'
pause=3.0;   %%%% seconds int
record=20.0; %%%% seconds of movie 

%% For on moving in that position, taking 10s video and wait 20s 
%%for manual focus refinement

for x=posx;
    for y=posy
        
        filename=strcat(obj,'X',num2str(x),'Y',num2str(y));
        
        fprintf(fid,strcat('\t<!-- Position x=',num2str(x),'y=',num2str(y),' -->'));fprintf(fid,'\n');
    
        fprintf(fid,'\t\t<!-- MOVE -->');fprintf(fid,'\n');
        fprintf(fid,'\t\t<microscope>');fprintf(fid,'\n');
 
        fprintf(fid,'\t\t\t<xystage axis="x">');fprintf(fid,'\n');
        fprintf(fid,strcat('\t\t\t\t<move_absolute>',num2str(x),' 5</move_absolute>'));fprintf(fid,'\n');
        fprintf(fid,'\t\t\t</xystage>');fprintf(fid,'\n');

        fprintf(fid,'\t\t\t<xystage axis="y">');fprintf(fid,'\n');
        fprintf(fid,strcat('\t\t\t\t<move_absolute>',num2str(y),' 5</move_absolute>'));fprintf(fid,'\n');
        fprintf(fid,'\t\t\t</xystage>');fprintf(fid,'\n');
        
        
        fprintf(fid,'\t\t\t<xystage axis="x">');fprintf(fid,'\n');
        fprintf(fid,'\t\t\t\t<wait_moving_end></wait_moving_end>');fprintf(fid,'\n');
        fprintf(fid,'\t\t\t</xystage>');fprintf(fid,'\n');
        
        
        fprintf(fid,'\t\t\t<xystage axis="y">');fprintf(fid,'\n');
        fprintf(fid,'\t\t\t\t<wait_moving_end></wait_moving_end>');fprintf(fid,'\n');
        fprintf(fid,'\t\t\t</xystage>');fprintf(fid,'\n');
        fprintf(fid,'\t\t</microscope>');fprintf(fid,'\n');
        
        
        fprintf(fid,'\t\t<!-- SET NAME -->');fprintf(fid,'\n');
        fprintf(fid,'\t\t<save>');fprintf(fid,'\n');
        fprintf(fid,strcat('\t\t\t<basename>',pathdir,filename,'</basename>'));fprintf(fid,'\n');
        fprintf(fid,'\t\t\t<append>DATE</append>');fprintf(fid,'\n');
        fprintf(fid,'\t\t</save>');fprintf(fid,'\n');
        
        fprintf(fid,'\t\t<!-- Turn on the light -->');fprintf(fid,'\n');
        fprintf(fid,strcat('\t\t<sleep>0:0:',num2str(pause),'</sleep>'));fprintf(fid,'\n');
        fprintf(fid,'\t\t<!-- Take a record -->');fprintf(fid,'\n');
        fprintf(fid,'\t\t<camera name="IIDC Point Grey Research Grasshopper3 GS3-U3-23S6M"><record>ON</record></camera>');fprintf(fid,'\n');
        fprintf(fid,strcat('\t\t<sleep>0:0:',num2str(record),'</sleep>'));fprintf(fid,'\n');
        fprintf(fid,'\t\t<camera name="IIDC Point Grey Research Grasshopper3 GS3-U3-23S6M"><record>OFF</record></camera>');fprintf(fid,'\n');
        fprintf(fid,'\t\t<sleep>0:0:1.0</sleep>');fprintf(fid,'\n');
        
        
        
        
        
        
        
        
        fprintf(fid,'\n');fprintf(fid,'\n');fprintf(fid,'\n');fprintf(fid,'\n');
        
    end
end
        
        
fprintf(fid,'</temika>'); fprintf(fid,'\n');
fclose(fid);