% Read Siemens PC data

clc
list=dir;
oselist=list(~ismember({list.name},{'.','..','.DS_Store','._.DS_Store'}));
list = list([list(:).isdir]);
names={list.name};

idxMag=find(contains(names,'_7'));%highest
idxVx=find(contains(names,'_5'));%second lowest
idxVy=find(contains(names,'_6'));%second highest
idxVz=find(contains(names,'_4'));%lowest

[Vx,Svx] = dicom2struct_siemens(names{1,idxVx},'test');
cd ..
[Vy,Svy] = dicom2struct_siemens(names{1,idxVy},'test');
cd ..
[Vz,Svz] = dicom2struct_siemens(names{1,idxVz},'test');
cd ..
[M,Sm] = dicom2struct_siemens(names{1,idxMag},'test');

triggerX=unique([Svx.trigger]);
triggerY=unique([Svy.trigger]);
triggerZ=unique([Svz.trigger]);
pxSp=unique([Svx.pxSpacing,Svy.pxSpacing,Svz.pxSpacing]);

venc=15;
SS=-4096/2/venc*(-1);
Vx=double(Vx)/SS-4096/2/SS;
Vy=double(Vy)/SS-4096/2/SS;
Vz=double(Vz)/SS-4096/2/SS;

Mg=mat2gray(double(M));
mask=Mg;
mask(mask<0.1)=0;
mask(mask>=0.1)=1;

VxM=Vx.*mask;
VyM=Vy.*mask;
VzM=Vz.*mask;

%phase shading
VxCorrection=mean(VxM,4);
VyCorrection=mean(VyM,4);
VzCorrection=mean(VzM,4);

VxMc=VxM-VxCorrection;
VyMc=VyM-VyCorrection;
VzMc=VzM-VzCorrection;

mask(mask==0)=NaN;
VxMc=VxMc.*mask;
VyMc=VyMc.*mask;
VzMc=VzMc.*mask;

% %% HISTOGRAM---------------------------------
% 
% VxMh=VxM;
% VxMh(VxM==0)=NaN;
% VyMh=VzM;
% VyMh(VyM==0)=NaN;
% VzMh=VzM;
% VzMh(VzM==0)=NaN;
% 
% % figure
% % histogram(VxMh,200,'EdgeColor','none'); hold on
% % histogram(VyMh,200,'EdgeColor','none')
% % histogram(VzMh,200,'EdgeColor','none')
% 
% VxMch=VxMc;
% VxMch(VxMc==0)=NaN;
% VyMch=VyMc;
% VyMch(VyMc==0)=NaN;
% VzMch=VzMc;
% VzMch(VzMc==0)=NaN;
% 
% % figure
% % histogram(VxMch,200,'EdgeColor','none'); hold on
% % histogram(VyMch,200,'EdgeColor','none')
% % histogram(VzMch,200,'EdgeColor','none')
% % %set(gca,'yscale','log')
% % legend({'$v_x$','$v_y$','$v_z$'},'Interpreter','latex','FontSize',18,'location','best');
% 
% ax = gca;
% ax.TickLabelInterpreter = 'latex';
% ax.FontSize = 18;
% xlabel('$\mathrm{velocity} \quad \left[ \mathrm{cm/s} \right] $', 'Interpreter','latex','FontSize',18);

%%
save('siemensProcessed.mat','VxMc','VyMc','VzMc','triggerX','triggerY','triggerZ','mask','Mg','pxSp')


%%

% Recon.seriesNum
% Recon.ID
% Recon.location
% Recon.MVC
% Recon.data
% Recon.M
% Recon.Vx
% Recon.Vy
% Recon.Vz




%%
num_iter=10; kappa=1; option=1; delta_t=3/44;
pixelspacing = [pxSp,pxSp,5.0000]';

series=1;
volume = 1;

VxMc(isnan(VxMc))=0;
VyMc(isnan(VyMc))=0;
VzMc(isnan(VzMc))=0;

Data(volume).vx = VxMc;
Data(volume).vy = VyMc;
Data(volume).vz = VzMc;
Data(volume).mask = mask;

Data(volume).m  = Mg;
frames=size(Data(volume).vx,4);
    
    
    h=waitbar(0,'Filtering');
    
    for k=1:frames

         Data(volume).vx_sm(:,:,:,k) = anisodiff3D(Data(volume).vx(:,:,:,k),...
                              num_iter,delta_t,kappa,option,pixelspacing);
         Data(volume).vy_sm(:,:,:,k) = anisodiff3D(Data(volume).vy(:,:,:,k),...
                              num_iter,delta_t,kappa,option,pixelspacing);
         Data(volume).vz_sm(:,:,:,k) = anisodiff3D(Data(volume).vz(:,:,:,k),...
                              num_iter,delta_t,kappa,option,pixelspacing);

         waitbar(k/frames,h)
         
    end
    
    close(h)
    
    
    
%%

slice=15;

a=figure;
imshow(mat2gray(Data(end).m(:,:,slice,1)),'Initialmagnification',400)
h_all = drawfreehand(gca);
mask_all = roiWait(h_all);
close(a)

a=figure;
imshow(mat2gray(Data(end).m(:,:,slice,1)),'Initialmagnification',400);
h_gm  = drawrectangle(gca,'Position',[50 50 5 15]);
mask_gm = roiWait(h_gm);
saveas(a,'roi_1_slice-8-ST.jpg')
close(a)

a=figure;
imshow(mat2gray(Data(end).m(:,:,slice,1)),'Initialmagnification',400);
h_sol  = drawrectangle(gca,'Position',[50 50 5 15]);
mask_sol = roiWait(h_sol);
close(a)



%%
trigger=triggerX-cat(2,0,triggerX(1:end-1));
trigger(1)=[];
trigger=trigger/1000;

mask=Data(series).mask(:,:,:,:);
mask(mask==0)=NaN;
[x,y,z]=meshgrid(1:size(Mg,2),1:size(Mg,1),1:size(Mg,3));

%                                                                                                                                                                                                                                 

x(isnan(mask(:,:,:,1)))=NaN;
y(isnan(mask(:,:,:,1)))=NaN;
z(isnan(mask(:,:,:,1)))=NaN;


XS=zeros([size(x),frames]);
YS=XS;
ZS=XS;
VX=XS;
VY=XS;
VZ=XS;
VR=XS;

    for yi=1:size(Mg,1)
            yi
          for xi=1:size(Mg,2)
              
              for zi=14:18%size(Mg,3)
              
                [xs,ys,zs,vx,vy,vz,vr]=track3dv2(x(yi,xi,zi),y(yi,xi,zi),z(yi,xi,zi),-Data(series).vx_sm,-Data(series).vy_sm,Data(series).vz_sm,trigger,pixelspacing);
          
                XS(yi,xi,zi,:)=xs;
                YS(yi,xi,zi,:)=ys;
                ZS(yi,xi,zi,:)=zs;
                VX(yi,xi,zi,:)=vx;
                VY(yi,xi,zi,:)=vy;
                VZ(yi,xi,zi,:)=vz;
                VR(yi,xi,zi,:)=vr;
                
              end
          
          end
    end


save('siemens_tracked.mat')    
    
%%

mask_sol_dynamic=zeros(size(XS));
mask_sol_dynamic(:,:,14:18,:)=repmat(mask_sol,[1,1,5,frames]);

mask_gm_dynamic=zeros(size(XS));
mask_gm_dynamic(:,:,14:18,:)=repmat(mask_gm,[1,1,5,frames]);

       %%                         
                         
% tensor calcs                                
Data(series).strain     =   strain_tensor_calcs_v2(VX,Data(series).vx_sm, VY, Data(series).vy_sm, VZ, Data(series).vz_sm,XS,YS,ZS,pixelspacing,trigger);
%%
% masking
Data(series).strain_SOL = strain_ROI_v3Siemens(Data(series).strain,double(mask_sol_dynamic));
Data(series).strain_GM  = strain_ROI_v3Siemens(Data(series).strain,double(mask_gm_dynamic));

% peak
Data(series).strain_SOL_peak = peak_ROI_SiemensV2(Data(series).strain_SOL);
Data(series).strain_GM_peak  = peak_ROI_SiemensV2(Data(series).strain_GM);

%colormaps
foldername='cmap';
mkdir(foldername)
cd(foldername)

for slice=15

%m=repmat(Data(series).m(:,:,slice,1),[1,1,1,size(Data(series).m,4)]);
m=repmat(Data(series).m(:,:,slice,:),[1,1,1,1]);


% displacements

limits=[0,7];
suffix=[foldername,num2str(slice),'_dx'];
units = '$|\Delta_{x}| \quad  [\mathrm{mm}]$';
d=10*squeeze(Data(series).strain.dx(:,:,slice,:));
colormap_strain(m,abs(d).*squeeze(mask(:,:,slice,:)),suffix,limits,units);
vepcGIFcmaps(m,abs(d).*squeeze(mask(:,:,slice,:)),suffix,limits,units);

suffix=[foldername,num2str(slice),'_dy'];
units = '$|\Delta_{y}| \quad  [\mathrm{mm}]$';
d=10*squeeze(Data(series).strain.dy(:,:,slice,:));
colormap_strain(m,abs(d).*squeeze(mask(:,:,slice,:)),suffix,limits,units);
vepcGIFcmaps(m,abs(d).*squeeze(mask(:,:,slice,:)),suffix,limits,units);

suffix=[foldername,num2str(slice),'_dz'];
units = '$|\Delta_{z}| \quad  [\mathrm{mm}]$';
d=10*squeeze(Data(series).strain.dz(:,:,slice,:));
colormap_strain(m,abs(d).*squeeze(mask(:,:,slice,:)),suffix,limits,units);
vepcGIFcmaps(m,abs(d).*squeeze(mask(:,:,slice,:)),suffix,limits,units);


% velocities
limits=[-3,3];
suffix=[foldername,num2str(slice),'_vx'];
units = '$v_{x} \quad  [\mathrm{cm/s}]$';
d=squeeze(Data(series).strain.vx(:,:,slice,:));
colormap_strain(m,d.*squeeze(mask(:,:,slice,:)),suffix,limits,units);
vepcGIFcmaps(m,d.*squeeze(mask(:,:,slice,:)),suffix,limits,units);

suffix=[foldername,num2str(slice),'_vy'];
units = '$v_{y} \quad  [\mathrm{cm/s}]$';
d=squeeze(Data(series).strain.vy(:,:,slice,:));
colormap_strain(m,d.*squeeze(mask(:,:,slice,:)),suffix,limits,units);
vepcGIFcmaps(m,d.*squeeze(mask(:,:,slice,:)),suffix,limits,units);

suffix=[foldername,num2str(slice),'_vz'];
units = '$v_{z} \quad  [\mathrm{cm/s}]$';
d=squeeze(Data(series).strain.vz(:,:,slice,:));
colormap_strain(m,d.*squeeze(mask(:,:,slice,:)),suffix,limits,units);
vepcGIFcmaps(m,d.*squeeze(mask(:,:,slice,:)),suffix,limits,units);


% Strain Euler
limits=[-0.7,0.7];
suffix=[foldername,num2str(slice),'_Eneg'];
units = '$E_{\lambda_{1}} \quad  [\mathrm{mm/mm}]$';
d=squeeze(Data(series).strain.E_lambda(:,:,slice,:,1));
colormap_strain(m,d.*squeeze(mask(:,:,slice,:)),suffix,limits,units);
vepcGIFcmaps(m,d.*squeeze(mask(:,:,slice,:)),suffix,limits,units);

suffix=[foldername,num2str(slice),'_Esum'];
units = '$E_{\lambda_{2}} \quad  [\mathrm{mm/mm}]$';
d=squeeze(Data(series).strain.E_lambda(:,:,slice,:,2));
colormap_strain(m,d.*squeeze(mask(:,:,slice,:)),suffix,limits,units);
vepcGIFcmaps(m,d.*squeeze(mask(:,:,slice,:)),suffix,limits,units);

suffix=[foldername,num2str(slice),'_Epos'];
units = '$E_{\lambda_{3}} \quad  [\mathrm{mm/mm}]$';
d=squeeze(Data(series).strain.E_lambda(:,:,slice,:,3));
colormap_strain(m,d.*squeeze(mask(:,:,slice,:)),suffix,limits,units);
vepcGIFcmaps(m,d.*squeeze(mask(:,:,slice,:)),suffix,limits,units);

suffix=[foldername,num2str(slice),'_Eshear'];
units = '$E_{max} \quad  [\mathrm{mm/mm}]$';
d=squeeze(Data(series).strain.ShearE_max(:,:,slice,:));
colormap_strain(m,d.*squeeze(mask(:,:,slice,:)),suffix,limits,units);
vepcGIFcmaps(m,d.*squeeze(mask(:,:,slice,:)),suffix,limits,units);

suffix=[foldername,num2str(slice),'_Evol'];
units = '$E_{vol}=\frac{\delta V}{V} \quad  [\mathrm{mm^3/mm^3}]$';
d=squeeze(Data(series).strain.E_Volumetric(:,:,slice,:));
colormap_strain(m,d.*squeeze(mask(:,:,slice,:)),suffix,limits,units);
vepcGIFcmaps(m,d.*squeeze(mask(:,:,slice,:)),suffix,limits,units);

% Strain Lagrange
limits=[-0.7,0.7];
suffix=[foldername,num2str(slice),'_Lneg'];
units = '$L_{\lambda_{1}} \quad  [\mathrm{mm/mm}]$';
d=squeeze(Data(series).strain.L_lambda(:,:,slice,:,1));
colormap_strain(m,d.*squeeze(mask(:,:,slice,:)),suffix,limits,units);
vepcGIFcmaps(m,d.*squeeze(mask(:,:,slice,:)),suffix,limits,units);

suffix=[foldername,num2str(slice),'_Lsum'];
units = '$L_{\lambda_{2}} \quad  [\mathrm{mm/mm}]$';
d=squeeze(Data(series).strain.L_lambda(:,:,slice,:,2));
colormap_strain(m,d.*squeeze(mask(:,:,slice,:)),suffix,limits,units);
vepcGIFcmaps(m,d.*squeeze(mask(:,:,slice,:)),suffix,limits,units);

suffix=[foldername,num2str(slice),'_Lpos'];
units = '$L_{\lambda_{3}} \quad  [\mathrm{mm/mm}]$';
d=squeeze(Data(series).strain.L_lambda(:,:,slice,:,3));
colormap_strain(m,d.*squeeze(mask(:,:,slice,:)),suffix,limits,units);
vepcGIFcmaps(m,d.*squeeze(mask(:,:,slice,:)),suffix,limits,units);

suffix=[foldername,num2str(slice),'_Lshear'];
units = '$L_{max} \quad  [\mathrm{mm/mm}]$';
d=squeeze(Data(series).strain.ShearL_max(:,:,slice,:));
colormap_strain(m,d.*squeeze(mask(:,:,slice,:)),suffix,limits,units);
vepcGIFcmaps(m,d.*squeeze(mask(:,:,slice,:)),suffix,limits,units);

suffix=[foldername,num2str(slice),'_Lvol'];
units = '$L_{vol}=\frac{\delta V}{V} \quad  [\mathrm{mm^3/mm^3}]$';
d=squeeze(Data(series).strain.L_Volumetric(:,:,slice,:));
colormap_strain(m,d.*squeeze(mask(:,:,slice,:)),suffix,limits,units);
vepcGIFcmaps(m,d.*squeeze(mask(:,:,slice,:)),suffix,limits,units);


% Strain Rate
limits=[-1500,1500];
suffix=[foldername,num2str(slice),'_SRneg'];
units = '$SR_{\lambda_{1}} \quad  [\mathrm{s^{-1}}]$';
d=squeeze(Data(series).strain.SR_lambda(:,:,slice,:,1));
colormap_strain(m,d,suffix,limits,units);
vepcGIFcmaps(m,d.*squeeze(mask(:,:,slice,:)),suffix,limits,units);

suffix=[foldername,num2str(slice),'_SRsum'];
units = '$SR_{\lambda_{2}} \quad  [\mathrm{s^{-1}}]$';
d=squeeze(Data(series).strain.SR_lambda(:,:,slice,:,2));
colormap_strain(m,d,suffix,limits,units);
vepcGIFcmaps(m,d.*squeeze(mask(:,:,slice,:)),suffix,limits,units);

suffix=[foldername,num2str(slice),'_SRpos'];
units = '$SR_{\lambda_{3}} \quad  [\mathrm{s^{-1}}]$';
d=squeeze(Data(series).strain.SR_lambda(:,:,slice,:,3));
colormap_strain(m,d,suffix,limits,units);
vepcGIFcmaps(m,d.*squeeze(mask(:,:,slice,:)),suffix,limits,units);

suffix=[foldername,num2str(slice),'_SRshear'];
units = '$SR_{max} \quad  [\mathrm{s^{-1}}]$';
d=squeeze(Data(series).strain.ShearSR_max(:,:,slice,:));
colormap_strain(m,d,suffix,limits,units);
vepcGIFcmaps(m,d.*squeeze(mask(:,:,slice,:)),suffix,limits,units);

% FA COLORMAPS 
limits = [0, 1]; % FA values range from 0 to 1
suffix = [foldername, num2str(slice), '_FAE']; % Name for FA_E
units = 'FAE'; % FA is unitless
d = squeeze(Data.strain.FAE(:,:,slice,:)); % Get FA_E data from the struct
colormap_strain(m, abs(d).*squeeze(mask(:,:,slice,:)), suffix, limits, units); % Display FA_E
vepcGIFcmaps(m, abs(d).*squeeze(mask(:,:,slice,:)), suffix, limits, units); % Create GIF for FA_E

limits = [0, 1]; % FA values range from 0 to 1
suffix = [foldername, num2str(slice), '_FAL']; % Name for FA_L
units = 'FAL'; % FA is unitless
d = squeeze(Data.strain.FAL(:,:,slice,:)); % Get FA_L data from the struct
colormap_strain(m, abs(d).*squeeze(mask(:,:,slice,:)), suffix, limits, units); % Display FA_l
vepcGIFcmaps(m, abs(d).*squeeze(mask(:,:,slice,:)), suffix, limits, units); % Create GIF for FA_l

limits = [0, 1]; % FA values range from 0 to 1
suffix = [foldername, num2str(slice), '_FASR']; % Name for FA_SR
units = 'FASR'; % FA is unitless
d = squeeze(Data.strain.FASR(:,:,slice,:)); % Get FA_SR data from the struct
colormap_strain(m, abs(d).*squeeze(mask(:,:,slice,:)), suffix, limits, units); % Display FA_SR
vepcGIFcmaps(m, abs(d).*squeeze(mask(:,:,slice,:)), suffix, limits, units); % Create GIF for FA_SR

end
cd ..






%%
strain_mvc_plots_v2Siemens(Data.strain_GM,'GM',trigger);
strain_mvc_plots_v2Siemens(Data.strain_SOL,'SOL',trigger);

save('siemens_final.mat','-v7.3')



