function status=scat_hic_interpolate()
% need to pad/cover the boundaries, interpolate and clean up.

% VARIABLES TO BE DEFINED
load('boundary_fields_Si3N4_v3_d_integrated')

dxy=[0.020,0.020];
dxy_nm=[20,20];

start_p_a_nm=220;  % do it manually to be sure to have integers and not messed up floats
end_p_a_nm=1000;
start_p_b_nm=200;
end_p_b_nm=1000;
%%%% END OF VARIABLES TO BE DEFINED MANUALLY

close all

nm=1e-9;
um=1e-6;

[mesh_p_b,mesh_p_a]=meshgrid(start_p_b_nm:dxy_nm(1):end_p_b_nm,start_p_a_nm:dxy_nm(2):end_p_a_nm);
mesh_p_a=mesh_p_a*nm; %for power
mesh_p_b=mesh_p_b*nm;

mesh_PowerHIC_TE_total=zeros(size(mesh_p_a)); % mesh offers rectangular data structure instead of the triangular data structures used in the mat file
mesh_PowerHIC_TM_total=zeros(size(mesh_p_a));
mesh_PowerHIC_TE_total_r=zeros(size(mesh_p_a)); % "r" for reference, used to clean up later
mesh_PowerHIC_TM_total_r=zeros(size(mesh_p_a)); % "r" for reference, used to clean up later
mesh_PowerHIC_TE_factor=zeros(size(mesh_p_a));
mesh_PowerHIC_TM_factor=zeros(size(mesh_p_a));
mesh_PowerHIC_TE_factor_r=zeros(size(mesh_p_a)); % "r" for reference, used to clean up later
mesh_PowerHIC_TM_factor_r=zeros(size(mesh_p_a)); % "r" for reference, used to clean up later
for ia=1:length(data(:,1))
    ib=1;
    exit_ib=0;
    first_row_value=0;
    while exit_ib==0
        % just changing name (utility removed)
        PowerHIC_TE_total_c(ia,ib)=PowerHIC_TE_total(ia,ib);
        PowerHIC_TM_total_c(ia,ib)=PowerHIC_TM_total(ia,ib);
        PowerHIC_TE_factor_c(ia,ib)=PowerHIC_TE_factor(ia,ib);
        PowerHIC_TM_factor_c(ia,ib)=PowerHIC_TM_factor(ia,ib);    
        
        % placing data in the mesh rectangular structure
        mpia=round((matrix_a(ia,ib)-start_p_a_nm*nm)/(dxy_nm(2)*nm))+1;
        mpib=round((matrix_b(ia,ib)-start_p_b_nm*nm)/(dxy_nm(1)*nm))+1;  
        mesh_PowerHIC_TE_total(mpia,mpib)=PowerHIC_TE_total_c(ia,ib);
        mesh_PowerHIC_TM_total(mpia,mpib)=PowerHIC_TM_total_c(ia,ib);
        mesh_PowerHIC_TE_factor(mpia,mpib)=PowerHIC_TE_factor_c(ia,ib);
        mesh_PowerHIC_TM_factor(mpia,mpib)=PowerHIC_TM_factor_c(ia,ib);
        
        % keeping one mesh as reference for some plots later
        mesh_PowerHIC_TE_total_r(mpia,mpib)=PowerHIC_TE_total_c(ia,ib);
        mesh_PowerHIC_TM_total_r(mpia,mpib)=PowerHIC_TM_total_c(ia,ib);
        mesh_PowerHIC_TE_factor_r(mpia,mpib)=PowerHIC_TE_factor_c(ia,ib);
        mesh_PowerHIC_TM_factor_r(mpia,mpib)=PowerHIC_TM_factor_c(ia,ib);
        
        %cover from beginning if this is first value (horizontal padding)
        if first_row_value==0
            first_row_value=1;
            mesh_PowerHIC_TE_total(mpia,1:mpib)=mesh_PowerHIC_TE_total(mpia,mpib);
            mesh_PowerHIC_TM_total(mpia,1:mpib)=mesh_PowerHIC_TM_total(mpia,mpib);
    
            mesh_PowerHIC_TE_factor(mpia,1:mpib)=mesh_PowerHIC_TE_factor(mpia,mpib);
            mesh_PowerHIC_TM_factor(mpia,1:mpib)=mesh_PowerHIC_TM_factor(mpia,mpib);
            
            %compile first value
            first_value_power(mpia)=mpib;
        end

        ib=ib+1;
        if (ib<=length(data(1,:)))
            if (data(ia,ib).a>0)
            else
                exit_ib=1;
            end
        else
            exit_ib=1;
        end
        if exit_ib==1
            %compile last values for future cleanup
            last_value_power(mpia)=mpib;
        end
    end
end
% cover vertically (some more padding, on the right of the data)
for mpib=1:length(mesh_p_a(1,:))
    mpia=max(find(mesh_PowerHIC_TE_total(:,mpib)>0));
    if mpia>=1
        mesh_PowerHIC_TE_total(mpia:length(mesh_p_a(:,1)),mpib)=mesh_PowerHIC_TE_total(mpia,mpib);
        mesh_PowerHIC_TM_total(mpia:length(mesh_p_a(:,1)),mpib)=mesh_PowerHIC_TM_total(mpia,mpib);  
        mesh_PowerHIC_TE_factor(mpia:length(mesh_p_a(:,1)),mpib)=mesh_PowerHIC_TE_factor(mpia,mpib);
        mesh_PowerHIC_TM_factor(mpia:length(mesh_p_a(:,1)),mpib)=mesh_PowerHIC_TM_factor(mpia,mpib);      
    end
end 
% change all the zeros into NaNs in the reference mesh
mesh_PowerHIC_TE_total_r(find(mesh_PowerHIC_TE_total_r==0))=NaN;
mesh_PowerHIC_TM_total_r(find(mesh_PowerHIC_TM_total_r==0))=NaN;
mesh_PowerHIC_TE_factor_r(find(mesh_PowerHIC_TE_factor_r==0))=NaN;
mesh_PowerHIC_TM_factor_r(find(mesh_PowerHIC_TM_factor_r==0))=NaN;

%Interpolate!
[mesh_b_i,mesh_a_i] = meshgrid(start_p_b_nm:dxy_nm(1)/4:end_p_b_nm,start_p_a_nm:dxy_nm(2)/4:end_p_a_nm); % the interpolated mesh is 4 times finer
mesh_b_i=mesh_b_i*nm;
mesh_a_i=mesh_a_i*nm;
mesh_PowerHIC_TE_total_i=interp2(mesh_p_b,mesh_p_a,mesh_PowerHIC_TE_total,mesh_b_i,mesh_a_i,'spline');
mesh_PowerHIC_TM_total_i=interp2(mesh_p_b,mesh_p_a,mesh_PowerHIC_TM_total,mesh_b_i,mesh_a_i,'spline');
mesh_PowerHIC_TE_factor_i=interp2(mesh_p_b,mesh_p_a,mesh_PowerHIC_TE_factor,mesh_b_i,mesh_a_i,'spline');
mesh_PowerHIC_TM_factor_i=interp2(mesh_p_b,mesh_p_a,mesh_PowerHIC_TM_factor,mesh_b_i,mesh_a_i,'spline');

%Clean up!
remove_pixels=0; %%% how much pixels at the data edge should be removed
for mai=1:length(mesh_a_i(:,1))
    for mbi=1:length(mesh_a_i(1,:))
        %find the corresponding coarse mesh index
        mpia=max(find(mesh_p_a(:,1)<=mesh_a_i(mai,1)));
        if (mesh_b_i(mai,mbi)<mesh_p_b(mpia,first_value_power(mpia)))||((mbi+remove_pixels)>=length(mesh_a_i(1,:)))||(mesh_b_i(mai,mbi+remove_pixels)>mesh_p_b(mpia,last_value_power(mpia)))
            mesh_PowerHIC_TE_total_i(mai,mbi)=NaN;
            mesh_PowerHIC_TM_total_i(mai,mbi)=NaN;
            mesh_PowerHIC_TE_factor_i(mai,mbi)=NaN;
            mesh_PowerHIC_TM_factor_i(mai,mbi)=NaN;    
        end
    end
end

save(mfilename)


% plot results
figure(3)

subplot(3,2,1)
surf(mesh_p_b,mesh_p_a,mesh_PowerHIC_TE_total_r);
title('Power TE Total')
ylabel('Height of the core')
xlabel('Width of the core')
view([0 0 1]);
axis equal

subplot(3,2,2)
surf(mesh_p_b,mesh_p_a,mesh_PowerHIC_TM_total_r);
title('Power TM Total')
ylabel('Height of the core')
xlabel('Width of the core')
view([0 0 1]);
axis equal

subplot(3,2,3)
surf(mesh_p_b,mesh_p_a,mesh_PowerHIC_TE_total);
title('Power TE Total')
ylabel('Height of the core')
xlabel('Width of the core')
view([0 0 1]);
axis equal

subplot(3,2,4)
surf(mesh_p_b,mesh_p_a,mesh_PowerHIC_TM_total);
title('Power TM Total')
ylabel('Height of the core')
xlabel('Width of the core')
view([0 0 1]);
axis equal

subplot(3,2,5)
surf(mesh_b_i,mesh_a_i,mesh_PowerHIC_TE_total_i);
title('Power TE total')
ylabel('Height of the core')
xlabel('Width of the core')
view([0 0 1]);
axis equal

subplot(3,2,6)
surf(mesh_b_i,mesh_a_i,mesh_PowerHIC_TM_total_i);
title('Power TM Total')
ylabel('Height of the core')
xlabel('Width of the core')
view([0 0 1]);
axis equal

figure(4)

subplot(3,2,1)
contourf(mesh_p_b,mesh_p_a,mesh_PowerHIC_TE_total_r);
title('Normalized overlap at the sidewalls, TE')
ylabel('Height of the core')
xlabel('Width of the core')
view([0 0 1]);
axis equal

subplot(3,2,2)
contourf(mesh_p_b,mesh_p_a,mesh_PowerHIC_TM_total_r);
title('Normalized overlap at the sidewalls, TM')
ylabel('Height of the core')
xlabel('Width of the core')
view([0 0 1]);
axis equal

subplot(3,2,3)
contourf(mesh_p_b,mesh_p_a,mesh_PowerHIC_TE_total);
title('Normalized overlap at the sidewalls, TE')
ylabel('Height of the core')
xlabel('Width of the core')
view([0 0 1]);
axis equal

subplot(3,2,4)
contourf(mesh_p_b,mesh_p_a,mesh_PowerHIC_TM_total);
title('Normalized overlap at the sidewalls, TM')
ylabel('Height of the core')
xlabel('Width of the core')
view([0 0 1]);
axis equal

subplot(3,2,5)
[C,h]=contourf(mesh_b_i,mesh_a_i,mesh_PowerHIC_TE_total_i);
%h = clabel(C,h);
title('Normalized overlap at the sidewalls, TE')
ylabel('Height of the core')
xlabel('Width of the core')
view([0 0 1]);
axis equal

subplot(3,2,6)
[C,h]=contourf(mesh_b_i,mesh_a_i,mesh_PowerHIC_TM_total_i);
%h = clabel(C,h);
title('Normalized overlap at the sidewalls, TM')
ylabel('Height of the core')
xlabel('Width of the core')
view([0 0 1]);
axis equal

figure(5)

subplot(3,2,1)
surf(mesh_p_b,mesh_p_a,mesh_PowerHIC_TE_factor_r);
title('Power TE factor')
ylabel('Height of the core')
xlabel('Width of the core')
view([0 0 1]);
axis equal

subplot(3,2,2)
surf(mesh_p_b,mesh_p_a,mesh_PowerHIC_TM_factor_r);
title('Power TM factor')
ylabel('Height of the core')
xlabel('Width of the core')
view([0 0 1]);
axis equal

subplot(3,2,3)
surf(mesh_p_b,mesh_p_a,mesh_PowerHIC_TE_factor);
title('Power TE factor')
ylabel('Height of the core')
xlabel('Width of the core')
view([0 0 1]);
axis equal

subplot(3,2,4)
surf(mesh_p_b,mesh_p_a,mesh_PowerHIC_TM_factor);
title('Power TM factor')
ylabel('Height of the core')
xlabel('Width of the core')
view([0 0 1]);
axis equal

subplot(3,2,5)
surf(mesh_b_i,mesh_a_i,mesh_PowerHIC_TE_factor_i);
title('Power TE factor')
ylabel('Height of the core')
xlabel('Width of the core')
view([0 0 1]);
axis equal

subplot(3,2,6)
surf(mesh_b_i,mesh_a_i,mesh_PowerHIC_TM_factor_i);
title('Power TM factor')
ylabel('Height of the core')
xlabel('Width of the core')
view([0 0 1]);
axis equal

figure(6)

subplot(3,2,1)
contourf(mesh_p_b,mesh_p_a,mesh_PowerHIC_TE_factor_r);
title('Normalized overlap at the sidewalls, TE')
ylabel('Height of the core')
xlabel('Width of the core')
view([0 0 1]);
axis equal

subplot(3,2,2)
contourf(mesh_p_b,mesh_p_a,mesh_PowerHIC_TM_factor_r);
title('Normalized overlap at the sidewalls, TM')
ylabel('Height of the core')
xlabel('Width of the core')
view([0 0 1]);
axis equal

subplot(3,2,3)
contourf(mesh_p_b,mesh_p_a,mesh_PowerHIC_TE_factor);
title('Normalized overlap at the sidewalls, TE')
ylabel('Height of the core')
xlabel('Width of the core')
view([0 0 1]);
axis equal

subplot(3,2,4)
contourf(mesh_p_b,mesh_p_a,mesh_PowerHIC_TM_factor);
title('Normalized overlap at the sidewalls, TM')
ylabel('Height of the core')
xlabel('Width of the core')
view([0 0 1]);
axis equal

subplot(3,2,5)
[C,h]=contourf(mesh_b_i,mesh_a_i,mesh_PowerHIC_TE_factor_i,-[0:0.05:0.7]);
%h = clabel(C,h);
title('Interpolated TE Factor')
ylabel('Height of the core')
xlabel('Width of the core')
view([0 0 1]);
axis equal

subplot(3,2,6)
[C,h]=contourf(mesh_b_i,mesh_a_i,mesh_PowerHIC_TM_factor_i,-[0:0.05:0.7]);
%h = clabel(C,h);
title('Interpolated TM Factor')
ylabel('Height of the core')
xlabel('Width of the core')
view([0 0 1]);
axis equal

figure(7)

subplot(2,2,1)
[C,h]=contourf(mesh_b_i,mesh_a_i,-mesh_PowerHIC_TE_total_i,-[0,0.005,0.01:0.01:0.23]);
h = clabel(C,h);
title('Interpolated Total Power TE')
ylabel('Height of the core')
xlabel('Width of the core')
view([0 0 1]);
axis([0.2,1,0.2,1]*um)
colormap gray

subplot(2,2,2)
[C,h]=contourf(mesh_b_i,mesh_a_i,-mesh_PowerHIC_TM_total_i,-[0,0.005,0.01:0.01:0.23]);
h = clabel(C,h);
title('Interpolated Total Power TM')
ylabel('Height of the core')
xlabel('Width of the core')
view([0 0 1]);
axis([0.2,1,0.2,1]*um)
colormap gray

subplot(2,2,3)
[C,h]=contourf(mesh_b_i,mesh_a_i,-mesh_PowerHIC_TE_factor_i,-[0.1:0.05:0.7]);
h = clabel(C,h);
title('Interpolated TE Factor')
ylabel('Height of the core')
xlabel('Width of the core')
view([0 0 1]);
axis([0.2,1,0.2,1]*um)
colormap gray

subplot(2,2,4)
[C,h]=contourf(mesh_b_i,mesh_a_i,-mesh_PowerHIC_TM_factor_i,-[0.1:0.05:0.7]);
h = clabel(C,h);
title('Interpolated TM Factor')
ylabel('Height of the core')
xlabel('Width of the core')
view([0 0 1]);
axis([0.2,1,0.2,1]*um)
colormap gray

figure(8)
[C,h]=contourf(mesh_b_i,mesh_a_i,-mesh_PowerHIC_TE_total_i,-[0,0.005,0.01:0.01:0.23]);%-[0,0.1,0.2,0.3,0.4,0.5:0.5:6.5]);
title('Interpolated Total Power, TE')
ylabel('Height of the core')
xlabel('Width of the core')
view([0 0 1]);
axis([0.2,1,0.2,1]*um)
colormap gray

figure(9)
[C,h]=contourf(mesh_b_i,mesh_a_i,-mesh_PowerHIC_TM_total_i,-[0,0.005,0.01:0.01:0.23]);%-[0,0.1,0.2,0.3,0.4,0.5:0.5:6.5]);
title('Interpolated Total Power, TM')
ylabel('Height of the core')
xlabel('Width of the core')
view([0 0 1]);
axis([0.2,1,0.2,1]*um)
colormap gray

figure(10)
[C,h]=contourf(mesh_b_i,mesh_a_i,-mesh_PowerHIC_TE_factor_i,-[0.1:0.05:0.7]);%-[0:0.03:0.45]);
title('TE Factor')
ylabel('Height of the core')
xlabel('Width of the core')
view([0 0 1]);
axis([0.2,1,0.2,1]*um)
colormap gray

figure(11)
[C,h]=contourf(mesh_b_i,mesh_a_i,-mesh_PowerHIC_TM_factor_i,-[0.1:0.05:0.7]);%-[0:0.03:0.45]);
title('TM Factor')
ylabel('Height of the core')
xlabel('Width of the core')
view([0 0 1]);
axis([0.2,1,0.2,1]*um)
colormap gray