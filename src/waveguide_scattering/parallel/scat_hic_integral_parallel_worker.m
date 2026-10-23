function []=scat_hic_integral_parallel_worker(filename, ia_list, shard_id)
% Parallel-shard worker for the stage-2 far-field integration.
%
% PROVENANCE: the per-cell body below is transcribed VERBATIM from Tymon
% Barwicz's scat_hic_integral_final_v3.m (which remains unmodified in the
% parent directory). The only changes are:
%   - the outer row loop runs over the caller-supplied ia_list instead of
%     1:length(data(:,1))
%   - results save to '<filename>_integrated_shard<id>' after every cell
% All physics calls (scat_3d_payne, scat_hic_integrant_*, adaptlob, dblquad),
% parameters, tolerances, and constants are identical to the serial version.
% Deterministic quadrature => per-cell results are identical to a serial run.
%
% Merge shards with scat_hic_integral_parallel_merge.m.
% GNW / Claude, 2026-08-19.

load(filename)
savefilename=strcat(filename,'_integrated_shard',num2str(shard_id));

% guard: drop requested rows beyond the sweep (row counts vary per material)
out_of_range = ia_list(ia_list > size(data,1));
if ~isempty(out_of_range)
  warning('ignoring out-of-range rows: %s (data has %d rows)', ...
          mat2str(out_of_range), size(data,1));
  ia_list = ia_list(ia_list <= size(data,1));
end

dxy=data(1,1).dxy;

um=1e-6;
nm=1e-9;

start=1e-10;
tol=1e-4;
Lc=50*nm;
sigma=1*nm;

for ia=ia_list(:).'
    ib=1;
    exit_ib=0;
    while exit_ib==0

        tic
        % required variables
        a=data(ia,ib).P.t; %%%
        b=data(ia,ib).P.wd; %%%
        ncore=data(ia,ib).P.nco;
        nclad=data(ia,ib).P.ncl;
        lambda=data(ia,ib).P.l0;
        ia
        ib

        % TE Polarization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % q = quad(fun,a,b,tol,trace,p1,p2,...)
        % the 2D reference
        PaynePower_TE(ia,ib)=adaptlob(@scat_3d_payne,start,2*pi,tol,[],data(ia,ib).TE.beta,lambda,nclad,ncore,Lc,sigma);

        % the X part
        R_YZ=0;
        PowerHIC_TE_Y_cos_1(ia,ib)=dblquad(@scat_hic_integrant_YZ_cos_1,start,2*pi,start,pi/2,tol,@adaptlob,a,b,data(ia,ib).TE.beta,R_YZ,lambda,nclad,ncore,Lc,sigma);
        PowerHIC_TE_Y_cos_2(ia,ib)=dblquad(@scat_hic_integrant_YZ_cos_2,start,2*pi,pi/2,pi,tol,@adaptlob,a,b,data(ia,ib).TE.beta,R_YZ,lambda,nclad,ncore,Lc,sigma);

        % the Z part
        R_YZ=1e6;
        PowerHIC_TE_Z_cos_1(ia,ib)=dblquad(@scat_hic_integrant_YZ_cos_1,start,2*pi,start,pi/2,tol,@adaptlob,a,b,data(ia,ib).TE.beta,R_YZ,lambda,nclad,ncore,Lc,sigma);
        PowerHIC_TE_Z_cos_2(ia,ib)=dblquad(@scat_hic_integrant_YZ_cos_2,start,2*pi,pi/2,pi,tol,@adaptlob,a,b,data(ia,ib).TE.beta,R_YZ,lambda,nclad,ncore,Lc,sigma);

        % combined integrals to double check the calculations
        R_YZ=sqrt(data(ia,ib).TE.R_YZ)*i; %sqrt compensates for changing from R in field to R in power
        PowerHIC_TE_YZ_cos_1(ia,ib)=dblquad(@scat_hic_integrant_YZ_cos_1,start,2*pi,start,pi/2,tol,@adaptlob,a,b,data(ia,ib).TE.beta,R_YZ,lambda,nclad,ncore,Lc,sigma);
        PowerHIC_TE_YZ_cos_2(ia,ib)=dblquad(@scat_hic_integrant_YZ_cos_2,start,2*pi,pi/2,pi,tol,@adaptlob,a,b,data(ia,ib).TE.beta,R_YZ,lambda,nclad,ncore,Lc,sigma);

        % combining everything
        PowerHIC_TE_norm(ia,ib)=(1/(1+data(ia,ib).TE.R_YZ))*((PowerHIC_TE_Y_cos_1(ia,ib)+PowerHIC_TE_Y_cos_2(ia,ib))+data(ia,ib).TE.R_YZ*(PowerHIC_TE_Z_cos_1(ia,ib)+PowerHIC_TE_Z_cos_2(ia,ib)));
        PowerHIC_TE_comb(ia,ib)=(PowerHIC_TE_YZ_cos_1(ia,ib)+PowerHIC_TE_YZ_cos_2(ia,ib));
        PowerHIC_TE_test(ia,ib)=PowerHIC_TE_comb(ia,ib)/PowerHIC_TE_norm(ia,ib); %%%
        PowerHIC_TE_factor(ia,ib)=PowerHIC_TE_norm(ia,ib)/PaynePower_TE(ia,ib);  %%%
        %in (dB/cm)/nm^2 % the field overlap is /um instead of /m % on top of that, phi squared should be the pic field normalized by the time averaged power... so factor of 2
        PowerHIC_TE_total(ia,ib)=PowerHIC_TE_norm(ia,ib)*(data(ia,ib).TE.Phi_squared*2*dxy(2)/um)*log10(exp(1))/10;


        % TM Polarization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % the 2D reference
        PaynePower_TM(ia,ib)=adaptlob(@scat_3d_payne,start,2*pi,tol,[],data(ia,ib).TM.beta,lambda,nclad,ncore,Lc,sigma);

        % the X part
        R_XZ=0;
        PowerHIC_TM_X_cos_1(ia,ib)=dblquad(@scat_hic_integrant_XZ_cos_1,start,2*pi,start,pi/2,tol,@adaptlob,a,b,data(ia,ib).TM.beta,R_XZ,lambda,nclad,ncore,Lc,sigma);
        PowerHIC_TM_X_cos_2(ia,ib)=dblquad(@scat_hic_integrant_XZ_cos_2,start,2*pi,pi/2,pi,tol,@adaptlob,a,b,data(ia,ib).TM.beta,R_XZ,lambda,nclad,ncore,Lc,sigma);

        % the Z part
        R_XZ=1e6;
        PowerHIC_TM_Z_sin_1(ia,ib)=dblquad(@scat_hic_integrant_XZ_sin_1,start,2*pi,start,pi/2,tol,@adaptlob,a,b,data(ia,ib).TM.beta,R_XZ,lambda,nclad,ncore,Lc,sigma);
        PowerHIC_TM_Z_sin_2(ia,ib)=dblquad(@scat_hic_integrant_XZ_sin_2,start,2*pi,pi/2,pi,tol,@adaptlob,a,b,data(ia,ib).TM.beta,R_XZ,lambda,nclad,ncore,Lc,sigma);

        PowerHIC_TM_norm(ia,ib)=(1/(1+data(ia,ib).TM.R_XZ))*((PowerHIC_TM_X_cos_1(ia,ib)+PowerHIC_TM_X_cos_2(ia,ib))+data(ia,ib).TM.R_XZ*(PowerHIC_TM_Z_sin_1(ia,ib)+PowerHIC_TM_Z_sin_2(ia,ib)));
        PowerHIC_TM_factor(ia,ib)=PowerHIC_TM_norm(ia,ib)/PaynePower_TM(ia,ib); %%%
        %in (dB/cm)/nm^2 % the field overlap is /um instead of /m % on top of that, phi squared should be the pic field normalized by the time averaged power... so factor of 2
        PowerHIC_TM_total(ia,ib)=PowerHIC_TM_norm(ia,ib)*(data(ia,ib).TM.Phi_squared*2*dxy(2)/um)*log10(exp(1))/10; %%%

        %%%%%%%%%%%%

        matrix_a(ia,ib)=a;
        matrix_b(ia,ib)=b;
        matrix_phi_squared_TM(ia,ib)=data(ia,ib).TM.Phi_squared*2*dxy(2)/um;
        matrix_phi_squared_TE(ia,ib)=data(ia,ib).TE.Phi_squared*2*dxy(2)/um;

        ib=ib+1;
        if (ib<=length(data(1,:)))
            if (data(ia,ib).a>0)
            else
                exit_ib=1;
            end
        else
            exit_ib=1;
        end

        toc
        disp('(seconds)')

        save(savefilename)
    end
end

end
