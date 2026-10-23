% STAGE 1 - boundary-field extraction sweep for SiON/SiO2 (n_core = 1.50)
%
% Sweeps the (height, width) grid of waveguide cross-sections used in
% Barwicz & Haus, J. Lightwave Technol. 23(9), 2719 (2005), solves each
% cross-section for its two lowest-order modes, and reduces each mode to
% the handful of sidewall quantities the scattering model needs.  Writes
% a struct array `data(ia,ib)` -- see docs/data_format.md.
%
% ORIGIN: derived by Tymon Barwicz from a mode-solver driver script he was
% given at MIT.  The physics (geometry sweep, TE/TM sorting, sidewall
% localization, boundary-field reduction) is his.
%
% *** THIS SCRIPT DOES NOT RUN AS SHIPPED. ***
% The 2D vector mode solver it was written against is not included with this
% project, and its call sites have been replaced by the stubs
% `modesolve_cross_section` and `mode_power`.  Supply your own solver in
% those two files to run stage 1 -- see docs/modesolver_interface.md.
% Precomputed stage-1 output for all three material systems is shipped in
% data/boundary_fields/, so stages 2 and 3 run without a mode solver.

% P - input structure; Q - output structure

close all;clear all;
global VMODE; VMODE = 1;
c = 299792458; um = 1e-6; nm = 1e-9;
t00 = clock;                            % Code timing

NMODES_CALC = 2;                        % How many modes to calculate?
NMODES_KEEP = 2;                        % How many of those modes to keep in the end?


% a is the total height of the waveguide in um
% b is the total width of the waveguide in um
max_a=4;
min_a=0.5;
max_b=4;
min_b=0.5;

% increments in width and height of the waveguide
% must be synchronized with dxy
da=0.05;
db=0.05;
domain=[18,17];       % calculation domaine(x,y) size
dxy = [0.05 0.05];    % Coarse Discretization


% load boundary_fields_SiON_4  %%% start from last row.

data(1,1).dxy=dxy;

for ia=length(data(:,1)):round((max_a-min_a)/da)
    a=max_a-(ia-1)*da;
    ib=0;
    exit_ib=0;
    
    while exit_ib==0
        ib=ib+1;
        b=max([min_b,round((2.5-a)/dxy(1))*dxy(1)])+(ib-1)*db; % the boundaries are defined here.
        if b>=min([max_b,5-a])
            exit_ib=1;
        end
        % Structure parameters for waveguide set P(k = 1..N)
		% P(k) contains data for the kth waveguide.
		P.l0 = 1550 * nm;                       % Freespace wavelength
		P.nco = 1.50; 
		P.ncl = 1.45;
		P.wd = um * b;
		P.t  = um * a;
	    
		P   % PRINT OUT THE P STRUCTURE
			
		% Generate waveguide index profile and computational domain
            x = [-domain(1)/2 : dxy(1)/2 : domain(1)/2]*um;           % Index points x coordinate with fine discretization 
            y = [-domain(2)/2 : dxy(2)/2 : domain(2)/2]*um;           % Index points y coordinate with fine discretization
                                                    % vectors go from -2.5 to 2.5 (for 5x5 domain) with dxy/2 steps in case 
                                                     
            % creating the actual refractive index points
            xx = x(:) * ones(1,length(y));            %the vector x in vertical with length(y) identical columns 
            yy = ones(length(x),1) * y(:).';          %the vector y in horizontal with length(x) identical rows
            
            n = ones(length(x), length(y)) * P.ncl;   % array with all points with the value of ncl  
            
            % Proper odd pixel size treatment
            if (mod(round(P.wd/dxy(1)/um),2)==0)&&(mod(round(P.t/dxy(2)/um),2)==0)
                n( find((abs(yy) <= P.t/2) .* (abs(xx) <=P.wd/2)) ) = P.nco;    % change all values of n to nco were the corresponding coordinate is within the core  
            elseif (mod(round(P.wd/dxy(1)/um),2)==0)&&(mod(round(P.t/dxy(2)/um),2)~=0)
                n( find((yy <=(P.t-dxy(2)*um)/2+dxy(2)*um).* (yy >=-(P.t-dxy(2)*um)/2) .* (abs(xx) <=P.wd/2)) ) = P.nco;
            elseif (mod(round(P.wd/dxy(1)/um),2)~=0)&&(mod(round(P.t/dxy(2)/um),2)==0)
                n( find( (abs(yy) <= P.t/2) .* (xx <=(P.wd-dxy(1)*um)/2+dxy(1)*um) .* (xx >=-(P.wd-dxy(1)*um)/2) ) ) = P.nco;
            elseif (mod(round(P.wd/dxy(1)/um),2)~=0)&&(mod(round(P.t/dxy(2)/um),2)~=0)
                n( find( ((yy <=(P.t-dxy(2)*um)/2+dxy(2)*um).* (yy >=-(P.t-dxy(2)*um)/2) .* (xx <=(P.wd-dxy(1)*um)/2+dxy(1)*um) .* (xx >=-(P.wd-dxy(1)*um)/2) ) ) )= P.nco;
            else
                error('ERROR in proper odd pixel size treatment')
            end   
                
            % Dielectric constant averaging for more accurate answers [MP] Oct 15, 2002
            P.MX = (length(x)-1)/2 +1; P.MY = (length(y)-1)/2 +1;                        % center of the array
            n(3:2:2*P.MX-3,:) = sqrt( (n(2:2:2*P.MX-4,:).^2 + n(4:2:2*P.MX-2,:).^2)/2 ); % for every other index-array points at some distance to the boundary, average it out in one dimension with the second other point to the left and to the right
            n(:,3:2:2*P.MY-3) = sqrt( (n(:,2:2:2*P.MY-4).^2 + n(:,4:2:2*P.MY-2).^2)/2 ); % for every other index-array points at some distance to the boundary, average it out in one dimension with the second other point to up and down
            % Create index data structure
            N.n = n; N.x = x(:); N.y = y(:);  
            
		% UNCOMMENT LINE BELOW TO SEE DIELECTRIC STRUCTURE
		% figure; imagesc(N.x,N.y,N.n.'); set(gca,'YDir','normal'); axis image; break; %the shown index profile is the transpose of the actual array
		
            % Set up and call modesolver
            OPTS.fieldmode = 'V';       % 'V'=vectorial, 'MX'=semivectorial Ex, 'MY'=semivec Ey
            k0 = 2*pi/P.l0;
		    %OPTS.sigma='LR'
            OPTS.coordmode = 'C'; OPTS.eigmode = 'b';         % Cartesian mode, beta eigenvalue
            
            %%% MODIFY THIS FOR FASTER CONVERGENCE
                        %%% TO ENSURE CONBERGENCE WE START FROM CORE INDEX ALWAYS 
            %if ib==1
            %   if ia==1
                   mu_guess = P.nco * k0;                            % Guess propagation constant (can guess core index)
                   %   else
                %   mu_guess = min([P.nco * k0,data(ia-1,ib).beta(1)*1.01]);
                   %   end
                   %else
               %mu_guess = min([P.nco * k0,data(ia,ib-1).beta(1)*1.01]);
               %end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            PML_WIDTH = [0 0 0 0]*um; PML_STRENGTH = [0.1/um 0.1/um]*4;     % Zero width, i.e. NO PMLs.
            
            % CALL TO MODESOLVER:
            [beta, F] = modesolve_cross_section(N, k0, mu_guess, OPTS, NMODES_CALC, PML_WIDTH, PML_STRENGTH);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
                
            alpha = imag(beta);
            dBum = 20E-6*alpha * log10(exp(1));
		
		neff=beta/k0
            
		% Reducing structure size and dumping to data file
		F.Ex = F.Ex(:,:,1:NMODES_KEEP); F.Ey = F.Ey(:,:,1:NMODES_KEEP); F.Ez = F.Ez(:,:,1:NMODES_KEEP); %just keep the modes you want, irrelevant if interested in all calculated modes
		F.Hx = F.Hx(:,:,1:NMODES_KEEP); F.Hy = F.Hy(:,:,1:NMODES_KEEP); F.Hz = F.Hz(:,:,1:NMODES_KEEP); %just keep the modes you want
		F.dx = dxy(1); F.dy = dxy(2);
		
		% SZ are the complete results, putting stuff togheter for saving purposes
		SZ.P    = P;
		SZ.dBum = dBum;
		SZ.beta = beta;
		SZ.F    = F;
		SZ.N    = N;            % Save index structure
		

		% For diagnostic purposes, plot dielectric and field profiles:
		%for m = 1:length(SZ)
        %    for k = 1:size(SZ(m).F.Ex,3)
        %        figure;
        %        subplot(2,2,1); imagesc(SZ(m).N.x /um,  SZ(m).N.y /um,  SZ(m).N.n.'); set(gca, 'ydir', 'normal'); ylabel('\mum');
        %        title(['Structure ' num2str(m) ', mode ' num2str(k)]);
        %        subplot(2,2,2); imagesc(SZ(m).F.Rr /um, SZ(m).F.Zr /um, real(SZ(m).F.Ex(:,:,k).')); set(gca, 'ydir', 'normal'); title('Re(Ex)');
        %        subplot(2,2,3); imagesc(SZ(m).F.Rz /um, SZ(m).F.Zz /um, real(SZ(m).F.Ey(:,:,k).')); set(gca, 'ydir', 'normal'); title('Re(Ey)'); xlabel('\mum'); ylabel('\mum');
        %        subplot(2,2,4); imagesc(SZ(m).F.Rz /um, SZ(m).F.Zr /um, imag(SZ(m).F.Ez(:,:,k).')); set(gca, 'ydir', 'normal'); title('Im(Ez)'); xlabel('\mum');
        %    end
        %end
		
		% Need to decide here on what is TE and what is TM
		% USE THE VARIABLES TE AND TM TO TRANSFER TO THE SCAT LOSS CALC
		% Sum the power in Ex and Ey and compare for each mode
		
		[P1, SzTM1, SzTE1] = mode_power(F, 1);
		[P2, SzTM2, SzTE2] = mode_power(F, 2);
		
		PEx1=sum(sum(abs(SzTE1)));
		PEx2=sum(sum(abs(SzTE2)));
		PEy1=sum(sum(abs(SzTM1)));
		PEy2=sum(sum(abs(SzTM2)));
		
		% Put the right fields in the following variables
		% OCTAVE PORT NOTE (2026-08-16): 'clear TE TM' added; see boundary_fields_Si_v3.m
		clear TE TM
		if (PEx1>PEy1)&(PEx2<PEy2)
            TE.Ex=F.Ex(:,round(length(F.Ex(1,:,1))/2),1); %use only a near-center line 
            TE.beta=beta(1);
            TE.F.Ex=F.Ex(:,:,1);
            TE.F.Ez=F.Ez(:,:,1);
            TE.F.Ey=F.Ey(:,:,1);
            TE.P=P1;
            TM.Ey=F.Ey(round(length(F.Ey(:,1,2))/2),:,2);
            TM.beta=beta(2);
            TM.F.Ey=F.Ey(:,:,2);
            TM.F.Ez=F.Ez(:,:,2);
            TM.F.Ex=F.Ex(:,:,2);
            TM.P=P2;
            disp('The foundamental mode is TE-like')
		elseif (PEx1<PEy1)&(PEx2>PEy2)
            TE.Ex=F.Ex(:,round(length(F.Ex(1,:,2))/2),2);
            TE.beta=beta(2);
            TE.F.Ex=F.Ex(:,:,2);
            TE.F.Ez=F.Ez(:,:,2);
            TE.F.Ey=F.Ey(:,:,2);
            TE.P=P2;
            TM.Ey=F.Ey(round(length(F.Ey(:,1,1))/2),:,1);
            TM.beta=beta(1);
            TM.F.Ey=F.Ey(:,:,1);
            TM.F.Ez=F.Ez(:,:,1);
            TM.F.Ex=F.Ex(:,:,1);
            TM.P=P1;
            disp('The foundamental mode is TM-like')
		else
            % If both modes are TE or TM ... send a warning and take the first mode
            % as TE and the second as TM
            warning('Both modes are either TE or TM. Degeneracy may create wrong Ez pattern.')
            P.wd
            P.t
            TE.Ex=F.Ex(:,round(length(F.Ex(1,:,1))/2),1); %use only a near-center line 
            TE.beta=beta(1);
            TE.F.Ex=F.Ex(:,:,1);
            TE.F.Ez=F.Ez(:,:,1);
            TE.F.Ey=F.Ey(:,:,1);
            TE.P=P1;
            TM.Ey=F.Ey(round(length(F.Ey(:,1,2))/2),:,2);
            TM.beta=beta(2);
            TM.F.Ey=F.Ey(:,:,2);
            TM.F.Ez=F.Ez(:,:,2);
            TM.F.Ex=F.Ex(:,:,2);
            TM.P=P2;
            disp('Taking first mode as TE-like')
            
		end    
		
		% Need to look for the field minimum to find the needed boundary
		%figure;
		%subplot(1,2,1)
		%plot(abs(TE.Ex))
		%title('Ex horizontal profile of TE mode')
		%xlabel('x grid points')
		%subplot(1,2,2)
		%plot(abs(TM.Ey))
		%title('Ey vertical profile of TM mode')
		%xlabel('y grid points')
		
		%Starting points for local minima search
		x_start=round((domain(1)/dxy(1)+1)/2+[1,-1]*(P.wd/(dxy(1)*um)-1)/4)
		y_start=round((domain(2)/dxy(2)+1)/2+[1,-1]*(P.t/(dxy(2)*um)-1)/4)
		
		%Looking for local minima
        data(ia,ib).error=0; 
		for ie=1:(length(x_start)+length(y_start))
            % set the starting index
            if ie <= length(x_start)
                E=TE.Ex;
                ix=x_start(ie);
            else
                E=TM.Ey;
                ix=y_start(ie-length(x_start));
            end
            %find the trend of the field profile
            if abs(E(ix))>=abs(E(ix-1))
                direction=-1;
            else
                direction=+1;
            end
            exit=0;
            %look for the first relative minimum
            while exit==0
                 if abs(E(ix))>=abs(E(ix+direction))
                    ix=ix+direction;
                 else
                    ibound=ix;
                    exit=1;
                 end
                 % convergence test
                 if ix==length(E) %edges of data---------
                    ie
                    warning('Upper edge of data encountered in field minimum search')
                    data(ia,ib).error=1;
                    exit=1;
                elseif ix==1
                    ie
                    warning('Lower edge of data encountered in field minimum search')
                    data(ia,ib).error=1;
                    exit=1;  
                end              
             end     
            clear exit
            % re-arrange the data
            if (ie <= length(x_start))&(data(ia,ib).error==0)
                x_bound(ie)=ibound;
            elseif data(ia,ib).error==0
                y_bound(ie-length(x_start))=ibound;
            end  
            if data(ia,ib).error==1
             	x_bound=round((domain(1)/dxy(1)+1)/2+[1,-1]*(P.wd/(dxy(1)*um)-1)/2);
        		y_bound=round((domain(2)/dxy(2)+1)/2+[1,-1]*(P.t/(dxy(2)*um)-1)/2);
            end  
		end
		
		%CAREFUL, THE MODE SOLVER'S COORDINATES ARE NOT LACEY'S COORDINATES X->Y Y->X 
		% Also, those are peak fields and peak powers.
		TE.R_YZ=(trapz(abs(TE.F.Ez(x_bound(1),y_bound(2):y_bound(1))).^2)+trapz(abs(TE.F.Ez(x_bound(2),y_bound(2):y_bound(1))).^2))/(trapz(abs(TE.F.Ex(x_bound(1),y_bound(2):y_bound(1))).^2)+trapz(abs(TE.F.Ex(x_bound(2),y_bound(2):y_bound(1))).^2));
		TM.R_XZ=(trapz(abs(TM.F.Ez(x_bound(1),y_bound(2):y_bound(1))).^2)+trapz(abs(TM.F.Ez(x_bound(2),y_bound(2):y_bound(1))).^2))/(trapz(abs(TM.F.Ey(x_bound(1),y_bound(2):y_bound(1))).^2)+trapz(abs(TM.F.Ey(x_bound(2),y_bound(2):y_bound(1))).^2));
		
		% Field overlap normalized to guided power. % those are peak fields and
		% peak powers. Divided by 2 to have overlap at one interface only: average
		% from the two interfaces.
        TE.Phi_squared=((trapz(abs(TE.F.Ez(x_bound(1),y_bound(2):y_bound(1))).^2)+trapz(abs(TE.F.Ez(x_bound(2),y_bound(2):y_bound(1))).^2))+(trapz(abs(TE.F.Ex(x_bound(1),y_bound(2):y_bound(1))).^2)+trapz(abs(TE.F.Ex(x_bound(2),y_bound(2):y_bound(1))).^2))+(trapz(abs(TE.F.Ey(x_bound(1),y_bound(2):y_bound(1))).^2)+trapz(abs(TE.F.Ey(x_bound(2),y_bound(2):y_bound(1))).^2)))/2/TE.P;
		TM.Phi_squared=((trapz(abs(TM.F.Ez(x_bound(1),y_bound(2):y_bound(1))).^2)+trapz(abs(TM.F.Ez(x_bound(2),y_bound(2):y_bound(1))).^2))+(trapz(abs(TM.F.Ex(x_bound(1),y_bound(2):y_bound(1))).^2)+trapz(abs(TM.F.Ex(x_bound(2),y_bound(2):y_bound(1))).^2))+(trapz(abs(TM.F.Ey(x_bound(1),y_bound(2):y_bound(1))).^2)+trapz(abs(TM.F.Ey(x_bound(2),y_bound(2):y_bound(1))).^2)))/2/TM.P;

        
        %clearing stuff
        TE.F=0;
        TM.F=0;
        TE.Ex=0;
        TM.Ey=0;
        
        data(ia,ib).P=P; % saving everything relevant in structure data               
        data(ia,ib).a=a*um;
        data(ia,ib).b=b*um;
        data(ia,ib).beta=beta;
        data(ia,ib).TE=TE;
        data(ia,ib).TM=TM;
        %data(ia,ib).N=N;
        data(ia,ib).x_bound=x_bound;
        data(ia,ib).y_bound=y_bound;   
        
   		% Save computed modes and overlaps
		save(mfilename,'data');   % Save all relevant data

    end
end