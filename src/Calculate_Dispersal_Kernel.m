%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Subfunction Seed Dispersal             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spec] = Calculate_Dispersal_Kernel(isp,spec,cellsize,doFFT,maxy,maxx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate and store dispersal kernel of species
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Inputs
%%% isp [-] the species index
%%% spec [-] species properties: kernel value, radius
%%% cellsize [m] the cellsize
%%% doFFT [-] whether dispersal is done by FFT or not
%%% maxy, maxx [] the rows and columns of the simulation region
%%% Outputs
%%% spec [] add kernel, radi (+ kernel_ft, if FFT wished)

firstKernel = true;
%%%% If same parameters already existed in other species, take the kernel and store it.
% for ii=1:isp-1
%     if (spec{isp}.alpha1==spec{ii}.alpha1)&&(spec{isp}.alpha2==spec{ii}.alpha2)&&(spec{isp}.dispfrac==spec{ii}.dispfrac)
%         spec{isp}.radi = spec{ii}.radi;
%         spec{isp}.kernel = spec{ii}.kernel;
%         firstKernel = false;
%         break
%     end
% end
%%%% If the parameters appear first time, kernel has to be calculated.
if firstKernel
    kerneltype = 1;
    if spec{isp}.dispfrac<1.0
        kerneltype = 2;
    end
    [spec] = Dispersal_Kernel(isp,spec,kerneltype,cellsize); %%% spec: + kernel, radi
end

% if doFFT  %%% current version no FFT
%     spec = Set_Dims_FFT(isp,spec,maxy,maxx);
%     spec = Calculate_Kernel_Transform(isp,spec); %%% spec: + kernel_ft
% end

end





%%%% spec  cell in structure
function [spec] = Dispersal_Kernel(isp,spec,kerneltype,cellsize)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate the general dispersal kernel of one species
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Inputs
%%% isp [-] the species index
%%% spec [-] species properties: dispfrac, alpha1, alpha2
%%% cellsize [m] the cellsize
%%% Outputs
%%% radi [] the radius of kernel (in grid cell number)
%%% spec [-] species properties: dispfrac, alpha1, alpha2 + kernel, radi

epskernel = 10^(-4);   %%% threshold for cut-off of kernel, all values smaller than this are set to 0

if kerneltype==1
    spro = spec{isp};
    [radi,kernel] = Single_Kernel(spro.alpha1,cellsize,epskernel);
elseif kerneltype==2
    spro = spec{isp};
    [radi,kernel] = Double_Kernel(spro.dispfrac,spro.alpha1,spro.alpha2,cellsize,epskernel);
end
clear spro
%%% store the kernel and radius in the "spec"
spec{isp}.kernel = kernel;
spec{isp}.radi = radi;
end


function [radi,kernel] = Double_Kernel(dispfrac,alpha1,alpha2,cellsize,epskernel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate the double dispersal of one species by adding two single kernels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Inputs
%%% dispfrac [-] the proportion of short distance dispersal
%%% alpha1 [m] the species' dispersal distance for short distance
%%% alpha2 [m] the species' dispersal distance for long distance
%%% cellsize [m] the cellsize
%%% epskernel [-] threshold for cut-off of kernel, all values smaller than this are set to 0
%%% Outputs
%%% radi [] the radius of kernel (in grid cell number)
%%% kernel [] the kernel value

if (1-dispfrac)>0
    %%% the kernel for long distance
    [rad2,kernel2] = Single_Kernel(alpha2,cellsize,epskernel);
    %%% the kernel for short distance
    [rad1,kernel1] = Single_Kernel(alpha1,cellsize,epskernel);
    %%% combine them
    radi = rad2;
    %central = [radi+1,radi+1];  %%% the row and colume of the central point in Matlabe, i.e., (0,0) in the Cartesian
    kernel = (1-dispfrac).*kernel2;     %%% long kernel
    kernel(radi+1-rad1:radi+1+rad1,radi+1-rad1:radi+1+rad1) = kernel(radi+1-rad1:radi+1+rad1,radi+1-rad1:radi+1+rad1) + dispfrac.*kernel1;  %%% long kernel+short kernel
else
    [radi,kernel] = Single_Kernel(alpha1,cellsize,epskernel);
end

end


function [radi,kernel] = Single_Kernel(alpha,cellsize,epskernel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate the single dispersal kernel of one species
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Inputs
%%% alpha [m] the species' dispersal distance
%%% cellsize [m] the cellsize
%%% epskernel [-] threshold for cut-off of kernel, all values smaller than this are set to 0
%%% Outputs
%%% radi [] the radius of kernel (in grid cell number)
%%% kernel [] the kernel value

%%%% Set the maxmium radius that the seeds can reach => set the minimum kernel value (epskernel) 
%%%% kernel (x) = exp(-x/alpha) => epskernel = exp(-x_eps/alpha) => x_eps = - alpha * ln(epskernel)

radfac = min(20, -log(epskernel)); %% [m] at which distance has the kernel function decreased down to the epskernel value?
radm = alpha*radfac;               %% [m] radius of actual discretization (in meters)
radi = round(radm/cellsize);     %% [-] radius of actual discretization (in grid cell number)
%%%%% Determine settings of fine discretization kernel
cellsizemin = cellsize;%100; %1000/40;      %% cell size of totally finest discretization, about 0.0001 of kernel cut off
cellsizemin = min(cellsizemin, cellsize); %% not larger than actual cell size
icellfrac = cellsize/cellsizemin;  %% ratio actual to finest discretization
radmax = round(radm/cellsizemin);   %% radius of finest discretization in number of grid cells (The 20 is sufficient at small gridsizes)

%%%%%%%% Fine discretization kernel %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% fine discretization kernel, calculate one quadrant
kernelfine1 = zeros(radmax+1,radmax+1);
for xwin=0:radmax
    for ywin=0:radmax
        kernelfine1(radmax-ywin+1,xwin+1) = Center_to_Center_Value(alpha,cellsizemin,xwin,ywin);
    end
end
%%%% fine discretization kernel, copy this one quadrant to other quadrants
kernelfine2 = fliplr(kernelfine1); %% the second quadrant
kernelfine1 = [kernelfine2(:,1:end-1) kernelfine1];
kernelfine = flipud(kernelfine1);
kernelfine = [kernelfine1; kernelfine(2:end,:)]; %%% four quadrants
clear kernelfine2

%%%%%%%% Actual kernel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% The kernel on the coarse (i.e., actual) grid is calculated by summing over the kernel values in all fine-grid cell in each coarse grid cell 
kernel = zeros(2*radi+1,2*radi+1);
for xwinr=-radi:radi
    for ywinr=-radi:radi
        xwinsta = max(min(radmax, xwinr*icellfrac - round(icellfrac/2)), -radmax);
        xwinend = max(min(radmax, xwinr*icellfrac + round(icellfrac/2)), -radmax);
        ywinsta = max(min(radmax, ywinr*icellfrac - round(icellfrac/2)), -radmax);
        ywinend = max(min(radmax, ywinr*icellfrac + round(icellfrac/2)), -radmax);
        kernel(radi+1+ywinr,radi+1+xwinr) = kernel(radi+1+ywinr,radi+1+xwinr) + sum(kernelfine(radmax+1+ywinsta:radmax+1+ywinend,radmax+1+xwinsta:radmax+1+xwinend),'all')/2;
        % for xwin=xwinstart:xwinend
        %     for ywin=ywinstart:ywinend
        %         kernel(radi+1+ywinr,radi+1+xwinr) = kernel(radi+1+ywinr,radi+1+xwinr) + kernelfine(radmax+1+ywin,radmax+1+xwin)/2;
        %     end
        % end
    end
end

%%%% Normalize kernel, sum = 1
krnlsm = sum(kernel,'all');
kernel = kernel./krnlsm;

end




function [ccvalue] = Center_to_Center_Value(alpha,cellsize,x1,y1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate kernel value center to center
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Inputs
%%% alpha [m] the species' dispersal distance
%%% cellsize [m] the cell size
%%% x1 [row] how many rows between the source and sink cells
%%% y1 [col] how many columns between the source and sink cells
%%% Outputs
%%% ccvalue [-] the kernel value center to center

dist = sqrt(x1.^2+y1.^2).*cellsize;
ccvalue = Kernel_Value(dist,alpha);
end



function [kv] = Kernel_Value(dist,alpha)
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Inputs
%%% dist [m] the distance from source to sink 
%%% alpha [m] the species' dispersal distance
%%% Outputs
%%% kv [-] the kernel value

kv = exp(-dist./alpha); %% It should be set to zero if smaller than epskernel
end















