%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Subfunction spatial dynamic in one time step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [B_H,B_L,ksv,ksvtm1,Ccrown_t] = Spatial_Dynamic_One_Time_Step(curYsta,curYend,curXsta,curXend,doDispersal,doFFT,...
                                outerseeds,stockability,B_H,B_L,ksv,bareId,Parameters_SD,cc_max,cellsize, ...
                                hc_H,hc_L,hc_Htm1,hc_Ltm1,Ccrown_t,Hveg,Lveg,doStochSeedDisp,PHE_S_H,PHE_S_L, ...
                                PHE_S_Htm1,PHE_S_Ltm1,Datam_S,values,fr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Spatial dynamics in one time step
%%% For non parallel, i.e., standard
%%%   1. Loop over the simulation area, for each cell
%%%        a. get the environmental factors
%%%        b. do the local dynamics
%%%   2. update the seeds
%%%   3. do the FFT dispersal, if wished
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Inputs
%%% curYsta,curYend,curXsta,curXend [] part of total simulation area currently simulated, e.g., for spin-up. Is not subarea of parallelisation!
%%% Ccrown [] vegetated fraction for each species type
%%% doDispersal [] flag, simulate dispersal?
%%% doFFT [] flags, dispersal by FFT? ------ Current version don't do FFT.
%%% outerseed [] the seeds from outsider simulation region
%%% stockability [] whether this grid cell can produce seeds or not, e.g., impervious surface
%%% spec [] the properties of this species
%%% seedrain [] store the number of very new seeds - size(maxy,maxx,Ccrown)
%%% B [] carbon pool biomass - size(maxy*maxx,Ccrown,8)
%%% ksv [] land cover in the simulation area - size(maxy*maxx,1)
%%% bareId [] flag - bare land
%%% cc_max [] the maximum species in one cell
%%% cellsize [m] the spatial resolution
%%% hc [m] vegetation height
%%% Ccrown_t [] the vegetation fraction in current time step
%%% Hveg, Lveg [] the location of high and low vegetation
%%% doStochSeedDisp [] whether do stochastic seed dispersal or not [true or false]
%%% PHE_S [] the phenology at this time step. To determine when disperses seeds: the transition from normal growth to senescence (3->4, even for evergreen vegetation)
%%% PHE_Stm1 [] the phenology at last time step.
%%% Datam_S [] the date of this time step
%%% values [] the priority of vegetation, e.g., [1 1 2:26]
%%% fr [] the fraction allocated to fine roots
%%% Outputs
%%% B [] new carbon pool biomass - size(maxy*maxx,Ccrown,8)
%%% ksv [] modified land cover in the simulation area - size(maxy*maxx,1)
%%% ksvtm1 [] the land cover in the simulation area (last time step) - size(maxy*maxx,1)

ksvtm1 = ksv;  %% store the last time step ksv

maxy = curYend-curYsta+1;   %% the number of rows of current simulation domain
maxx = curXend-curXsta+1;   %% the number of columns of current simulation domain
Ccrown_t = reshape(Ccrown_t,maxy,maxx,cc_max);

seedrain = zeros(maxy,maxx,cc_max);

seedprod = reshape(B_H(:,:,5)+B_L(:,:,5),maxy,maxx,cc_max);   %%% seedprod: seed production - size(maxy*maxx,Ccrown)  [[[seedprod = B_H/L(:,Ccrown,5)]]]
ksv = reshape(ksv,maxy,maxx);
hc = reshape(hc_H+hc_L,maxy,maxx,cc_max);
hctm1 = reshape(hc_Htm1+hc_Ltm1,maxy,maxx,cc_max);

PHE_S = reshape(PHE_S_H+PHE_S_L,maxy,maxx,cc_max);
PHE_Stm1 = reshape(PHE_S_Htm1+PHE_S_Ltm1,maxy,maxx,cc_max);

count = 1;
for x=curXsta:curXend
    for y=curYsta:curYend
		if (x>maxx)||(x<1)||(y>maxy)||(y<1)
			disp('Wrong coordinates!')
			disp(x)
			disp(y)
        end
        
        if stockability(y,x)>1e-2
            % disp(ksv(y,x))
            % disp(y)
            % disp(x)
		    [Ccrown,spec,ZR95_H,~,hcM,~,~,~,~] = Parameters_SD(ksv(y,x)); 
		    for isp=1:length(Ccrown)
                doDisDate = (PHE_S(y,x,isp)==4 && PHE_Stm1(y,x,isp)==3) | (isp==length(Ccrown)&&(Datam_S(2)==10)&&(Datam_S(3)==20));
			    if seedprod(y,x,isp)>1e-2 && doDisDate   %% if this cell is stockable, has seeds and is the proper date, do dispersal
				    if doDispersal && ~doFFT && (hc(y,x,isp)>hcM(isp)) && (Ccrown_t(y,x,isp)>0)  %%% when tree becomes mature, it can do seed dispersal.
                                                                                           %%% when this vegetation type exists, seed disperses.
					    [spec] = Calculate_Dispersal_Kernel(isp,spec,cellsize,doFFT,maxy,maxx);
					    [seedrain,ksv,curDisp] = Seed_Dispersal_From_This_Cell_Spec(isp,y,x,stockability,spec,maxy,maxx,seedrain,seedprod,ksv,bareId,cc_max,doStochSeedDisp,values);
					    if doStochSeedDisp && (sum(curDisp,'all')>seedprod(y,x,isp)) %%% modify the seeds if seedrain > seedprod, when doing stochastic seed dispersal. 
                            overseed = (sum(curDisp,'all')-seedprod(y,x,isp))/sum(curDisp,'all').*curDisp;
                            seedrain(:,:,isp) = seedrain(:,:,isp) - overseed;

                            if ZR95_H(isp)>0
                                B_H(count,isp,5) = B_H(count,isp,5) - seedprod(y,x,isp);  %% x is the first for-loop, because in matlab, store in column
                            else
                                B_L(count,isp,5) = B_L(count,isp,5) - seedprod(y,x,isp);
                            end
                        else
                            if ZR95_H(isp)>0
                                B_H(count,isp,5) = B_H(count,isp,5) - sum(curDisp,'all');  %% x is the first for-loop, because in matlab, store in column
                            else
                                B_L(count,isp,5) = B_L(count,isp,5) - sum(curDisp,'all');
                            end
                        end
				    end
			    end
            end
        end
		count = count + 1;
    end
end

%%% change the land cover value after dispersing seeds.
[ksv,Ccrown_t] = Change_land_cover_value(curYsta,curYend,curXsta,curXend,ksv,bareId,hc,hctm1,Parameters_SD,Ccrown_t,cc_max,stockability,Datam_S,values);

ksv = reshape(ksv,maxy*maxx,1);
Ccrown_t = reshape(Ccrown_t,maxy*maxx,cc_max);

%%% Do seed dispersal by FFT is wished.
%if doDispersal && doFFT
%    seedrain = Seedrain_by_FFT(cc_max,maxx,maxy,spec,seedrain,stockability,seedprod);
%end

%%% Update the new seeds ending up in this simulation area
[seedrain,newSeeds] = Update_New_Seeds(cc_max,curYsta,curYend,curXsta,curXend,maxy,maxx,doDispersal,seedrain);

%%% Get seeds from outside the simulation region (newSeeds: finally, arriving seeds of all species in each cell (the entire grid))
newSeeds = Seeds_From_OutSim_Region(cc_max,curYsta,curYend,curXsta,curXend,doDispersal,outerseeds,newSeeds);

%%% seedbank(B(3/4)) adds new seeds.
finerootH = reshape(B_H(:,:,3),maxy,maxx,cc_max);
carbonhyH = reshape(B_H(:,:,4),maxy,maxx,cc_max);
finerootL = reshape(B_L(:,:,3),maxy,maxx,cc_max);
carbonhyL = reshape(B_L(:,:,4),maxy,maxx,cc_max);

fr = repmat(reshape(fr, 1, 1, cc_max), maxy, maxx, 1);
% fr_g = fr(:,:,4);
% fr_g(contains(ksv,'D')) = 0.45;
% fr(:,:,4) = fr_g;

finerootH(:,:,Hveg) = finerootH(:,:,Hveg) + fr(:,:,Hveg).*newSeeds(:,:,Hveg); 
carbonhyH(:,:,Hveg) = carbonhyH(:,:,Hveg) + (1-fr(:,:,Hveg)).*newSeeds(:,:,Hveg);
finerootL(:,:,Lveg) = finerootL(:,:,Lveg) + fr(:,:,Lveg).*newSeeds(:,:,Lveg); 
carbonhyL(:,:,Lveg) = carbonhyL(:,:,Lveg) + (1-fr(:,:,Lveg)).*newSeeds(:,:,Lveg);

B_H(:,:,3) = reshape(finerootH,maxy*maxx,cc_max);
B_H(:,:,4) = reshape(carbonhyH,maxy*maxx,cc_max);
B_L(:,:,3) = reshape(finerootL,maxy*maxx,cc_max);
B_L(:,:,4) = reshape(carbonhyL,maxy*maxx,cc_max);

end




function [seedrain,newSeeds] = Update_New_Seeds(cc_max,ys,ye,xs,xe,maxy,maxx,doDispersal,seedrain)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Update the seeds in the sink cells in a given subarea; in verynewseeds
%%% the current seeds landing in the subarea are collected. Once this is
%%% finished, all the seeds are transferred to the final variable newseeds.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Inputs
%%% cc_max [] the maximum species in one cell
%%% ys, ye, xs, xe [] the start and end of coordinates (row and column) (relative) of current grid subarea
%%% doDispersal [] if explicit dispersal has to be done done; normally the case; corresponds to idisp
%%% seedrain [] store the number of very new seeds - size(maxy,maxx,Ccrown)
%%% Outputs
%%% seedrain [] store the number of very new seeds - size(maxy,maxx,Ccrown)
%%% newSeeds [] finally, arriving seeds of all species in each cell (the entire grid)

if (ys<1)||(ys>maxy) || (ye<1)||(ye>maxy) || (xs<1)||(xs>maxx) || (xe<1)||(xe>maxx)
    disp('Coordinates out of bounds in UpdateNewSeeds!')
end

newSeeds = zeros(ye-ys+1,xe-xs+1,cc_max);
%%% ignore the immigration in the codes ********
if doDispersal
    newSeeds(:,:,:) = seedrain(:,:,:);
    seedrain(:,:,:) = 0;   %% verynewseeds is refreshed.
end

end




function [newSeeds] = Seeds_From_OutSim_Region(cc_max,ys,ye,xs,xe,doDispersal,outerseeds,newSeeds)
%%%%%%%%%%%%%%%%%%%%%%%
%%% Add a given number of seeds to all grid cell in seedrain_newSeeds
%%%%%%%%%%%%%%%%%%%%%%%
%%% Inputs
%%% cc_max [] the maximum species in one cell
%%% ys, ye, xs, xe [] the start and end of coordinates (row and column) (relative) of current grid subarea
%%% doDispersal [] if explicit dispersal has to be done done; normally the case; corresponds to idisp
%%% outerseed [INTEGER] the number of seeds/patch entering cells from outside the total simulation area 
%%% newSeeds [] finally, arriving seeds of all species in each cell (the entire grid) - size(maxy,maxx,Ccrown)
%%% Outputs
%%% newSeeds [] finally, arriving seeds of all species in each cell (the entire grid)

%%% ignore the immigration *******
if doDispersal
    if outerseeds>0   %% outerspace boundary conditions
        newSeeds(:,:,:) = newSeeds(:,:,:) + outerseeds;
    end
end

end








