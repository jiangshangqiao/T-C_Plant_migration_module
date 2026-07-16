%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Main parameters %%%%%
%%% Just use high vegetation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Ccrown,spec,ZR95_H,ZR95_L,hcM,SAI_H,SAI_L,TBio_H,TBio_L] = Parameters_SD_2D2(ANSWER)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% PARAMETERS SOILS AND VEGETATION  %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ANSWER;
%%%% 100 Forest seedling
%%%%%%%%%%%%%%%%%%%%%%%%%%
switch upper(ANSWER)
    case {'A--'}
        %%%% land cover partition  Forest (Larix decidua -- a deciduous coniferous tree)
        Cwat = 0; Curb = 0.0 ; Crock = 0.0;
        Cbare = 0.0; Ccrown = [0.25 0.0 0.0];
        cc = length(Ccrown);%% Crown area
        II = [1 1 0 1 0]>0;
        III_H = [1 0 0];
        III_L = [0 0 0];
        n = 1;
    case {'-B-'}
        %%%% land cover partition  Eve. shrub (Rhododendron ferrugineum)
        Cwat = 0; Curb = 0.0 ; Crock = 0.0;
        Cbare = 0.0; Ccrown = [0.0 0.25 0.0];
        cc = length(Ccrown);%% Crown area
        II = [1 1 0 1 0]>0;
        III_H = [0 0 0];
        III_L = [0 1 0];
        n = 2;
    case {'3'}
        %%%% land cover partition  Deci. shrub (Vaccinium myrtillus)
        Cwat = 0; Curb = 0.0 ; Crock = 0.0;
        Cbare = 0.0; Ccrown = [0.25];
        cc = length(Ccrown);%% Crown area
        II = [0 0 1 0 0]>0;
        n = 3;
    case {'--C'}
        %%%% land cover partition  Grass (Poa alpina -- C3)
        Cwat = 0; Curb = 0.0 ; Crock = 0.0;
        Cbare = 0.0; Ccrown = [0.0 0.0 0.25];
        cc = length(Ccrown);%% Crown area
        II = [1 1 0 1 0]>0;
        III_H = [0 0 0];
        III_L = [0 0 1];
        n = 4;
    case '---'
        %%%% land cover partition  Bare land
        Cwat = 0; Curb = 0.0 ; Crock = 0.0;
        Cbare = 1.0; Ccrown = [0.0 0.0 0.0];
        cc = length(Ccrown);%% Crown area
        II = [1 1 0 1 0]>0;
        III_H = [0 0 0];
        III_L = [0 0 0];
        n = 5;
    case {'AB-'}
        %%% tree/tree seedling + shrub/shrub seedling
        Cwat = 0; Curb = 0.0 ; Crock = 0.0;
        Cbare = 0.0; Ccrown = [0.25 0.25 0.0];
        cc = length(Ccrown);%% Crown area
        II = [1 1 0 1 0]>0;
        III_H = [1 0 0];
        III_L = [0 1 0];
        n = 6;
    case {'A-C'}
        %%% tree/tree seedling + grass/grass seedling
        Cwat = 0; Curb = 0.0 ; Crock = 0.0;
        Cbare = 0.0; Ccrown = [0.25 0.0 0.25];
        cc = length(Ccrown);%% Crown area
        II = [1 1 0 1 0]>0;
        III_H = [1 0 0];
        III_L = [0 0 1];
        n = 7;
    case {'-BC'}
        %%% shrub/shrub seedling + grass/grass seedling
        Cwat = 0; Curb = 0.0 ; Crock = 0.0;
        Cbare = 0.0; Ccrown = [0.0 0.25 0.25];
        cc = length(Ccrown);%% Crown area
        II = [1 1 0 1 0]>0;
        III_H = [0 0 0];
        III_L = [0 1 1];
        n = 8;
    case {'ABC'}
        %%% tree/tree seelding + shrub/shrub seedling + grass/grass seedling
        Cwat = 0; Curb = 0.0 ; Crock = 0.0;
        Cbare = 0.0; Ccrown = [0.25 0.25 0.25];
        cc = length(Ccrown);%% Crown area
        II = [1 1 0 1 0]>0;
        III_H = [1 0 0];
        III_L = [0 1 1];
        n = 9;
    otherwise
        disp('INDEX FOR SOIL VEGETATION PARAMETER INCONSISTENT')
        return
end
%%%%%%%%%%%%%%
%%%%% root depth
ZR95_H = [800 0 0];   %%[mm]
ZR95_L = [0 800 250];   %%[mm]
ZR95_H = ZR95_H.*III_H; ZR95_L = ZR95_L.*III_L; 

%%%%% mature height
hcM = [10 0.48 0 0.005 0];
hcM = hcM(II);

%%%% SAI
SAI_H = [0.1 0 0 0 0];
SAI_L = [0 0.1 0.1 0.001 0];
SAI_H = SAI_H(II); SAI_L = SAI_L(II);

TBio_H = [20 0 0 0 0]; 
TBio_L = [0 1 1 1 0];
TBio_H = TBio_H(II); TBio_L = TBio_L(II);
%%%% species properties
spec = cell(5,1);  %% *******************
spec{1}.dispfrac = 0.99;
spec{1}.alpha1 = 25;
spec{1}.alpha2 = 200;
spec{2}.dispfrac = 0.99;
spec{2}.alpha1 = 10;
spec{2}.alpha2 = 50;
spec{3}.dispfrac = 0;
spec{3}.alpha1 = 0;
spec{3}.alpha2 = 0;
spec{4}.dispfrac = 0.99;
spec{4}.alpha1 = 10;
spec{4}.alpha2 = 50;
spec{5}.dispfrac = 0;
spec{5}.alpha1 = 0;
spec{5}.alpha2 = 0;
spec = spec(II);

end



