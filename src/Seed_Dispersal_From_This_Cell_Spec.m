%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Subfunction Seed dispersal without FFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% stockability according to the land cover
function [seedrain,ksv,curDisp] = Seed_Dispersal_From_This_Cell_Spec(isp,y,x,stockability,spec,maxy,maxx,seedrain,seedprod,ksv,bareId, ...
                                            cc_max,doStochSeedDisp,values)
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Do the explicit seed dispersal for one species from one cell (y,x) to all sink cells
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Inputs
%%% isp [] index of species
%%% x, y [] coordinates of current grid cell
%%% stockability [] whether this grid cell can produce seeds or not, e.g., impervious surface
%%% spec [] the properties of this species
%%% maxy, maxx [] the number of rows and columns of simulation domain
%%% seedrain [] store the number of very new seeds - size(maxy,maxx,Ccrown)
%%% seedprod [] seed production - size(maxy,maxx,Ccrown)
%%% ksv [] land cover in the simulation area - size(maxy,maxx)
%%% bareId [] flag - bare land
%%% cc_max [] the maximum of the vegetation type
%%% doStochSeedDisp [] whether do stochastic seed dispersal or not [true or false]
%%% values [] the priority of vegetation, e.g., [1 1 2:26]
%%% Outputs
%%% seedrain [] store the number of very new seeds - size(maxy,maxx,Ccrown)
%%% ksv [] modified land cover in the simulation area - size(maxy,maxx)
%%% curDisp [] store seed dispersal just from this cell

producedseeds = seedprod(y,x,isp);     %%% need to change *************, and seedrain **********

curDisp = zeros(maxy,maxx);

if stockability(y,x)==0  %% if cell is not stockable, no seed dispersal from this cell, e.g., urban
    return
end
radi = spec{isp}.radi;
if producedseeds>0
    for ywin=-radi:radi
        for xwin=-radi:radi
            wy = y + ywin;  % wy = y - ywin;
            wx = x + xwin;  % wx = x - xwin;
            if (wy<1) || (wy>maxy) || (wx<1) || (wx>maxx)
                continue
            else
                %%% what happens if seeds cross edge of simulation domain? Check whether this sink cell is allowed by the boundary conditions
                [dodisp,wx,wy] = Boundary_Conditions(wy,maxy,wx,maxx,'a');
                if dodisp
                    if stockability(wy,wx)>0 %%% if sink cell is stockable
                        prob = spec{isp}.kernel(radi+1+ywin,radi+1+xwin); %%% probability of dispersal between source and sink cell
                        if doStochSeedDisp        %% Stochastic dispersal (Poisson distribution)
                            seedsGround = Draw_Seeds_From_Poisson(producedseeds,prob); %Uniform, Poisson
                            contrib = seedsGround;  
                        else                      %% Deterministic dispersal, probability interpreted as fraction
                            contrib = prob*producedseeds;
                        end
                        curDisp(wy,wx) = contrib;   %% store seed dispersal just from this cell
                        seedrain(wy,wx,isp) = seedrain(wy,wx,isp) + contrib;

                        if seedrain(wy,wx,isp)<=1e-2
                            seedrain(wy,wx,isp) = 0;
                            curDisp(wy,wx) = 0;
					    end

						%%% convert land use to seedling
						if seedrain(wy,wx,isp)>1e-2 %%% The sink cell has seeds
                            %%% lower number = higher priority
                            priority = containers.Map({'A','B','C','D','E','F','G','H','I','G','K',...
                                'L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','-'}, values);
                            %%% skip if ksv{wy,wx} already has a non-'-' value
                            if ksv{wy,wx}(isp)~='-'
                                continue;  %ksv{wy,wx}(isp) = ksv{wy,wx}(isp);
                            end 
                            %%% skip if ksv{y,x} has a lowercase or '-'
                            if ~isstrprop(ksv{y,x}(isp),'upper')
                                continue;  %ksv{wy,wx}(isp) = ksv{wy,wx}(isp);
                            end
                            %%% check the priority: only update if no higher or equal priority in ksv{wy,wx}
                            update = true;
                            for j=1:cc_max
                                bj = ksv{wy,wx}(j);
                                if isstrprop(bj,'upper') || bj=='-'
                                    if bj~='-' && priority(bj)<=priority(ksv{y,x}(isp))
                                        update = false;
                                        break;
                                    end
                                end
                            end
                            %%% update 
                            if update
                                ksv{wy,wx}(isp) = lower(ksv{y,x}(isp));
                            end
						end
                    end
                end
            end
        end
    end
end

end


function [dodisp, wx, wy] = Boundary_Conditions(wy,maxy,wx,maxx,boundaries)
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Do the explicit seed dispersal for one species from one cell (y,x) to
%%% all sink cells
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Inputs
%%% wy, wx [] absolute coordinates of sink cell, potentially modified 
%%% maxy, maxx [] number of rows and columns of simulation domain
%%% boundaries [] boundary condition, absorbing (a), cyclic (c),
%%%               half-cyclic south-north (s), half cyclic west-east (w)
%%% Outputs
%%% dodisp [] do dispersal into the originally determined sink cell
%%% wy, wx [] 

dodisp = true;

switch boundaries
    case 'a'
        if (wx<1)||(wx>maxx)||(wy<1)||(wy>maxy)
            dodisp = false;
        end
    case 'w'
        if wx>=0
            wx = mod(wx,maxx);
        else
            wx = maxx - mod(abs(wx),maxx);
        end
        if(wx==0)
            wx = maxx;
        end
        if (wy<1)||(wy>maxy)
            dodisp = false;
        end
    case 's'
        if wy>=0
            wy = mod(wy,maxy);
        else
            wy = maxy - mod(abs(wy),maxy);
        end
        if (wy==0)
            wy = maxy;
        end
        if (wx<1)||(wx>maxx)
            dodisp=false;
        end
    case 'c'
        if wy>=0
            wy = mod(wy,maxy);
        else
            wy = maxy - mod(abs(wy),maxy);
        end
        if (wy==0)
            wy = maxy;
        end
        if wx>=0
            wx = mod(wx,maxx);
        else
            wx = maxx - mod(abs(wx),maxx);
        end
        if(wx==0)
            wx = maxx;
        end
    otherwise
        disp('Wrong boundary conditions! Choose one of a, c, s, w!')
        if (wx<1)||(wx>maxx)||(wy<1)||(wy>maxy)
            dodisp = false;
        end
end
end


function [ms] = Draw_Seeds_From_Poisson(n,prob)
%%%%%%%%%%%%%%%%%%%%%%%%
%%% Samples from a poisson distribution
%%%      Since the expected number of seeds deposited at a sink cell (ms) is small
%%%      relative to the seed source, it is reasonable to assume that "ms" is
%%%      Poisson distributed.
%%%%%%%%%%%%%%%%%%%%%%%%
%%% Inputs
%%% n [] seeds emitted from the source cell
%%% prob [] probability
%%% Outputs
%%% ms [] seeds deposited at a sink cell


if prob>=1
    ms = n;
    return
end
if (prob<=0)||(n<=0)
    ms = 0;
    return
end

%%%%
if n<1
    n = 1;
end

%%%%%
scale = 1;  %%% As Poisson distribution is designed for integer, carbon pool of flowers is around 5-20.
n1 = n*scale;

randomvalue = rand;

lambda = prob*n1;
distsum = exp(-lambda);

for i=1:n1
    m = i;
    if distsum>=randomvalue || distsum>=1
        break
    end
    distsum = distsum + poisspdf(m,lambda); 
end
ms = m-1;
ms = ms/scale;

end


function [randSeeds] = Draw_Seeds_From_Uniform(number,prob)
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Samples from a uniformaly distributed random number
%%% U(0,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Inputs
%%% number [] sample size
%%% prob [] probability
%%% Outputs
%%% randSeeds [] drawn value

if prob>=1
    randSeeds = number;
    return
end
if (prob<=0)||(number<=0)
    randSeeds = 0;
    return
end

randValue = rand;
if prob>randValue
    randSeeds = number*prob;
else
    randSeeds = 0;
end

end


% function [ix] = Draw_From_Dis_Tri(n,prob)
% %%%%%%%%%%%%%%%%%%%%%%%%
% %%% Samples from a binomial distribution
% %%% 1. for large sample size, the binomial distribution is approximated by
% %%%       a. a Poisson distribution, if small probability
% %%%       b. a normal distribution, if large probability
% %%% 2. otherwise the binomial distribution is explicitly calculated
% %%%%%%%%%%%%%%%%%%%%%%%%
% %%% Inputs
% %%% n [] sample size
% %%% prob [] probability
% %%% Outputs
% %%% ix [] drawn value
% 
% 
% if prob>=1
%     ix = n;
%     return
% end
% if (prob<=0)||(n<=0)
%     ix = 0;
%     return
% end
% 
% randomvalue = rand;
% startdist = (1-prob)^n;
% if (startdist==0)||(startdist==1.0)||(n>50)
%     mu = prob*n;
%     if mu<20   %% for small "prob", a Poisson distribution
%         distsum = exp(-mu);
%         for i=1:10000
%             ri = i;
%             if (distsum>=randomvalue)||(distsum>=1)
%                 break
%             end
%             distsum = distsum + poisspdf(ri,mu);
%         end
%         ix = ri-1;
%     else   %% otherwise a normal distribution
%         sigma = sqrt(mu);
%         ix = Draw_from_Normal_Dist(mu,sigma,randomvalue);
%     end
% end
% 
% end
% 
% 
% function [rx] = Draw_from_Normal_Dist(mu,sigma,randomvalue)
% %%%%%%%%%%%%%%%%%%%%%%%%
% %%% Samples from a normal distribution
% %%%%%%%%%%%%%%%%%%%%%%%%
% %%% Inputs
% %%% mu [] mean
% %%% sigma [] standard deviation
% %%% randomvalue [] random value
% %%% Outputs
% %%% rx [] drawn value
% 
% %%% First, initialize the standard normal distribution (mu=0,sigma=1) array
% %%% in the range of -5 to +5. It has 202 elements. Element 100 corresponds
% %%% to the median and mean, namely 1 and 202 correspond to -5 and 5, respectively.
% nd = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0000106885, 0.0000133457, ...
%            0.0000166238, 0.0000206575, 0.0000256088, 0.0000316712, 0.0000390756, 0.0000480963, 0.0000590589, ...
%            0.000072348, 0.0000884173, 0.0001078, 0.00013112, 0.000159109, 0.000192616, 0.000232629, ...
%            0.000280293, 0.000336929, 0.000404058, 0.000483424, 0.000577025, 0.000687138, 0.000816352, ...
%            0.000967603, 0.00114421, 0.0013499, 0.00158887, 0.00186581, 0.00218596, 0.00255513, 0.00297976, ...
%            0.00346697, 0.00402459, 0.00466119, 0.00538615, 0.00620967, 0.00714281, 0.00819754, 0.00938671, ...
%            0.0107241, 0.0122245, 0.0139034, 0.0157776, 0.0178644, 0.0201822, 0.0227501, 0.0255881, 0.0287166, ...
%            0.0321568, 0.0359303, 0.0400592, 0.0445655, 0.0494715, 0.0547993, 0.0605708, 0.0668072, 0.0735293, ...
%            0.0807567, 0.088508, 0.0968005, 0.10565, 0.11507, 0.125072, 0.135666, 0.146859, 0.158655, 0.171056, ...
%            0.18406, 0.197663, 0.211855, 0.226627, 0.241964, 0.257846, 0.274253, 0.29116, 0.308538, 0.326355, ...
%            0.344578, 0.363169, 0.382089, 0.401294, 0.42074, 0.440382, 0.460172, 0.480061, 0.5, 0.519939, ...
%            0.539828, 0.559618, 0.57926, 0.598706, 0.617911, 0.636831, 0.655422, 0.673645, 0.691462, 0.70884, ...
%            0.725747, 0.742154, 0.758036, 0.773373, 0.788145, 0.802337, 0.81594, 0.828944, 0.841345, 0.853141, ...
%            0.864334, 0.874928, 0.88493, 0.89435, 0.9032, 0.911492, 0.919243, 0.926471, 0.933193, 0.939429, ...
%            0.945201, 0.950529, 0.955435, 0.959941, 0.96407, 0.967843, 0.971283, 0.974412, 0.97725, 0.979818, ...
%            0.982136, 0.984222, 0.986097, 0.987776, 0.989276, 0.990613, 0.991802, 0.992857, 0.99379, 0.994614, ...
%            0.995339, 0.995975, 0.996533, 0.99702, 0.997445, 0.997814, 0.998134, 0.998411, 0.99865, 0.998856, ...
%            0.999032, 0.999184, 0.999313, 0.999423, 0.999517, 0.999596, 0.999663, 0.99972, 0.999767, 0.999807, ...
%            0.999841, 0.999869, 0.999892, 0.999912, 0.999928, 0.999941, 0.999952, 0.999961, 0.999968, 0.999974, ...
%            0.999979, 0.999983, 0.999987, 0.999989, 0.999991, 0.999993, 0.999995, 0.999996, 0.999997, 0.999997, ...
%            0.999998, 0.999998, 0.999999, 0.999999, 0.999999, 0.999999, 1.0, 1.0, 1.0, 1.0];
% %%% if random value=0, then a value at the very left tail of the distribution is chosen
% if (randomvalue == 0.0)
%     rx = -5.0 * s + mu;
%     return
% end
% %%% 0.5 is the median of the distribution. Therefore for smaller random values the first half, for larger random values the second half of the array can be searched
% if (randomvalue>0.5)
%     istart = 101;
% else
%     istart = 1;
% end 
% %%% search in the table where the random value belongs to 
% for i=istart:istart+101
%     if(nd(i)>randomvalue)
%         index = i;
%         break
%     end
% end
% 
% f1 = nd(index-1);
% f2 = nd(index);
% %%% transform the indices from 1 to 200 to (0,1) normal distribution
% v1 = (index - 1)/20.0 - 5.0;
% v2 = (index)/20.0 - 5.0;
% xtransf = v1 + (v2 - v1)*(randomvalue - f1)/(f1 - f2);
% 
% rx = xtransf*sigma + mu;
% end

% function [ms] = Draw_Seeds_From_Poisson(n,prob)
% %%%%%%%%%%%%%%%%%%%%%%%%
% %%% Samples from a poisson distribution
% %%%      Since the expected number of seeds deposited at a sink cell (ms) is small
% %%%      relative to the seed source, it is reasonable to assume that "ms" is
% %%%      Poisson distributed.
% %%%%%%%%%%%%%%%%%%%%%%%%
% %%% Inputs
% %%% n [] seeds emitted from the source cell
% %%% prob [] probability
% %%% Outputs
% %%% ms [] seeds deposited at a sink cell
% 
% 
% if prob>=1
%     ms = n;
%     return
% end
% if (prob<=0)||(n<=0)
%     ms = 0;
%     return
% end
% 
% random1 = rand;  %% this to randomly select which grid may have seeds.
% if prob>random1
%     %%%%%
%     scale = 100;  %% As Poisson distribution is designed for integer, carbon pool of flowers is around 5-20.
%     n1 = n*scale;
% 
%     randomvalue = rand;  %% this is to see how many seeds this grid may have, if this grid can have seeds.
% 
%     lambda = prob*n1;
%     distsum = exp(-lambda);
% 
%     for i=1:n1
%         m = i;
%         if distsum>=randomvalue || distsum>=1
%             break
%         end
%         distsum = distsum + poisspdf(m,lambda); 
%     end
%     ms = m-1;
%     ms = ms/scale;
% else
%     ms = 0;
% end
% 
% end













