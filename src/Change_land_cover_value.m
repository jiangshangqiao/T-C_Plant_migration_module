%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Subfunction - change the land cover value and fraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ksv,Ccrown_t] = Change_land_cover_value(curYsta,curYend,curXsta,curXend,ksv,...
                                        bareId,hc,hctm1,Parameters_SD,Ccrown_t,cc_max,stockability, ...
                                        Datam_S,values)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Change land cover value for each time step
%%% Seedling: a, b, c ...
%%% Maturity: 1, 2, 3 ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Input data
%%% values [] the priority of vegetation, e.g., [1 1 2:26]
%%% Output data
%%% Ccrown_t [] the changed fraction

fraThre = 0.95;  %% the threshold to change the value

maxy = curYend-curYsta+1;  
maxx = curXend-curXsta+1; 
% Ccrown_t = reshape(Ccrown_t,maxy,maxx,cc_max);

for x=curXsta:curXend
    for y=curYsta:curYend
        if stockability(y,x)>1e-4
            % if all(hc(y,x,:)<1e-4)
            %     %%% if all hc<1e-4, this cell is set to bare land.
            %     ksv(y,x) = bareId;
            % end
            [Ccrown,~,~,~,hcM,~,~,~,~] = Parameters_SD(ksv(y,x));   %%% Ccrown in last time step 
            %%% change the land cover value
            for isp=1:length(Ccrown)
                if hc(y,x,isp)>hcM(isp) && Ccrown_t(y,x,isp)>0
                    lc = ksv{y,x}(isp);  %%% get the current vegetation species 
                    if lc~='-'
                        ksv{y,x}(isp) = upper(lc);
                    end
                end
                %%% For shrub and tree
                %%% when hc is decreasing, this plant may die.
                if ((hc(y,x,isp)<0.001)&&(hctm1(y,x,isp)>=0.001)) 
                    ksv{y,x}(isp) = '-';
                end
                %%% For grass, shrub and tree
                %%% In summer, if hc<0.006, this plant die.
                if (Datam_S(2)==7) && (Datam_S(3)==2) && (hc(y,x,isp)<0.006)
                    ksv{y,x}(isp) = '-';
                end
            end
            
            priority = containers.Map({'A','B','C','D','E','F','G','H','I','G','K',...
                                'L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','-'}, values);
            %%%%%%%% change the vegetation fraction %%%%%%%%%%%%%%%%
            cum_sum = 0;
            Ccrown_new = zeros(size(Ccrown_t(y,x,:)));
            %%% group indices by priority level
            [~,~,grp] = unique(values(1:cc_max),'stable');
            groups  = cell(1,max(grp));
            for g=1:max(grp)
                groups{g} = find(grp==g);
            end
            %%% fraction change
            for g=1:length(groups)
                idx = groups{g};
                group_sum = sum(Ccrown_t(y,x,idx),'all');
                if cum_sum+group_sum>=fraThre
                    remaining = 1 - cum_sum;
                    %%% check if any single is >=fraThre
                    big_enough = Ccrown_t(y,x,idx)>=fraThre;
                    if any(big_enough)
                        %%% allocate all to the first such member
                        % j = idx(find(big_enough,1,'first'));
                        % Ccrown_new(j) = remaining;
                        j95 = idx(find(big_enough,1,'first'));
                        j50 = idx(Ccrown_t(y,x,idx)>0.5);
                        j50 = j50(j50~=j95);
                        Ccrown_t(y,x,j50) = 0.5;
                        group_sum = sum(Ccrown_t(y,x,idx),'all');
                        Ccrown_new(idx) = Ccrown_t(y,x,idx)/group_sum * remaining;
                    else
                        %%% proportional allocation
                        Ccrown_new(idx) = Ccrown_t(y,x,idx)/group_sum * remaining;
                    end
                    break
                else
                    Ccrown_new(idx) = Ccrown_t(y,x,idx);
                    cum_sum = cum_sum + group_sum;
                end
            end
            Ccrown_t(y,x,:) = Ccrown_new;

            %%%%%%%% change the land cover value %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Fc > fraThre, just keep the first few vegetation types
            cum_f = 0;
            for j=1:cc_max
                frac = Ccrown_t(y,x,j);
                pj = priority(upper(ksv{y,x}(j)));
                cum_f = cum_f + frac;
                %%% if equal priority exists, alone vegetation can reach fraThre, just keep this vegetation in case this vegetation is not the first one.
                if frac>=fraThre
                    pr = arrayfun(@(c) priority(c), upper(ksv{y,x}));
                    idx = find((pr>=pj) & ((1:cc_max)~=j));
                    ksv{y,x}(idx) = '-';
                    break
                elseif cum_f>=fraThre
                    ksv{y,x}(j+1:end) = '-';
                    break;
                end
            end

        end
    end
end

end