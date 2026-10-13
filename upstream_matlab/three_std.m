%modas search box cal std
% load('/data2/winter/obs_noarc.mat')
clearvars -except ori AOTD
lat_aim = 60:0.25:90;
lon_aim = -180:0.25:180-0.25;

mon = [AOTD.mo];
t_upper = nan(4,121,1440,102);
t_lower = nan(4,121,1440,102);
s_upper = nan(4,121,1440,102);
s_lower = nan(4,121,1440,102);
parpool(8);
for T = 1:4
    tic
    if T == 1
        season = [4,5,6];
    elseif T==2
        season = [7,8,9];
    elseif T==3
        season = [10,11,12];
    else
        season = [1,2,3];
    end
    ind0 = ismember(mon,season);
    ori0 = AOTD(ind0);
    lat0 = [ori0.lat];
    lon0 = [ori0.lon];
    for i = 1:length(lat_aim)
        lat = lat_aim(i);
%         Lx = 2*120000/(0.35*lat*lat+300);
%         Ly = 2*120000/(0.35*lat*lat+400);
%      
        parfor j = 1:length(lon_aim)
            lon = lon_aim(j);
            if lon==-180
                lonc = lon*ones(1,length(lon0));
                lonc(lon0>=0) = -lonc(lon0>=0);
                del_y = abs(lon0-lonc);
            else
                del_y = abs(lon0-lon);
            end
            del_x = abs(lat0-lat);
            dis = abs(del_x+del_y*1i);
            idx = dis<=1;
            
%             distx = abs(distance(lat,lon,lat0,lon,6371));
%             disty = abs(distance(lat,lon,lat,lon0,6371));
%             idx = (distx<=Lx)&(disty<=Ly);
            if sum(idx) >= 1
                orii = ori0(idx);
                t = [orii.temp];
                s = [orii.salt];
                for k = 2%1:102
                    tk = t(k,:);
                    sk = s(k,:);
                    if (sum(~isnan(tk))>=5)&(std(~isnan(tk))~=0)
                        tk = rmoutliers(tk); 
                        t_upper(T,i,j,k) = nanmean(tk) + 3*std(tk(~isnan(tk)));
                        t_lower(T,i,j,k) = nanmean(tk) - 3*std(tk(~isnan(tk)));
                        %                     ind = (t>=t_lower)&(t<=t_upper);
                        %                     tk(~ind) = nan;
                    end
                    if (sum(~isnan(sk))>=5)&(std(~isnan(sk))~=0)
                        sk = rmoutliers(sk(~isnan(sk))); 
                        s_upper(T,i,j,k) = nanmean(sk) + 3*std(sk(~isnan(sk)));
                        s_lower(T,i,j,k) = nanmean(sk) - 3*std(sk(~isnan(sk)));
                        %                     ind = (s>=s_lower)&(s<=s_upper);
                        %                     sk(~ind) = nan;
                    end
                end
            end
        end
    end
    toc
end