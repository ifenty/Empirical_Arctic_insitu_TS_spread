clearvars -except ori AOTD
% load('../obs_noarc_v2.mat')
lat_aim = 60:0.25:90;
lon_aim = -180:0.25:180-0.25;

mon = [AOTD.mo];
t_upper = nan(4,121,1440,102);
t_lower = nan(4,121,1440,102);
s_upper = nan(4,121,1440,102);
s_lower = nan(4,121,1440,102);
load('/data2/winter/ts_range_v5.mat')
% parpool(16);
fi = 1;
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
    ind0 = find(ismember(mon,season));
    ori0 = AOTD(ind0);
    latt = [ori0.lat];
    lonn = [ori0.lon];
    for i = 1:length(ori0)
        lat0 = latt(i);
        lon0 = lonn(i);
        distx = distance(lat0,lon0,lat_aim,lon0,6371);
        disty = distance(lat0,lon0,lat0,lon_aim,6371);
        indx = find(abs(distx) == min(abs(distx)));
        indy = find(abs(disty) == min(abs(disty)));
        indx = indx(1);
        indy = indy(1);
        t = ori0(i).temp;
        s = ori0(i).salt;
        if isempty(t)||isempty(s)
            continue
        end
        t1 = nan(102,1);
        s1 = nan(102,1);
        for k = 1:102
            tk = t(k,1);
            sk = s(k,1);
            t_range_lower = squeeze(t_lower(T,indx,indy,k));
            t_range_upper = squeeze(t_upper(T,indx,indy,k));
            s_range_lower = squeeze(s_lower(T,indx,indy,k));
            s_range_upper = squeeze(s_upper(T,indx,indy,k));
            if ~isnan(t_range_lower)
                if ~isnan(tk)
                    ind0 = (tk>=t_range_lower)&(tk<=t_range_upper);
                    if ind0
                        t1(k,1) = tk;
                    else
                        t1(k,1) = nan;
                    end
                end
            else
                t1(k,1) = tk; 
            end
            if ~isnan(s_range_lower)
                if ~isnan(sk)
                    ind0 = (sk>=s_range_lower)&(sk<=s_range_upper);
                    if ind0
                        s1(k,1) = sk;
                    else
                        s1(k,1) = nan;
                    end
                end
            else
                s1(k,1) = sk;
            end
        end
%         AOTD(fi) = ori0(i);
%         AOTD(fi).temp = t1;
%         AOTD(fi).salt = s1;
        AOTD1(fi).lat  = ori0(i).lat;
        AOTD1(fi).lon  = ori0(i).lon;
        AOTD1(fi).ye   = ori0(i).ye;
        AOTD1(fi).mo   = ori0(i).mo;
        AOTD1(fi).da   = ori0(i).da;
        AOTD1(fi).ho   = ori0(i).ho;
        AOTD1(fi).mi   = ori0(i).mi;
        AOTD1(fi).se   = ori0(i).se;
        AOTD1(fi).source = ori0(i).source;
        AOTD1(fi).temp = t1;
        AOTD1(fi).salt = s1;
        fi = fi+1;
    end
    toc
end
