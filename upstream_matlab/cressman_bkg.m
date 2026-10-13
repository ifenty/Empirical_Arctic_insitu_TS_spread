% load('G:\���ݼ�\����\��ֵ�����ݼ�\arctic_data.mat')
% load('/data2/ljl/master_arctic/arctic_data_mon.mat')
clearvars -except arctic_data arctic_data_climatology_qc observ_mon_final dataset_mon data_review ori AOTD
MONTH = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', ...
    'Sep', 'Oct', 'Nov', 'Dec'};
%% ʱ�ռ�Ȩƽ������̬
% ori = data_review;
lon = -180:0.25:179.75;
lat = 60:0.25:90;
[LAT_AIM,LON_AIM] = meshgrid(lat,lon);
Lt = 1000;
load('/data2/ljl/master_arctic/topographic_land.mat')
% load('std_data.mat')
% a = dataset_mon.May;
% depth = a(1).depth;
depth = AOTD.pres;
mkdir('/data2/ljl//master_arctic//barnes_test9');
% load('topographic_land.mat');
% idx_land = isnan(MASK);
fi = 0;
% arctic_data_1 = [];
% for i = 1:12
%     obs = dataset_mon.(MONTH{i});
%     arctic_data_1 = [arctic_data_1,obs];
% end
arctic_data_1 = AOTD;
climatology_temp = nan(1440,121,57);
climatology_salt = nan(1440,121,57);
% %     idx_TT = find(mon == TT);
%     arctic_data_1 = observ_mon_final.(MONTH{TT});
lon_obs = [arctic_data_1.lon];
lat_obs = [arctic_data_1.lat];
temp = [arctic_data_1.temp];
salt = [arctic_data_1.salt];
fi = fi+1;
file_name = sprintf('/data2//ljl//master_arctic//barnes_test9//climatology_annual.nc');
nccreate(file_name, 'geolon', 'Dimensions', {'lon', 1440,'lat',121});
nccreate(file_name, 'geolat', 'Dimensions', {'lon', 1440,'lat',121});
nccreate(file_name, 'temp', ...
    'Dimensions', {'lon', 1440, 'lat', 121, 'depth', 57});
nccreate(file_name, 'salt', ...
    'Dimensions', {'lon', 1440, 'lat', 121, 'depth', 57});
nccreate(file_name, 'depth', ...
    'Dimensions', {'depth', 57});
% parpool(32)
for i = 1:57
    tic
    t = temp(i,:);
    %                 if ~all(isnan(t))
    s = salt(i,:);
    idx_ii = find(~isnan(t)|~isnan(s));
    arctic_data_11 = arctic_data_1(idx_ii);
    lonn = [arctic_data_11.lon];
    latt = [arctic_data_11.lat];
    t = t(idx_ii);
    s = s(idx_ii);
%     [t,s] = obs_modas(latt,lonn,t,s,LAT_AIM,LON_AIM);
%     idx_land = MASK==-1;
%     t(idx_land) = nan;
%     s(idx_land) = nan;
    [t,s] = obs_merge(latt,lonn,t,s,LAT_AIM,LON_AIM);
    idx_land = MASK==-1;
    t(idx_land) = nan;
    s(idx_land) = nan;
    z1 = Cressman_bkg(LON_AIM,LAT_AIM,t,MASK,depth(i));
%     h = fspecial('average', [3 3]);
    idx = isnan(z1);
    z1(idx)=nan;
%     z1 = Barnes_sc(LON_AIM,LAT_AIM,t,zo,2,MASK,depth(i));
    z2 = Cressman_bkg(LON_AIM,LAT_AIM,s,MASK,depth(i));
    idx = isnan(z2);
    z2(idx)=nan;
    climatology_temp(:,:,i) = z1;
    climatology_salt(:,:,i) = z2;
    toc
end
ncwrite(file_name,'geolon',LON_AIM);
ncwrite(file_name,'geolat',LAT_AIM);
ncwrite(file_name, 'temp', climatology_temp);
ncwrite(file_name, 'salt', climatology_salt);
ncwrite(file_name, 'depth', depth(1:57));




function [t,s] = obs_modas(latt,lonn,t,s,LAT_AIM,LON_AIM)
temp_merge = nan(1440,121);
salt_merge = nan(1440,121);
% mon = month(time(1));
% years = year(time);
% timm = datetime(1998,7,1,12,00,00)*nan(length(years),1);
parfor i = 1:1440
    for j = 1:121
        lon = LON_AIM(i,j);
        lat = LAT_AIM(i,j);
        k=2;
        d = distance(lat,lon,latt,lonn,6371);     
%         for k = 1:2
            Lx = k*120000/(0.35*lat*lat+300);
            Ly = k*120000/(0.35*lat*lat+400);
            Lt = 1000;
            L = sqrt(Lx*Lx+Ly*Ly);
            idx = d<=L;
            %             if sum(idx)<100&&k==1
            %                 continue
            %             else
            t_idx = t(idx);
            s_idx = s(idx);
            idx_t = ~isnan(t_idx);
            idx_s = ~isnan(s_idx);
            t_idx = t_idx(idx_t);
            s_idx = s_idx(idx_s);
            %                 time_idx = time(idx)';
            lat_idx = latt(idx);
            lat_idx_t = lat_idx(idx_t);
            lat_idx_s = lat_idx(idx_s);
            lon_idx = lonn(idx);
            lon_idx_t = lon_idx(idx_t);
            lon_idx_s = lon_idx(idx_s);
            %                 timmm = timm(idx)';
            bt = exp(-power((lat-lat_idx_t)/Lx,2)-power((lon-lon_idx_t)/Ly,2));
            bs = exp(-power((lat-lat_idx_s)/Lx,2)-power((lon-lon_idx_s)/Ly,2));
            temp_merge(i,j) = sum(bt.*t_idx)/sum(bt);
            salt_merge(i,j) = sum(bs.*s_idx)/sum(bs);
            %             end
%             t_idx = t(idx);
%             s_idx = s(idx);
%             lat_idx = latt(idx);
%             lon_idx = lonn(idx);
%             b = exp(-power((lat-lat_idx)/Lx,2)-power((lon-lon_idx)/Ly,2));
%             temp_merge(i,j) = sum(b.*t_idx)/sum(b);
%             salt_merge(i,j) = sum(b.*s_idx)/sum(b);
%         end

    end
end
t = temp_merge;
s = salt_merge;
end





function [t,s] = obs_merge(latt,lonn,t,s,LAT_AIM,LON_AIM)
temp_merge = nan(1440,121);
salt_merge = nan(1440,121);
parfor i = 1:1440
    for j = 1:121
        lon = LON_AIM(i,j);
        lat = LAT_AIM(i,j);
        del_x = abs(latt-lat);
        % %         del_y = min(abs(lonn-lon-180),abs(lonn-lon+180),abs(lonn-lon))
        if lon==-180
            lonc = lon*ones(1,length(lonn));
            lonc(lonn>=0) = -lonc(lonn>=0);
            del_y = abs(lonn-lonc);
        else
            del_y = abs(lonn-lon);
        end
        dis = abs(del_x+del_y*1i);
        temp_merge(i,j) = nanmean(t(dis<=0.25));
        salt_merge(i,j) = nanmean(s(dis<=0.25));
    end
end
t = temp_merge;
s = salt_merge;
end
function z = Cressman_bkg(xo,yo,zo,mask,depth)
z0 = nan(1440,121);
x = xo(~isnan(zo));
y = yo(~isnan(zo));
zz = zo(~isnan(zo));

parfor i = 1:1440
    for j = 1:121
        r2 = power(distance(y,x,yo(i,j),xo(i,j),6371),2);
%         Lx = 5*120000/(0.35*yo(i,j)*yo(i,j)+300);
%         Ly = 5*120000/(0.35*yo(i,j)*yo(i,j)+400);
%         R2 = Lx*Lx+Ly*Ly;
        R2 = power(321,2);
%         R2 = power(999,2);
        tt = zz(r2<R2);
        if (mask(i,j)~=-1)&&mask(i,j)>=depth
            if ~isempty(tt)
                ww = (R2-r2(r2<R2))./(R2+r2(r2<R2));
                if numel(tt)==1
                    z0(i,j) = sum(ww(~isnan(tt)).*tt(~isnan(tt)))/(sum(ww(~isnan(tt)))+1);
                elseif numel(tt)>1
                    z0(i,j) = sum(ww(~isnan(tt)).*tt(~isnan(tt)))/sum(ww(~isnan(tt)));
                end
            end
        end
    end
end
idx_land = mask==-1;
z1 = z0;
z1(idx_land) = nan;

parfor i = 1:1440
    for j = 1:121
        if ~isnan(zo(i,j))
        r2 = power(distance(y,x,yo(i,j),xo(i,j),6371),2);
        R2 = power(267,2);
%         R2 = power(666,2);
%         Lx = 3*120000/(0.35*yo(i,j)*yo(i,j)+300);
%         Ly = 3*120000/(0.35*yo(i,j)*yo(i,j)+400);
%         R2 = Lx*Lx+Ly*Ly;
        tt = zz(r2<R2);
        if (mask(i,j)~=-1)&&mask(i,j)>=depth
            if ~all(isnan(tt))
                ww = (R2-r2(r2<R2))./(R2+r2(r2<R2));
                if numel(tt)==1
                    z1(i,j) = z1(i,j) + sum(ww.*(tt-z1(i,j)))/(sum(ww(~isnan(tt)))+1);
                elseif numel(tt)>1
                    z1(i,j) = z1(i,j) + sum(ww.*(tt-z1(i,j)))/sum(ww(~isnan(tt)));
                end
            end
        end
    end
    end
end

z2 = z1;
z2(idx_land) = nan;


parfor i = 1:1440
    for j = 1:121
        if ~isnan(zo(i,j))
        r2 = power(distance(y,x,yo(i,j),xo(i,j),6371),2);
%         Lx = 2*120000/(0.35*yo(i,j)*yo(i,j)+300);
%         Ly = 2*120000/(0.35*yo(i,j)*yo(i,j)+400);
%         R2 = Lx*Lx+Ly*Ly;
        R2 = power(214,2);
% R2 = power(333,2);
        tt = zz(r2<R2);
        if (mask(i,j)~=-1)&&mask(i,j)>=depth
            if ~all(isnan(tt))
                ww = (R2-r2(r2<R2))./(R2+r2(r2<R2));
                if numel(tt)==1
                    z2(i,j) = z2(i,j) + sum(ww.*(tt-z2(i,j)))/(sum(ww(~isnan(tt)))+1);
                elseif numel(tt)>1
                    z2(i,j) = z2(i,j) + sum(ww.*(tt-z2(i,j)))/sum(ww(~isnan(tt)));
                end
            end
        end
        end
    end
end
z = z2;
z(idx_land) = nan;


end
function z3 = Barnes_sc(xo,yo,z,zo,k,mask,depth)
x = xo(~isnan(z));
y = yo(~isnan(z));
zz = z(~isnan(z));
z3 = successive_correction(x,y,xo,yo,zz,zo,k,mask,depth);
%%
    function z2 = successive_correction(x,y,xo,yo,zz,zo,k,mask,depth)
        m = 1;
        for i = 1:k
            z1 = single_correction(x,y,xo,yo,zz,zo,m,mask,depth);
            m = m +1;
            zo = z1;
        end
        z2 = zo;
    end
%%
    function z1 = single_correction(x,y,xo,yo,zz,zo,m,mask,depth)
        if m == 1
            R = 999;
            n = 8e4;
        elseif m == 2
            R = 999;
            n = 16e3;
        end
        parfor i = 1:1440
            for j = 1:121
                r2 = power(distance(y,x,yo(i,j),xo(i,j),6371),2);
                %                 Lx =  120000/(0.35*yo(i,j)*yo(i,j)+300);
                %                 Ly = 120000/(0.35*yo(i,j)*yo(i,j)+400);
                %                 R2 = power(abs(Lx+Ly*1i),2);
                R2 = power(R,2);
                tt = zz(r2<R2);
                if (mask(i,j)~=-1)&&mask(i,j)>=depth
                    if ~all(isnan(tt))
                        w = exp(-1.*r2(r2<R2)./(0.2*n));
                        if numel(tt)==1
                            z1(i,j) = zo(i,j)+sum(w.*(tt-zo(i,j)))/(sum(w)+1);
                        elseif numel(tt)>1
                            z1(i,j) = zo(i,j)+sum(w.*(tt-zo(i,j)))/sum(w);
                        else
                            z1(i,j) = zo(i,j);
                        end
                    else
                        z1(i,j) = zo(i,j);
                    end
                else
                    z1(i,j) = nan;
                end
            end
        end
    end
end









function z = Barnes(x,y,zo,xo,yo,k,filename)
z0 = nan(1440,121);
parfor i = 1:1440
    for j = 1:121
        
        r2 = power(distance(y,x,yo(i,j),xo(i,j),6371),2);
        Lx =  max(2*120000/(0.35*yo(i,j)*yo(i,j)+300),1000);
        Ly = max(2*120000/(0.35*yo(i,j)*yo(i,j)+400),1000);
        R2 = power(abs(Lx+Ly*1i),2);
        % R2 = power(555,2);
        tt = zo(r2<R2);
        if ~isempty(tt)
            %             w = exp(-4.*r2(r2<R2)./R2);
            ww = exp((R2-r2(r2<R2))./(R2+r2(r2<R2)));
            if numel(tt)==1
                z0(i,j) = sum(ww.*tt)/(sum(ww)+1);
            elseif numel(tt)>1
                z0(i,j) = sum(ww.*tt)/sum(ww);
            end
        end
    end
end
parfor i = 1:1440
    for j = 1:121
        r2 = power(distance(y,x,yo(i,j),xo(i,j),6371),2);
        Lx =  120000/(0.35*yo(i,j)*yo(i,j)+300);
        Ly = 120000/(0.35*yo(i,j)*yo(i,j)+400);
        R2 = power(abs(Lx+Ly*1i),2);
        tt = zo(r2<R2);
        if ~isempty(tt)
            w = exp(-4.*r2(r2<R2)./(0.2*8e4));
            if numel(tt)==1
                z(i,j) = z0(i,j)+sum(w.*(tt-z0(i,j)))/(sum(w)+1);
            elseif numel(tt)>=10
                z(i,j) = z0(i,j)+sum(w.*(tt-z0(i,j)))/sum(w);
            else
                z(i,j) = z0(i,j);
            end
        else
            w=0;
            z(i,j) = z0(i,j);
        end
    end
end
end