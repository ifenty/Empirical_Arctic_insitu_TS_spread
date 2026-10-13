% load('G:\���ݼ�\����\��ֵ�����ݼ�\arctic_data.mat')
% load('/data2/ljl/master_arctic/arctic_data_mon.mat')
clearvars -except arctic_data arctic_data_climatology_qc observ_mon_final dataset_mon ori data_review AOTD
MONTH = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', ...
    'Sep', 'Oct', 'Nov', 'Dec'};
%% ʱ�ռ�Ȩƽ������̬
lon = -180:0.25:179.75;
lat = 60:0.25:90;
[LAT_AIM,LON_AIM] = meshgrid(lat,lon);
Lt = 1000;
load('/data2/ljl/master_arctic/topographic_land.mat')
load('std_data.mat')
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
% arctic_data_1 = ori;
climatology_temp = nan(1440,121,57);
climatology_salt = nan(1440,121,57);
% %     idx_TT = find(mon == TT);
%     arctic_data_1 = observ_mon_final.(MONTH{TT});
% time = {ori.time}';
% mon = nan(length(time),1);
% for i = 1:length(time)
%    timm = time{i,1};
%    mon(i,1) = str2double(timm(6:7)); 
% end
mon = [AOTD.mo];
fi = fi+1;
parpool(52)
for j = 1:4
    if j == 1
    idx = mon==3|mon==2|mon==1;
    elseif j == 2
        idx = mon==4|mon==5|mon==6;
    elseif j == 3
        idx = mon==7|mon==8|mon==9;
    elseif j == 4
        idx = mon==10|mon==11|mon==12;
    end
    arctic_data_1 = AOTD(idx);
    file_name = sprintf('/data2//ljl//master_arctic//barnes_test9//climatology_seasonal_%02d.nc',j);
    nccreate(file_name, 'geolon', 'Dimensions', {'lon', 1440,'lat',121});
    nccreate(file_name, 'geolat', 'Dimensions', {'lon', 1440,'lat',121});
    nccreate(file_name, 'temp', ...
        'Dimensions', {'lon', 1440, 'lat', 121, 'depth', 57});
    nccreate(file_name, 'salt', ...
        'Dimensions', {'lon', 1440, 'lat', 121, 'depth', 57});
    nccreate(file_name, 'depth', ...
        'Dimensions', {'depth', 57});

lon_obs = [arctic_data_1.lon];
lat_obs = [arctic_data_1.lat];
temp = [arctic_data_1.temp];
salt = [arctic_data_1.salt];
for i = 1:57
    tic
    t = temp(i,:);
    %                 if ~all(isnan(t))
    s = salt(i,:);
    idx_ii = find(~isnan(t));
    arctic_data_11 = arctic_data_1(idx_ii);
    lonn = [arctic_data_11.lon];
    latt = [arctic_data_11.lat];
%     time = datetime({arctic_data_11.time})';
%     [t,s] = obs_modas(latt,lonn,t(~isnan(t)),s(~isnan(t)),LAT_AIM,LON_AIM);
%     idx_land = MASK==-1;
%     t(idx_land) = nan;
%     s(idx_land) = nan;
    [t,s] = obs_merge(latt,lonn,t(~isnan(t)),s(~isnan(t)),LAT_AIM,LON_AIM);
    idx_land = MASK==-1;
    t(idx_land) = nan;
    s(idx_land) = nan;
    zo = ncread('/data2/ljl/master_arctic/barnes_test9/climatology_annual.nc','temp',[1,1,i],[inf,inf,1]);
    z1 = Cressman_bkg(LON_AIM,LAT_AIM,zo,t,MASK,depth(i));
%     h = fspecial('average', [3 3]);
%     z1 = fillmissing(z1,'previous',1);
%     z1 = fillmissing(z1,'previous',2);
%     z1 = imfilter(z1, h, 'replicate');
%     z1 = fillmissing(z1,'movmedian',9,1);
%     z1 = fillmissing(z1,'movmedian',9,2);
    z1(isnan(zo))=nan;
%     z1 = Barnes_sc(LON_AIM,LAT_AIM,t,zo,2,MASK,depth(i));
    zo = ncread('/data2/ljl/master_arctic/barnes_test9/climatology_annual.nc','salt',[1,1,i],[inf,inf,1]);
    z2 = Cressman_bkg(LON_AIM,LAT_AIM,zo,s,MASK,depth(i));
%     z2 = fillmissing(z2,'linear',2);
%     z2 = fillmissing(z2,'linear',1);
%     z2 = imfilter(z2, h,'symmetric');
%     z2 = fillmissing(z2,'movmedian',9,1);
%     z2 = fillmissing(z2,'movmedian',9,2);
    z2(isnan(zo))=nan;
    climatology_temp(:,:,i) = z1;
    climatology_salt(:,:,i) = z2;
    toc
end
ncwrite(file_name,'geolon',LON_AIM);
ncwrite(file_name,'geolat',LAT_AIM);
ncwrite(file_name, 'temp', climatology_temp);
ncwrite(file_name, 'salt', climatology_salt);
ncwrite(file_name, 'depth', depth(1:57));

end


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
        d = distance(lat,lon,latt,lonn,6371);     
        for k = 1:2
            Lx = 2*k*120000/(0.35*lat*lat+300);
            Ly = 2*k*120000/(0.35*lat*lat+400);
            Lt = 1000;
            L = sqrt(Lx*Lx+Ly*Ly);
            idx = d<=L;
%             if sum(idx)<100&&k==1
%                 continue
%             else
                t_idx = t(idx);
                s_idx = s(idx);
%                 time_idx = time(idx)';
                lat_idx = latt(idx);
                lon_idx = lonn(idx);
%                 timmm = timm(idx)';
                b = exp(-power((lat-lat_idx)/Lx,2)-power((lon-lon_idx)/Ly,2));
                temp_merge(i,j) = sum(b.*t_idx)/sum(b);
                salt_merge(i,j) = sum(b.*s_idx)/sum(b);
%             end
            t_idx = t(idx);
            s_idx = s(idx);
%             time_idx = time(idx)';
            lat_idx = latt(idx);
            lon_idx = lonn(idx);
%             timmm = timm(idx)';
            b = exp(-power((lat-lat_idx)/Lx,2)-power((lon-lon_idx)/Ly,2));
            temp_merge(i,j) = sum(b.*t_idx)/sum(b);
            salt_merge(i,j) = sum(b.*s_idx)/sum(b);
        end

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
        temp_merge(i,j) = mean(t(dis<=0.25));
        salt_merge(i,j) = mean(s(dis<=0.25));
    end
end
t = temp_merge;
s = salt_merge;
end



function z = Cressman_bkg(xo,yo,zzz,zo,mask,depth)
x = xo(~isnan(zo));
y = yo(~isnan(zo));
zz = zo(~isnan(zo));
z0 = zzz;
idx_land = mask==-1;
z1 = z0;
z1(idx_land) = nan;

parfor i = 1:1440
    for j = 1:121
        if ~isnan(zo(i,j))
            r2 = power(distance(y,x,yo(i,j),xo(i,j),6371),2);
            R2 = power(321,2);
            tt = zz(r2<R2);
            %         ttt = zz(r2<R3);
            if (mask(i,j)~=-1)&&mask(i,j)>=depth
                if ~all(isnan(tt))
                    ww = (R2-r2(r2<R2))./(R2+r2(r2<R2));
                    if numel(tt)==1
                        z1(i,j) = z1(i,j) + sum(ww.*(tt-z1(i,j)))/(sum(ww(~isnan(tt)))+1);
                        %                     continue
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
            R2 = power(267,2);
            %         Lx = min(1*120000/(0.35*yo(i,j)*yo(i,j)+300),111);
            %         Ly = min(1*120000/(0.35*yo(i,j)*yo(i,j)+400),111);
            %         R2 = Lx*Lx+Ly*Ly;
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
z3 = z2;
z3(idx_land) = nan;
parfor i = 1:1440
    for j = 1:121
        if ~isnan(zo(i,j))
            r2 = power(distance(y,x,yo(i,j),xo(i,j),6371),2);
            %         Lx = 2*120000/(0.35*yo(i,j)*yo(i,j)+300);
            %         Ly = 2*120000/(0.35*yo(i,j)*yo(i,j)+400);
            %         R2 = Lx*Lx+Ly*Ly;
            R2 = power(214,2);
            % R2 = power(55,2);
            tt = zz(r2<R2);
            if (mask(i,j)~=-1)&&mask(i,j)>=depth
                if ~all(isnan(tt))
                    ww = (R2-r2(r2<R2))./(R2+r2(r2<R2));
                    if numel(tt)==1
                        z3(i,j) = z3(i,j) + sum(ww.*(tt-z3(i,j)))/(sum(ww(~isnan(tt)))+1);
                    elseif numel(tt)>1
                        z3(i,j) = z3(i,j) + sum(ww.*(tt-z3(i,j)))/sum(ww(~isnan(tt)));
                    end
                end
            end
        end
    end
end

z = z3;
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