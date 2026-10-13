clearvars -except arctic_data arctic_data_climatology_qc ...
    observ_mon_final dataset_mon observ_mon_final arc_season arctic_data_1 ori arc_mon ori data_review AOTD
MONTH = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', ...
    'Sep', 'Oct', 'Nov', 'Dec'};
SEASON = {'Spr','Sum','Aut','Win'};
%% ʱ ռ Ȩƽ      ̬
lon = -180:0.25:179.75;
lat = 60:0.25:90;
[LAT_AIM,LON_AIM] = meshgrid(lat,lon);
Lt = 1000;
load('/data2/ljl/master_arctic/topographic_land.mat')
%a = dataset_mon.May;
%%depth = arc_mon.depth;
mkdir('/data2/ljl//master_arctic//barnes_test9');
% load('topographic_land.mat');
% idx_land = isnan(MASK);
fi = 0;
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
depth = AOTD.pres;
mon = [AOTD.mo];
parpool(32)
for j = 1:12
    idx = (mon==j);
    arctic_data_1 = AOTD(idx);
    lon_obs = [arctic_data_1.lon];
    lat_obs = [arctic_data_1.lat];
    temp = [arctic_data_1.temp];
    salt = [arctic_data_1.salt];
%     month_name = sprintf('/data2//ljl//master_arctic//barnes_test7//climatology_seasonal_%02d.nc',j/3);
    file_name = sprintf('/data2//ljl//master_arctic//barnes_test9//climatology_%02d.nc',j);
    if j ==1||j ==2||j ==3
        season_name = sprintf('/data2/ljl/master_arctic/barnes_test9/climatology_seasonal_%02d.nc',1);
    elseif j ==4||j ==5||j ==6
        season_name = sprintf('/data2/ljl/master_arctic/barnes_test9/climatology_seasonal_%02d.nc',2);
    elseif j ==7||j ==8||j ==9
        season_name = sprintf('/data2/ljl/master_arctic/barnes_test9/climatology_seasonal_%02d.nc',3);
    else
        season_name = sprintf('/data2/ljl/master_arctic/barnes_test9/climatology_seasonal_%02d.nc',4);
    end
    nccreate(file_name, 'geolon', 'Dimensions', {'lon', 1440,'lat',121});
    nccreate(file_name, 'geolat', 'Dimensions', {'lon', 1440,'lat',121});
    nccreate(file_name, 'temp', ...
        'Dimensions', {'lon', 1440, 'lat', 121, 'depth', 57});
    nccreate(file_name, 'salt', ...
        'Dimensions', {'lon', 1440, 'lat', 121, 'depth', 57});
    nccreate(file_name, 'depth', ...
        'Dimensions', {'depth', 57});
    for i = 1:57
        tic
        t = temp(i,:);
        %                 if ~all(isnan(t))
        s = salt(i,:);
        idx_ii = find(~isnan(t));
        arctic_data_11 = arctic_data_1(idx_ii);
        lonn = [arctic_data_11.lon];
        latt = [arctic_data_11.lat];
%         time = datetime({arctic_data_11.time})';
        %     [t,s] = obs_modas(latt,lonn,t(~isnan(t)),s(~isnan(t)),LAT_AIM,LON_AIM,time);
        %     idx_land = MASK==-1;
        %     t(idx_land) = nan;
        %     s(idx_land) = nan;
        [t,s] = obs_merge(latt,lonn,t(~isnan(t)),s(~isnan(t)),LAT_AIM,LON_AIM);
        idx_land = MASK==-1;
        t(idx_land) = nan;
        s(idx_land) = nan;
        h = fspecial('average', [3 3]);
        zo = ncread(season_name,'temp',[1,1,i],[inf,inf,1]);
%         z1 = Cressman_bkg(LON_AIM,LAT_AIM,zo,t,MASK,depth(i));
z1 = Barnes_sc(LON_AIM,LAT_AIM,t,zo,2,MASK,depth(i));
        idx = isnan(z1);
%         z1 = fillmissing(z1,'previous',1);
%         z1 = fillmissing(z1,'previous',2);
%         z1 = fillmissing(z1,'movmedian',9,1);
%         z1 = fillmissing(z1,'movmedian',9,2);
%         z1 = imfilter(z1, h, 'replicate');
%         z1(idx)=nan;
                
        
        zo = ncread(season_name,'salt',[1,1,i],[inf,inf,1]);
%         z2 = Cressman_bkg(LON_AIM,LAT_AIM,zo,s,MASK,depth(i));
                z2 = Barnes_sc(LON_AIM,LAT_AIM,s,zo,2,MASK,depth(i));
%         z2 = fillmissing(z2,'previous',1);
%         z2 = fillmissing(z2,'previous',2);
%         z2 = fillmissing(z2,'movmedian',9,1);
%         z2 = fillmissing(z2,'movmedian',9,2);
        %         h = fspecial('average', [11 11]);
%         z2 = imfilter(z2, h, 'replicate');
%         z2(idx)=nan;
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
function [t,s,num] = obs_merge(latt,lonn,t,s,LAT_AIM,LON_AIM)
temp_merge = nan(1440,121);
salt_merge = nan(1440,121);
num = nan(1440,121);
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
        num(i,j) = numel(t(dis<=0.25));
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
zzmask = z;
zz = z(~isnan(z));
z3 = successive_correction(x,y,xo,yo,zz,zo,k,mask,depth,zzmask);
%%
    function z2 = successive_correction(x,y,xo,yo,zz,zo,k,mask,depth,zzmask)
        m = 1;
        for i = 1:k
            z1 = single_correction(x,y,xo,yo,zz,zo,m,mask,depth,zzmask);
            m = m +1;
            zo = z1;
        end
        z2 = zo;
    end
%%
    function z1 = single_correction(x,y,xo,yo,zz,zo,m,mask,depth,zzmask)
        if m == 1
            R = 555;
            n = 8e4;
        elseif m == 2
            R = 555;
            n = 16e3;
        end
        parfor i = 1:1440
            for j = 1:121
%                 if ~isnan(zzmask(i,j))
                    r2 = power(distance(y,x,yo(i,j),xo(i,j),6371),2);
                    %                 Lx =  120000/(0.35*yo(i,j)*yo(i,j)+300);
                    %                 Ly = 120000/(0.35*yo(i,j)*yo(i,j)+400);
                    %                 R2 = power(abs(Lx+Ly*1i),2);
                    R2 = power(R,2);
                    tt = zz(r2<R2);
                    if ~isnan(mask(i,j))&&mask(i,j)>=depth
                        if ~all(isnan(tt))
                            w = exp(-1.*r2(r2<R2)./(0.2*n));
                            if numel(tt)==1
                                z1(i,j) = zo(i,j)+sum(w.*(tt-zo(i,j)))/(sum(w)+1);
                            elseif numel(tt)>1
                                z1(i,j) = zo(i,j)+sum(w.*(tt-zo(i,j)))/sum(w);
                            else
%                                 z1(i,j) = zo(i,j)+sum(w.*(tt-zo(i,j)))/sum(w);
                                z1(i,j) = zo(i,j);
                            end
                        else
                            z1(i,j) = zo(i,j);
                        end
                    else
                        z1(i,j) = nan;
                    end
%                 else
%                     z1(i,j) = zo(i,j);
%                 end
            end
        end
    end
end