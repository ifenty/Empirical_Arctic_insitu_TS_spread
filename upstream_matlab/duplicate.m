% parpool(8);
fi = 0;
for ii = 1:7
    if ii==1
        filepath = '../ITP7';
        source = 1;
    elseif ii == 2
        filepath = '../Argo7';
        source = 2;
    elseif ii ==3
        filepath = '../BOL7';
        source = 3;
    elseif ii == 4
        filepath = '../Nabos7';
        source = 4;
    elseif ii == 5
        filepath = '../CNarc7';
        source = 5;
    elseif ii == 6
        filepath = '../WOD7';
        source = 6; 
    elseif ii == 7
        filepath = '../CMEMS7';
        source = 7;
    end
    % filepath = '../AOTD_backup/WOD';
    filelist = dir(fullfile(filepath,'**','*.nc'));
    % outpath = '../WOD_tttest';
    % parpool(8);
    for i = 1:length(filelist)
        tic
        filename = fullfile(filelist(i).folder,filelist(i).name);
        
        lon0 = ncread(filename,'lon');
        lat0 = ncread(filename,'lat');
        time0 = ncread(filename,'time');
        temp0 = ncread(filename,'temp');
        salt0 = ncread(filename,'salt');
        pres0 = ncread(filename,'pres');
        dens0 = ncread(filename,'dens');
        lat_qc = ncread(filename,'lat_qc');
        lon_qc = ncread(filename,'lon_qc');
        lat0(lat_qc~=1|lon_qc~=1)=nan;
        lon0(lat_qc~=1|lon_qc~=1)=nan;
        lat_qc = repmat(lat_qc,1,size(temp0,1))';
        lon_qc = repmat(lon_qc,1,size(temp0,1))';
        temp_qc = ncread(filename,'temp_qc');
        salt_qc = ncread(filename,'salt_qc');
        pres_qc = ncread(filename,'pres_qc');
        dens_qc = ncread(filename,'dens_qc');
        time_qc = ncread(filename,'time_qc');
        temp0(temp_qc~=2|pres_qc~=1|lat_qc~=1|lon_qc~=1)=nan;
        salt0(salt_qc~=2|pres_qc~=1|lat_qc~=1|lon_qc~=1)=nan;
        pres0(pres_qc~=1|lat_qc~=1|lon_qc~=1)=nan;
        [~,name,~] = fileparts(filename);
        ts = strsplit(name,'_');
        ts = strcat(ts{2},'01');
        time = datetime(time0+datenum(datetime(ts,'InputFormat','yyyyMMdd')),'ConvertFrom','datenum');
        ye = year(time);
        mo = month(time);
        da = day(time);
        ho = hour(time);
        mi = minute(time);
        se = second(time);
        if ii == 1
            for j = 1:length(lat0)
                fi = fi +1;
                ori(fi).source = source;
                ori(fi).lat = lat0(j);
                ori(fi).lon = lon0(j);
                ori(fi).ye = ye(j);
                ori(fi).mo = mo(j);
                ori(fi).da = da(j);
                ori(fi).ho = ho(j);
                ori(fi).mi = mi(j);
                ori(fi).se = se(j);
                ori(fi).temp = temp0(:,j);
                ori(fi).salt = salt0(:,j);
                ori(fi).pres = pres0(:,j);
                ori(fi).dens = dens0(:,j);
                
            end
        else
            latt = [ori.lat];
            lonn = [ori.lon];
            yee = [ori.ye];
            moo = [ori.mo];
            daa = [ori.da];
            hoo = [ori.ho];
            mii = [ori.mi];
            see = [ori.se];
            timm = datetime(yee, moo, daa, hoo, mii, see);
            for k = 1:length(lat0)
                dt = hours(abs(time(k)-timm));
                dlat = abs(lat0(k) - latt);
                dlon = abs(lon0(k) - lonn);
                ind = dt<6&dlat<0.1&dlon<0.1;
                idx = find(ind==1);
                if isempty(idx)
                    fi = fi +1;
                    ori(fi).source = source;
                    ori(fi).lat = lat0(k);
                    ori(fi).lon = lon0(k);
                    ori(fi).ye = ye(k);
                    ori(fi).mo = mo(k);
                    ori(fi).da = da(k);
                    ori(fi).ho = ho(k);
                    ori(fi).mi = mi(k);
                    ori(fi).se = se(k);
                    ori(fi).temp = temp0(:,k);
                    ori(fi).salt = salt0(:,k);
                    ori(fi).pres = pres0(:,k);
                    ori(fi).dens = dens0(:,k);
                end

            end
        end
    end
    clearvars -except ori fi ii
end

