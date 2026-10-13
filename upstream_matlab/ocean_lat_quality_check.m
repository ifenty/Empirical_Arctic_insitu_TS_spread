function lat_qc = ocean_lat_quality_check(lat,lat_range)
% lat_rangeÊÇlatµÄÆðÊ¼·¶Î§
lat_qc = lat>lat_range(1,1)&lat<lat_range(1,2);
end