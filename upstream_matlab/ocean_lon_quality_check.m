function lon_qc = ocean_lon_quality_check(lon,lon_range)
% lon_rangeÊÇlonµÄÆðÊ¼·¶Î§
lon_qc = lon>lon_range(1,1)&lon<lon_range(1,2);
end