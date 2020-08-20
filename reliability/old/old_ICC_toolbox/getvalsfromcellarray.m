function outvec=getvalsfromcellarray(celldata,x,y,z)
% length(outvec)=length(celldata)

outvec=cellfun(@(v) v(x,y,z), celldata);

if iscell(outvec)
    outvec=cell2mat(outvec); % esp for ICC data and stats
end