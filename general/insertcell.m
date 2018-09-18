function c = insertcell(c,ins,idx)
% from http://www.mathworks.com/matlabcentral/newsreader/view_thread/323580
    c = [c(1:idx-1) {ins} c(idx:end)];
end