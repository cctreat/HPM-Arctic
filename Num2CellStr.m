function [outcell] = Num2CellStr(vect)
if size(vect,1) == 1
    vect=vect';
end

outcell = cell(length(vect),1);
for i=1:length(vect)
    outcell{i} = num2str(vect(i),3);
end