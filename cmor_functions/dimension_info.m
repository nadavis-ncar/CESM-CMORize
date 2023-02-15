%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Returns dimension information from specification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function info=dimension_info(specification,dimension,dictionary)

info=strings(0);
names=fieldnames(specification.coordinate.axis_entry);
for i=1:length(names)
   if strcmp(names{i},dimension)
      info=eval(['specification.coordinate.axis_entry.',names{i}]);
   end
end

if isempty(info)
   info=translate_cesm(dictionary,dimension,'Dimension',0);
end
