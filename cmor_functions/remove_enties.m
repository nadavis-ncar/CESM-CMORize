%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Remove entries from a list
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function list=remove_entries(list,entries)

delete_list=[];
for i=1:length(entries)
   for j=1:length(list)
      if strcmp([list(j).folder,'/',list(j).name],entries{i})
         delete_list=cat(1,delete_list,j);
      end
   end
end
list(delete_list)=[];
