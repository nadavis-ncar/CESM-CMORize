%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Determines the index number of the given index name
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function index=get_dim_index(dim,dim_name,closeness)

index=[];
for i=1:length(dim)
   if strcmp(closeness,'exact')
      if strcmp(dim(i).name,dim_name)
         index=i;
      end
   elseif strcmp(closeness,'contains')
      if contains(dim(i).name,dim_name)
         index=i;
      end
   else
      error('unrecognized closeness parameter - must be "exact" or "contains"')
   end
end

if isempty(index)
   error('index not found')
end
