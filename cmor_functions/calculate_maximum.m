%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate maximum value along arbitrary dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data=calculate_maximum(data,dim_name)

dim_num=data.dim_to_remove(find(strcmp(data.dim_names_to_remove,dim_name),1,'first'));
data.native.value=squeeze(max(data.native.value,[],dim_num));

