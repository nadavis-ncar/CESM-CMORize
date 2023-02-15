%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate sum along arbitrary dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data=calculate_sum(data,dim_name,varargin)

dim_num=data.dim_to_remove(find(strcmp(data.dim_names_to_remove,dim_name),1,'first'));
data.native.value=squeeze(sum(data.native.value,dim_num));

