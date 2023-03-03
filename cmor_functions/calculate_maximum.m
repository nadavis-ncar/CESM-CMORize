%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate maximum value along arbitrary dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data=calculate_maximum(data,dim_name)

dim_index=get_dim_index(data.native_dim,dim_name,'exact');
data.native.value=squeeze(max(data.native.value,[],dim_index));
data.native_dim(dim_index)=[];
