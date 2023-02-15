%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate age of air
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tracer=age_of_air(tracer,p_ref,y_ref)

aoa=zeros(size(tracer.native.value));

%Setup indexing for eval
dim_string='(';
for i=1:length(size(tracer.native.value))
   dim_name=tracer.native_dim(i).name;
   switch dim_name
      case 'time'
         string_add=':';
         time_dim=i;
      case 'lev'
         string_add='p';
         lev_dim=i;
      case 'lat'
         string_add='y';
         lat_dim=i;
      case 'lon'
         string_add='x';
         lon_dim=i;
      otherwise
         error('Unrecognized dimension name, cannot compute age of air')
   end
   if i<length(size(tracer.native.value))
      string_add=[string_add,','];
   else
      string_add=[string_add,')'];
   end
   dim_string=[dim_string,string_add];
end
time_string=strrep(dim_string,':','t');
ref_string=strrep(dim_string,'x',':');
ref_string=strrep(ref_string,'y','y_ref');
ref_string=strrep(ref_string,'p','p_ref');

%Get reference value
[val,y_ref]=min(abs(tracer.native_dim(lat_dim).value-y_ref));
[val,p_ref]=min(abs(tracer.native_dim(lev_dim).value-p_ref));
disp(['Reference latitude ',sprintf('%2.1f',tracer.native_dim(lat_dim).value(y_ref))])
disp(['Reference pressure ',sprintf('%2.1f',tracer.native_dim(lev_dim).value(p_ref))])
tracer.native.value=movmean(tracer.native.value,12,time_dim);
reference_timeseries=eval(['squeeze(mean(tracer.native.value',ref_string,',',num2str(lon_dim),'));']);
reference_timeseries=repmat(reference_timeseries(:),[1 length(reference_timeseries)]);

nt=size(tracer.native.value,time_dim);

%Calculate aoa using the matrix method
for y=1:length(tracer.native_dim(lat_dim).value)
   disp([sprintf('%2.2f',100*(y-1)/length(tracer.native_dim(lat_dim).value)),'% complete'])
   for x=1:length(tracer.native_dim(lon_dim).value)
      for p=1:length(tracer.native_dim(lev_dim).value)
         timeseries=eval(['squeeze(tracer.native.value',dim_string,')']);
         timeseries=repmat(timeseries(:),[1 length(timeseries)])';
         [val,time_index]=min(abs(timeseries-reference_timeseries),[],2);
         eval(['aoa',dim_string,'=(time_index''-[1:nt])/12;']);
      end
   end
end

tracer.native.value=aoa;
