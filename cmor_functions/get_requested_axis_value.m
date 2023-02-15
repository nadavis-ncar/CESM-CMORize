%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get requested axis values or single value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dim=get_requested_axis_value(dim,axis_entries)

axes=struct2cell(axis_entries);
axes_entry=find(strcmp(fieldnames(axis_entries),dim.out.name),1,'first');
if ~isempty(axes{axes_entry}.requested)
   requested=axes{axes_entry}.requested;
   for i=1:length(requested)
      dim.interp(i)=str2double(requested{i});
   end
   dim.interp=dim.interp(:);
elseif ~isempty(axes{axes_entry}.value)
   dim.interp=str2double(axes{axes_entry}.value);
else
   error('No value found for axis requesting specified values')
end

