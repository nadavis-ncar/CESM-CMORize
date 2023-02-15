%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Scan dictionary and return field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=translate_cesm(dictionary,in,field,closest_match)

out=[];

dictionary_entries=fieldnames(eval(['dictionary.',field]));
for i=1:length(dictionary_entries)
   if strcmp(dictionary_entries{i},in)
      out=eval(['dictionary.',field,'.',dictionary_entries{i}]);
   end
end

if closest_match==1
   if isempty(out)
      for i=1:length(dictionary_entries)
         if contains(in,dictionary_entries{i})
            out=eval(['dictionary.',field,'.',dictionary_entries{i}]);
         end
      end
   end
end

if isempty(out)
   out='no_match';
end
