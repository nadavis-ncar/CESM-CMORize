%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get index variables for model run (r#i#f#p#)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function globals=set_index_variables(globals,case_name)

index_vars={'realization_index';'initialization_index';'physics_index';'forcing_index';'variant_label'};
case_array=convertStringsToChars(case_name);

for i=1:length(index_vars)
   for j=1:length(globals)
      if strcmp(globals(j).name,index_vars{i})
         if ~strcmp(index_vars{i},'variant_label')
            globals(j).value=case_array(i*2);
         else
            globals(j).value=case_name;
         end
      end
   end
end
