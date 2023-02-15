%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Translate variable details, make actionable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [variable_out]=translate_cesm_variable(in)

opt_count=1;
in_components=fieldnames(in);
variables_components={'var';'operation';'eval';'p0';'lat0';'axis';'option'};
variable_out=struct('var',strings(0),'operation',strings(0),'eval',strings(0),'p0',[],'lat0',[],'axis',strings(0));

for i=1:length(in_components)
   for j=1:length(variables_components)-1
      if contains(in_components{i}, variables_components{j})
         eval(['variable_out.',variables_components{j},'=cat(1,variable_out.',variables_components{j},',in.',in_components{i},');'])
      end
   end
   if contains(in_components{i}, 'option')
      variable_out.option(opt_count)=struct('var',eval(['in.',in_components{i},'.var']),'string',eval(['in.',in_components{i},'.string']));
      opt_count=opt_count+1;
   end
end
