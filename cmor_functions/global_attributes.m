%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Gather required global attributes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function globals=global_attributes(required_attributes,CV_file,output,cmor_specification,specification,frequency,var,cesm_globals)

globals=struct;
cesm_globals_name=fieldnames(cesm_globals);
cesm_globals=struct2cell(cesm_globals);
mip_globals_name=fieldnames(CV_file);
mip_globals=struct2cell(CV_file);

for i=1:length(required_attributes)
   globals(i).name=required_attributes{i};
   overwrite=0;

   %Pull values from cesm spec file
   for j=1:length(cesm_globals)
      if strcmp(globals(i).name,cesm_globals_name{j})
         input=cesm_globals{j};
         if isstruct(input)
            input=eval([input.eval]);
         end
         overwrite=1;
      end
   end

   %Otherwise pull from mip file
   if overwrite==0
      for j=1:length(mip_globals)
         if strcmp(globals(i).name,mip_globals_name{j})
            input=mip_globals{j};
            if iscell(input)
               input=input{1};
            end
            if isstruct(input)
               input=struct2cell(input);
               input=input{1};
            end
         end
      end
   end
   globals(i).value=input;
end
