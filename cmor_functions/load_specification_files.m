%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load specification files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [specification,specification_files]=load_specification_files(file_preamble)

specification=struct;
specs={'CV';'fx';'grids';'coordinate';'formula_terms'};
for i=1:length(specs)
   specification_files{i}=[file_preamble,specs{i},'.json'];
   eval(['specification.',specs{i},'=json_load(''',specification_files{i},''');']);
end
