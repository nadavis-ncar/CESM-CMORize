%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get hybrid coefficient formula attributes 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function info=get_formula_variable(formula_terms,var)

names=fieldnames(formula_terms);
formula_terms=struct2cell(formula_terms);
for i=1:length(formula_terms)
   if strcmp(names{i},var)
      info=formula_terms(i);
   end
end
