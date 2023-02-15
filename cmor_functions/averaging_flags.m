%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Determine which dims need to be averaged
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function averaging=averaging_flags(cell_methods)

parsed=parse_string(cell_methods);
index=1;
ave_count=1;

while ~isempty(parsed)
   %Dimension
   if contains(parsed{index},':')
      if length(parsed)>index
         %Do no averaging, remove entry 
         if contains(parsed{index+1},':')
            parsed(index)=[];
         else
            %Add averaging flag, remove entry pair
            averaging{ave_count}=strrep(parsed{index},':','');
            ave_count=ave_count+1;
            parsed(index:index+1)=[];
         end
      else
         %Remove entry
         parsed(index)=[];
      end
   else
      %Move along (probably unnecessary?)
      index=index+1;
   end
end

