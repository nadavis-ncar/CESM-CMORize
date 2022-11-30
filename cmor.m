%Script to CMORize CESM output
%Currently only functioning for the atmosphere component
%nadavis@ucar.edu

%Load user-supplied specification details 
cmor_specification_file = 'cmor_specification.json'; 
cmor_specification=json_load(cmor_specification_file);
 
%Load mip-specific details
[specification,specification_files]=load_specification_files([cmor_specification.cmor_specification_dir,cmor_specification.cmor_specification_format]);

%Load output file structures
output_specification_files=dir([cmor_specification.cmor_specification_dir,cmor_specification.cmor_specification_format,'*.json']);
output_specification_files=remove_entries(output_specification_files,specification_files);

%Load CESM->MIP variable dictionary
%Includes information for translating MIP requests for frequency, realm, etc.
cesm_dictionary=json_load(cmor_specification.cmor_variable_dictionary);

%Create output files
for i=4%1:length(output_specification_files)

   output_specification_file=[output_specification_files(i).folder,'/',output_specification_files(i).name];
   output_file=json_load(output_specification_file);

   output=strrep(output_specification_files(i).name,cmor_specification.cmor_specification_format,'');
   output=strrep(output,'.json','');

   vars=fieldnames(output_file.variable_entry);

   ps_init=0;
   interpolant=struct('latitude',[],'longitude',[],'lev',[]);

   %Output file for all variables
   for v=30%1:length(vars)

      %eval converts a string to dynamic code
      %most straightforward way to deal with a decoded json
      local_var_spec=eval(['output_file.variable_entry.',vars{v}]);

      %Create directory
      dir_output=[cmor_specification.cmor_output_dir,cmor_specification.case_name,'/postprocess/output/',output,'/',vars{v}];
      if ~exist(dir_output)
         mkdir(dir_output)
      end

      dir_input_main=[cmor_specification.case_output_dir,cmor_specification.case_name,'/'];

      %Realm, frequency information
      realm=translate_cesm(cesm_dictionary,local_var_spec.modeling_realm,'Realms');
      frequency=translate_cesm(cesm_dictionary,local_var_spec.frequency,'Frequency');
      variables_list=translate_cesm(cesm_dictionary,vars{v},'Variables');
  
      %Model output location
      dir_input=[dir_input_main,realm,'/proc/tseries/',frequency,'_1/',cmor_specification.case_name]; 

      %Perform special operations, scale_factors, or just load the variable 
      if length(fieldnames(variables_list))>1
         [variable]=translate_cesm_variable(variables_list);
      else
         variable{1}=variables_list;
      end
         
      load_ps=0;

      %Load variable(s) 
      for vin=1:length(variable.var)
         var_files=dir([dir_input,'*.',variable.var{vin},'.*']);
         file_name=[var_files(1).folder,'/',var_files(1).name];

         dims=eval(['parse_string(output_file.variable_entry.',vars{v},'.dimensions);']);
         for j=1:length(dims)
            dim{vin}{j}.native.name=translate_cesm(cesm_dictionary,dims{j},'Dimension');    
            dim{vin}{j}.native.value=ncread(file_name,dim{vin}{j}.native.name);
            dim{vin}{j}.out.name=dims{j};
            dim{vin}{j}.out.info=dimension_info(specification,dims{j},cesm_dictionary);
            dim{vin}{j}.interp_special='';

            %is there a requested grid for this dimension?
            if isstruct(dim{vin}{j}.out.info)
               dim{vin}{j}.interp=dim{vin}{j}.out.info.requested;
            end
 
            %do we need to interpolate the vertical grid?
            if strcmp(dim{vin}{j}.native.name,'lev')
               if strcmp(dim{vin}{j}.out.name,'alevel')
                  dim{vin}{j}.out.name='standard_hybrid_sigma';
                  dim{vin}{j}.interp='';
               end
               dim{vin}{j}.interp_special='vertical';
               load_ps=1;
               a=ncread(file_name,'hyam');
               b=ncread(file_name,'hybm');
               dim{vin}{j}.interp='';
            end
         end
         var{vin}.native.value=ncread(file_name,variable.var{vin});
         var{vin}.native.name=variable.var{vin};

         %Only load once for this output format
         if load_ps & ~ps_init
            var_files=dir([dir_input,'*.PS.*']);
            file_name=[var_files(1).folder,'/',var_files(1).name];
            ps=ncread(file_name,'PS');
         end 
      end

      %Inherit dimensions of input var
      var_out.dim=dim{1};

      %Perform arithmetic operation
      if ~isempty(variable.eval) 
         eval(['var_out.native.value=',variable.eval{1},';'])
      else
         var_out=var{1};
      end

      %Special operations - omega to wa, age of air, integrate, max value
         
      %Do any interpolation
      for j=1:length(var_out.dim) 
         if strcmp(var_out.dim{j}.interp,'interpolate')
            if strcmp(var_out.dim{j}.interp_special,'vertical')
               [var_out,interpolant.lev]=interpolate_field(var_out.value,j,var_out.dim{j},a,b,ps,interpolant.lev);
            else
               eval(['[var_out,interpolant.',dims{j},']=interpolate_field(var_out.value,j,var_out.dim{j},interpolant.',dims{j},');']);
            end
         end
      end

      %Gather global attributes
   %globals=global_attributes(specification.CV.CV.required_global_attributes,specification.CV.CV,output_file);   

      %Gather variable attributes
      %bnd variables

      %Save to file
           
   end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sets up the interpolant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function interpolant=initialize_interpolant(dim_in,dim_out)

%Scan each grid point and determine nearest input points for given output


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Interpolates along specified dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [field_out,interpolant]=interpolate_field(field,j,dim,varargin)

a=[];
b=[];
ps=[];
interpolant=[];
if nargin>3
   a=varargin{1};
   b=varargin{2};
   ps=varargin{3};
end
if nargin==7
   interpolant=varargin{4};
end

dimsizes=size(field);
dims=1:length(dimsizes);

%Reshape field to have j on right, arbitrary 5D array
dims_reshaped=dims;
dims_reshaped(j)=[];
field=permute(field,[dims_reshaped j]);
dimsizes=permute(dimsizes,[dims_reshaped j]);

%Set up Interpolant
dimsizes(end)=length(dim.out.info.requested);
field_out=zeros(dimsizes);

%Expand to arbirary rank 
field=expand_field(field);

%Set up interpolant
if isempty(interpolant)
   if isempty(ps)
      interpolant=initialize_interpolant(dim_in,dim_out);
   else
      interpolant=initialize_interpolant(dim_in,dim_out,a,b,ps);
   end
end

%Reshape field back to original order

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Expand a field to an arbitrary matrix rank of 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function field=expand_field(field)

dimsizes=size(field);
padding=4-length(dimsizes);
if padding>0
   operator='(';
   for i=1:padding
      operator=cat(1,operator,':,');
   end
   eval(['field',operator,':)=field']);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Returns dimension information from spec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function info=dimension_info(specification,dimension,dictionary)

info=strings(0);
names=fieldnames(specification.coordinate.axis_entry);
for i=1:length(names)
   if strcmp(names{i},dimension)
      info=eval(['specification.coordinate.axis_entry.',names{i}]);
   end
end

if isempty(info)
   info=translate_cesm(dictionary,dimension,'Dimension'); 
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Returns the index of the matching field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ind=return_index(structure,field)
 
names=fieldnames(structure);
for i=1:length(names)
   if strcmp(field,names{i})
      ind=i;
   end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parses a string into individual words
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function parsed=parse_string(in)

parsed = strings(0);
while (in ~= "")
   [token,in] = strtok(in);
   parsed = [parsed ; token];
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Translate variable details, make actionable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [variable_out]=translate_cesm_variable(in)

in_components=fieldnames(in);
variables_components={'var';'operation';'eval'};
variable_out=struct('var',strings(0),'operation',strings(0),'eval',strings(0));

for i=1:length(in_components)
   for j=1:length(variables_components)
      if contains(in_components{i}, variables_components{j})
         eval(['variable_out.',variables_components{j},'=cat(1,variable_out.',variables_components{j},',in.',in_components{i},');']) 
      end
   end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Scan dictionary and return field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=translate_cesm(dictionary,in,field)

dictionary_entries=fieldnames(eval(['dictionary.',field]));
for i=1:length(dictionary_entries)
   if strcmp(dictionary_entries{i},in)
      out=eval(['dictionary.',field,'.',dictionary_entries{i}]);
   end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Gather required global attributes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function globals=global_attributes(required_attributes,CV_file,output_file)
   
globals=struct;
for i=1:length(required_attributes)
   globals(i).name=required_attributes{i};
   for j=1:length(CV_file)
   end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Remove entries from a list
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function list=remove_entries(list,entries)

delete_list=[];
for i=1:length(entries)
   for j=1:length(list)
      if strcmp([list(j).folder,'/',list(j).name],entries{i})
         delete_list=cat(1,delete_list,j);
      end
   end
end
list(delete_list)=[];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load file and decode json
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function decoded=json_load(file)
fid = fopen(file);
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
decoded=jsondecode(str);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load specification files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [specification,specification_files]=load_specification_files(file_preamble)

specification=struct;
specs={'CV';'fx';'grids';'coordinate';'formula_terms'};
for i=1:length(specs)
   specification_files{i}=[file_preamble,specs{i},'.json'];
   eval(['specification.',specs{i},'=json_load(''',specification_files{i},''');']);
end

end


