%Script to CMORize CESM output
%Currently only functioning for the atmosphere component
%Assumes lat/lon grid, single variable timeseries (convert SE/MPAS to gridded output first)
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

cesm_globals=json_load(['cesm_global_attributes_Amon.json']);
cesm_globals_names=fieldnames(cesm_globals);

%Load/transform data, collect metadata, create output files, for each output spec
for i=4%1:length(output_specification_files)

   output_specification_file=[output_specification_files(i).folder,'/',output_specification_files(i).name];
   output_file=json_load(output_specification_file);

   output=strrep(output_specification_files(i).name,cmor_specification.cmor_specification_format,'');
   output=strrep(output,'.json','');

   vars=fieldnames(output_file.variable_entry);
   vars_info=struct2cell(output_file.variable_entry); 
  
   ps_init=0;

   %Output file for each variable
   for v=30%1:length(vars)

      local_var_spec=vars_info{v};
 
      local_var_spec.averaging=averaging_flags(local_var_spec.cell_methods);

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

      %Gather variable info, special operations, scale_factors; or just load the variable 
      if length(fieldnames(variables_list))>1
         [variable]=translate_cesm_variable(variables_list);
      else
         variable{1}=variables_list;
      end
         
      load_ps=0;

      %Grab global variables for the file
      globals=global_attributes(specification.CV.CV.required_global_attributes,specification.CV.CV,output,cmor_specification,specification,frequency,vars{v});

      %Load variable(s) 
      for vin=1:length(variable.var)
         var_files=dir([dir_input,'*.',variable.var{vin},'.*']);
         file_name=[var_files(1).folder,'/',var_files(1).name];

         dims=parse_string(local_var_spec.dimensions);
         for j=1:length(dims)
            dim{vin}{j}.native.name=translate_cesm(cesm_dictionary,dims{j},'Dimension');    
            dim{vin}{j}.native.value=ncread(file_name,dim{vin}{j}.native.name);
            dim{vin}{j}.out.name=dims{j};
            dim{vin}{j}.out.info=dimension_info(specification,dims{j},cesm_dictionary);
            dim{vin}{j}.interp_special='';
            
            %is there a requested grid for this dimension?
            if isstruct(dim{vin}{j}.out.info)
               dim{vin}{j}.interp=dim{vin}{j}.out.info.requested;
            elseif strcmp(dim{vin}{j}.out.info,'lev')
               dim{vin}{j}.out.info=struct('axis','Z','positive','down','long_name','hybrid sigma pressure coordinate','standard_name','atmosphere_hybrid_sigma_pressure_coordinate','formula','p = ap + b*ps','formula_terms','ap: ap b: b ps: ps','valid_max','','valid_min','');
            end
 
            %do we need to interpolate the vertical grid?
            if strcmp(dim{vin}{j}.native.name,'lev')
               if strcmp(dim{vin}{j}.out.name,'alevel')
                  dim{vin}{j}.out.name='lev';
                  dim{vin}{j}.interp='';
               end
               dim{vin}{j}.interp_special='vertical';
               load_ps=1;
               a=ncread(file_name,'hyam')*1e5;
               b=ncread(file_name,'hybm');
               ai=ncread(file_name,'hyai')*1e5;
               bi=ncread(file_name,'hybi');
               ilev=ncread(file_name,'ilev');
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
      %Eval converts a string to code and executes
      if ~isempty(variable.eval) 
         eval(['var_out.native.value=',variable.eval{1},';'])
      else
         var_out.native=var{1}.native;
      end

      %Special operations - omega to wa, age of air, integrate, max value
      %%%Todo        
 
      %Do any interpolation
      for j=1:length(var_out.dim) 
         if strcmp(var_out.dim{j}.interp,'interpolate')
            if strcmp(var_out.dim{j}.interp_special,'vertical')
               var_out=interpolate_field(var_out.value,j,var_out.dim{j},a,b,ps);
            else
               var_out=interpolate_field(var_out.value,j,var_out.dim{j});
            end
            var_out.dim{j}.out.value=dim{vin}{j}.interp;
         else
            var_out.dim{j}.out.value=ncread(file_name,var_out.dim{j}.native.name);
         end
      end

      %Averaging, preserves singular dimensions
      for j=1:length(dims)
         for k=1:length(local_var_spec.averaging)
            if strcmp(local_var_spec.averaging{k},dims{j}) && ~strcmp(local_var_spec.averaging{k},'time')
               var_out.native.value=mean(var_out.native.value,j);
            end
         end
      end 

      do_ab=0;
      %Bnd variables
      for j=1:length(var_out.dim)
         if strcmp(var_out.dim{j}.out.name,'lev')
            var_out.dim{j}.out.bnds = create_bnds(var_out.dim{j}.out,ilev);
            do_ab=j;
         elseif strcmp(var_out.dim{j}.out.name,'time')
            var_out.dim{j}.out.bnds=ncread(file_name,'time_bnds');
         else
            var_out.dim{j}.out.bnds = create_bnds(var_out.dim{j}.out);
         end
      end
      if do_ab~=0
         var_out.dim{do_ab}.out.a_bnds = create_bnds(a,ai);
         var_out.dim{do_ab}.out.b_bnds = create_bnds(b,bi);
         var_out.dim{do_ab}.out.a=a;
         var_out.dim{do_ab}.out.b=b;
      end

      %Save to file
           
   end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create NetCDF output file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Determine which dims need to be averaged
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
      
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create bounding variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bnds=create_bnds(dim_in,varargin)

if isstruct(dim_in)
   dim=dim_in.value;
else
   dim=dim_in;
end
bnds=zeros(2,length(dim));

%Passed interface field (lev, a, b, etc.)
if nargin>1
   dim_interface=varargin{1};
   for i=1:size(bnds,2)
      bnds(1,i)=dim_interface(i);
      bnds(2,i)=dim_interface(i+1);
   end
else
   %Valid dim bounds exist (lat, lon)
   if ~isempty(dim_in.info.valid_max)
      maxbnd=str2num(dim_in.info.valid_max);
      minbnd=str2num(dim_in.info.valid_min);

      if xor(dim(1)==minbnd,dim(end)==maxbnd)
         if dim(1)==minbnd
            loop_start=1;
            loop_end=size(bnds,2)-1;
            bnds(1,size(bnds,2))=dim(end);
            bnds(2,size(bnds,2))=maxbnd;
         else
            loop_start=2;
            loop_end=size(bnds,2);
            bnds(1,1)=minbnd;
            bnds(2,1)=dim(1);
         end
         for j=loop_start:loop_end
            bnds(1,j)=dim(j);
            bnds(2,j)=dim(j+1);
         end
      else 
         bnds(1,1)=minbnd;
         bnds(2,end)=maxbnd;

         bnds(2,1)=bnds(1,1)+(dim(2)-dim(1))/2;
         bnds(1,end)=bnds(2,end)-(dim(end)-dim(end-1))/2;
         for j=2:size(bnds,2)-1
            bnds(1,j)=dim(j)-(dim(j)-dim(j-1))/2;
            bnds(2,j)=dim(j)+(dim(j+1)-dim(j))/2;
         end
      end
   else
      error('Unable to construct bnds variable - no valid min/max in mip spec and no interface field')
   end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Interpolates along specified dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function field_out=interpolate_field(field,interp_dim,dim,varargin)

%Optional: load in surface pressure, a/b's for vertical interp, generate pressure array
if nargin>3
   a=varargin{1};
   b=varargin{2};
   ps=varargin{3};

   ps=repmat(ps,[ones(length(size(ps)),1) length(a)]);
   a=permute(repmat(a(:),[1 size(ps)]),[2:length(size(ps)) 1]);
   b=permute(repmat(b(:),[1 size(ps)]),[2:length(size(ps)) 1]);

   dim.native.value=1e5*a+ps.*b;
   dim_reshape=1;
else
   dim_reshape=0;
end

dimsizes=size(field);
dims=1:length(dimsizes);

%Reshape field to have j on right, arbitrary 4D array
dims_reshaped=dims;
dims_reshaped(interp_dim)=[];
field=permute(field,[dims_reshaped interp_dim]);
rank_permute=length(size(field));

%We need to reshape so the pressure array matches
if dim_reshape==1
   dim.native.value=permute(dim.native.value,[dims_reshaped interp_dim]);
end
dimsizes=permute(dimsizes,[dims_reshaped interp_dim]);

%Set up interpolant
dimsizes(end)=length(dim.out.info.requested);
field_out=zeros(dimsizes);

%Expand to arbirary rank with interp on right
field=expand_field(field);

%Interpolate
for i=1:size(field,1)
   for j=1:size(field,2)
      parfor k=1:size(field,3)
         if dim_reshape==1
            field_out(i,j,k,:)=interp1(squeeze(dim.native.value(i,j,k,:)),...
                                  squeeze(field(i,j,k,:)),dim.out.info.requested);
         else
            field_out(i,j,k,:)=interp1(dim.native.value,...
                                  squeeze(field(i,j,k,:)),dim.out.info.requested);
         end
      end
   end
end

%Collapse any extra singular dimensions, ensure any preexisting singular remains
field_out=collapse_field(field_out);

%Reshape field back to original order
j_index=length(size(field_out));
if interp_dim==1
   dims_reshaped=[j_index 1:j_index-1];
elseif interp_dim==j_index
   dims_reshaped=[1:j_index];
else
   dims_reshaped=[1:interp_dim-1 j_index interp_dim:j_index-1];
end
field_out=permute(field_out,dims_reshaped); 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Collapse a field to a specified matrix rank
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function field=collapse_field(field,rank)

dimsizes=length(size(field));

if length(field)>rank

   padding=rank;
   operator='(';
   for i=1:padding-1
      operator=cat(1,operator,':,');
   end
   operator=cat(operator,':)');
   
   singulars=length(dimsizes)-rank; 
   singular_operator='(';
   for i=1:singulars
      singular_operator=cat(1,singular_operator,'1,');
   end
   for i=1:4-singulars+1
     singular_operator=cat(1,singular_operator,':,');
   end
   singular_operator=cat(1,singular_operator,':)');

   eval(['field',operator,'=field',singular_operator,';']);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Expand a field to an arbitrary matrix rank of 4
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
function globals=global_attributes(required_attributes,CV_file,output,cmor_specification,specification,frequency,var)
   
globals=struct;
cesm_globals=json_load(['cesm_global_attributes_',output,'.json']);
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
