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

   output_var_details=struct2cell(output_file.variable_entry);  
   %Output file for each variable
   for v=23%1:length(vars)
      tic
      disp(vars{v})

      local_var_spec=vars_info{v};
 
      local_var_spec.averaging=averaging_flags(local_var_spec.cell_methods);

      dir_input_main=[cmor_specification.case_output_dir,cmor_specification.case_name,'/'];

      %Realm, frequency information
      realm=translate_cesm(cesm_dictionary,local_var_spec.modeling_realm,'Realms');
      frequency=translate_cesm(cesm_dictionary,local_var_spec.frequency,'Frequency');
      variables_list=translate_cesm(cesm_dictionary,vars{v},'Variables');
  
      %Model output location
      dir_input=[dir_input_main,realm,'/proc/tseries/',frequency,'_1/',cmor_specification.case_name]; 

      %Gather variable info, special operations, scale_factors; or just load the variable 
      if length(variables_list)==1 & ~isstruct(variables_list) 
         variables.var{1}=variables_list;
         variables.operation=[];
         variables.eval=[];
      else
         [variables]=translate_cesm_variable(variables_list);
      end
      load_ps=0;
      save_ps=0;
      load_lev=0;

      %Grab global variables for the file
      globals=global_attributes(specification.CV.CV.required_global_attributes,specification.CV.CV,output,cmor_specification,specification,frequency,vars{v});
      globals=set_index_variables(globals,cmor_specification.cmor_case_name);
    
      %Load variable(s) 
      for vin=1:length(variables.var)
         var_files=dir([dir_input,'*.',variables.var{vin},'.*']);
         file_name=[var_files(1).folder,'/',var_files(1).name];

         date=ncread(file_name,'date');
         time=ncread(file_name,'time');
         dims=parse_string(local_var_spec.dimensions);

         %Get dimensions of input variable
         file_structure=ncinfo(file_name);
         for j=1:length(file_structure.Variables)
            if strcmp(file_structure.Variables(j).Name,variables.var{vin})
               for k=1:length(ncinfo(file_name).Variables(j).Dimensions)
                  dim_native_names{k}=ncinfo(file_name).Variables(j).Dimensions(k).Name;
                  dim_native(k)=k;
               end       
            end
         end 

         %Load dimension info, connect output dims to input dims
         for j=1:length(dims)
            dim{vin}{j}.native.name=translate_cesm(cesm_dictionary,dims{j},'Dimension');    
            dim{vin}{j}.native.ind=find(strcmp(dim_native,dim{vin}{j}.native.name),1,'first');
            dim_native(strcmp(dim{vin}{j}.native.name,dim_native_names))=[];
            dim_native_names(strcmp(dim{vin}{j}.native.name,dim_native_names))=[]; 
            dim{vin}{j}.native.value=ncread(file_name,dim{vin}{j}.native.name);
            dim{vin}{j}.out.name=dims{j};
            dim{vin}{j}.out.info=dimension_info(specification,dims{j},cesm_dictionary);
            dim{vin}{j}.interp_special='';
            
            %Assimilate MIP lev info
            if strcmp(dim{vin}{j}.out.info,'lev')
               dim{vin}{j}.out.info=dimension_info(specification,'standard_hybrid_sigma',cesm_dictionary);
            end
            
            %do we need to interpolate the vertical grid?
            if strcmp(dim{vin}{j}.native.name,'lev')                
               lev_dim=j;
               load_lev=1;
               load_ps=1;
               save_ps=1;
               if strcmp(dim{vin}{j}.out.name,'alevel')
                  dim{vin}{j}.out.name='lev';
                  dim{vin}{j}.interp='';
               end
               dim{vin}{j}.interp_special='vertical';
               dim{vin}{j}.interp='';
            end

            %Add consistent time units
            if strcmp(dim{vin}{j}.native.name,'time')
               dim{vin}{j}.out.info.units=ncreadatt(file_name,'time','units');
            end
         end

         %Is lev dimension being removed? Probably need to load fields
         if ~isempty(dim_native)
            if strcmp(dim_native_names,'lev')
               load_lev=1;
               load_ps=1;
            end 
         end
    
         %Load lev formula fields
         if load_lev==1  
            a=ncread(file_name,'hyam')*1e5;
            b=ncread(file_name,'hybm');
            ai=ncread(file_name,'hyai')*1e5;
            bi=ncread(file_name,'hybi');
            ilev=ncread(file_name,'ilev');
         end
 
         %Load variable
         var{vin}.native.value=ncread(file_name,variables.var{vin});
         var{vin}.native.name=variables.var{vin};

         %Load ps
         if load_ps==1
            var_files=dir([dir_input,'*.PS.*']);
            file_name=[var_files(1).folder,'/',var_files(1).name];
            ps.value=ncread(file_name,'PS');
            for j=1:length(vars)
               if strcmp(vars{j},'ps')
                  ps.info=vars_info{j};
               end
            end
             
        end 
      end

      %Inherit dimensions of input var
      var_out.dim=dim{1};
      ps_out=ps;

      %Note dimensions to be removed (must be empty to save variable)
      var_out.dim_names_to_remove=dim_native_names;
      var_out.dim_to_remove=dim_native;

      %Perform arithmetic operation
      %Eval converts a string to code and executes
      if isstruct('variables')
         if ~isempty(variables.eval) 
            eval(['var_out.native.value=',variable.eval{1},';'])
         end
      else
         var_out.native=var{1}.native;
      end

      %Special operations - omega to wa, age of air, integrate, max value
      %Need to calculate three based on ps-pres field
      if ~isempty(variables.operation)
         operation=variables.operation;
         switch operation
            case 'age_of_air'
               disp('calculating age of air')
               var_out=age_of_air(var_out,variables.p0,variables.lat0);
            case 'integral'
               var_out=calculate_integral(var_out,variables.axis,ps_out,a,b);
            case 'max_value'
               var_out=calculate_maximum(var_out,variables.axis);
            case 'omega_to_w'
               var_out.temp=var{2}.native.value;
               var_out=omega_to_w(var_out,ps_out,a,b);
         end
      end
 
      grid_label='gn';
      %Do any interpolation
      for j=1:length(var_out.dim) 
         if ~isempty(var_out.dim{j}.out.info.requested)
            if strcmp(var_out.dim{j}.interp_special,'vertical')
               var_out=interpolate_field(var_out.native.value,j,var_out.dim{j},a,b,ps);
            else
               var_out=interpolate_field(var_out.native.value,j,var_out.dim{j});
               grid_label='gr';
               if load_ps==1
                  ps_out.value=interpolate_field(ps_out.value,j,var_out.dim{j});
               end
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
               ps_out.value=mean(ps_out.value,j);
               grid_label=cat(1,grid_label,'z');
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
            var_out.dim{j}.out.bnds= ncread(file_name,'time_bnds');
         else
            var_out.dim{j}.out.bnds = create_bnds(var_out.dim{j}.out);
         end
      end
      if do_ab~=0
         var_out.dim{do_ab}.out.formula.a_bnds.value = create_bnds(a,ai);
         var_out.dim{do_ab}.out.formula.a_bnds.info = get_formula_variable(specification.formula_terms.formula_entry,'a_bnds');
         var_out.dim{do_ab}.out.formula.b_bnds.value = create_bnds(b,bi);
         var_out.dim{do_ab}.out.formula.b_bnds.info = get_formula_variable(specification.formula_terms.formula_entry,'b_bnds');
         var_out.dim{do_ab}.out.formula.a.value=a;
         var_out.dim{do_ab}.out.formula.a.info = get_formula_variable(specification.formula_terms.formula_entry,'a');
         var_out.dim{do_ab}.out.formula.b.value=b;
         var_out.dim{do_ab}.out.formula.b.info = get_formula_variable(specification.formula_terms.formula_entry,'b');
      end
      var_out.info=output_var_details{v};

      %Save to file
      for j=1:length(globals)
         if strcmp(globals(j).name,'source_id')
            id_index=j;
         end
      end

      %Create directory
      dir_output=[cmor_specification.cmor_output_dir,cmor_specification.case_name,'/postprocess/output/',output,'/',vars{v},'/',grid_label,'/'];
      if ~exist(dir_output)
         mkdir(dir_output)
      end

      outfile=[dir_output,vars{v},'_',output,'_',globals(id_index).value,'_',cmor_specification.cmor_experiment,...
              '_',cmor_specification.cmor_case_name,'_',grid_label,'_',cmor_specification.case_dates,'.nc'];

      if save_ps==1
         create_netcdf(var_out,outfile,globals,ps_out);
      else
         create_netcdf(var_out,outfile,globals);
      end 
     
      clear variables variable dim dim_native
     
      toc
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Convert omega to w
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data=omega_to_w(data,ps,a,b)

ps=repmat(ps,[ones(1,length(size(ps))) length(a)]);
a=permute(repmat(a(:),[1 size(ps)]),[2:length(size(ps))+1 1]);
b=permute(repmat(b(:),[1 size(ps)]),[2:length(size(ps))+1 1]);
size(a)
size(b)
size(ps)
pres=a+b.*ps;

%Reshape field back to original order
size(pres)
rho=pres./(287*data.temp);
data.native.value=rho*9.81*data.native.value;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate maximum value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data=calculate_maximum(data,dim_name)

dim_num=data.dim_to_remove(find(strcmp(data.dim_names_to_remove,dim_name),1,'first'));
data.native.value=squeeze(max(data.native.value,[],dim_num));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate integral
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data=calculate_integral(data,dim_name,varargin)

dim_num=data.dim_to_remove(find(strcmp(data.dim_names_to_remove,dim_name),1,'first'));

if strcmp(data.dim_names_to_remove,'lev')
   ps=varargin{1};
   a=varargin{2};
   b=varargin{3};

   ps=repmat(ps,[ones(length(size(ps)),1) length(a)]);
   a=permute(repmat(a(:),[1 size(ps)]),[2:length(size(ps)) 1]);
   b=permute(repmat(b(:),[1 size(ps)]),[2:length(size(ps)) 1]);

   pres=a+ps.*b;
   for i=1:size(ps,1)
      for j=1:size(ps,2)
         for t=1:size(data.native.value,4)
            data_out(i,j,t)=trapz(squeeze(pres(i,j,t,:)),dim.native.value(i,j,:,t)); 
         end
      end
   end
   data.native.value=data_out; 
else
   data.native.value=squeeze(trapz(data.dim{dim_num}.native.value,data.native.value,dim_num));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate age of air
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tracer=age_of_air(tracer,p_ref,y_ref)

aoa=zeros(size(tracer.native.value));

%Setup indexing for eval
dim_string='(';
for i=1:length(size(tracer.native.value))
   dim_name=tracer.dim{i}.native.name;
   switch dim_name
      case 'time'
         string_add=':';
         time_dim=i;
      case 'lev'
         string_add='p';
         lev_dim=i;
      case 'lat'
         string_add='y';
         lat_dim=i;
      case 'lon'
         string_add='x';
         lon_dim=i;
      otherwise
         error('Unrecognized dimension name, cannot compute age of air')
   end
   if i<length(size(tracer.native.value))
      string_add=[string_add,','];
   else
      string_add=[string_add,')'];
   end
   dim_string=[dim_string,string_add];
end
time_string=strrep(dim_string,':','t');
ref_string=strrep(dim_string,'x',':');
ref_string=strrep(ref_string,'y','y_ref');
ref_string=strrep(ref_string,'p','p_ref');

%Get reference value
[val,y_ref]=min(abs(tracer.dim{lat_dim}.native.value-y_ref));
[val,p_ref]=min(abs(tracer.dim{lev_dim}.native.value-p_ref));

tracer.native.value=movmean(tracer.native.value,12,time_dim);
reference_timeseries=eval(['squeeze(mean(tracer.native.value',ref_string,',',num2str(lon_dim),'));']);
reference_timeseries=repmat(reference_timeseries(:),[1 length(reference_timeseries)]);

nt=size(tracer.native.value,time_dim);

%Calculate aoa
for y=1:length(tracer.dim{lat_dim}.native.value)
   disp([sprintf('%2.2f',100*(y-1)/length(tracer.dim{lat_dim}.native.value)),'% complete'])
   for x=1:length(tracer.dim{lon_dim}.native.value)
      for p=1:length(tracer.dim{lev_dim}.native.value)
         timeseries=eval(['squeeze(tracer.native.value',dim_string,')']);
         timeseries=repmat(timeseries(:),[1 length(timeseries)])';
         [val,time_index]=min(abs(timeseries-reference_timeseries),[],2);
         eval(['aoa',dim_string,'=(time_index''-[1:nt])/12;']);
      end
   end
end

tracer.native.value=aoa;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get formula variable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function info=get_formula_variable(formula_terms,var)

names=fieldnames(formula_terms);
formula_terms=struct2cell(formula_terms);
for i=1:length(formula_terms)
   if strcmp(names{i},var)
      info=formula_terms(i);
   end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get index variables (r#i#f#p#)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create NetCDF output file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function create_netcdf(var_out,filename,globals,varargin)

lev_ind=[];
do_ps=0;
if nargin>3
   ps=varargin{1};
   ps_dim=[];
   do_ps=1;
end

if exist(filename,'file')==2
   delete(filename);
   disp('Overwriting existing file')
end
ncid=netcdf.create(filename,'NETCDF4');
varid = netcdf.getConstant('GLOBAL');
for i=1:length(globals)
   netcdf.putAtt(ncid,varid,globals(i).name,globals(i).value);
end

for i=1:length(var_out.dim)
   ndim(i)=length(var_out.dim{i}.out.value);
   if strcmp(var_out.dim{i}.out.name,'time')
      dim(i)=netcdf.defDim(ncid,var_out.dim{i}.out.info.out_name,netcdf.getConstant('NC_UNLIMITED'));
   else
      dim(i)=netcdf.defDim(ncid,var_out.dim{i}.out.info.out_name,ndim(i));
   end
   if do_ps==1
      if ~strcmp(var_out.dim{i}.out.name,'lev')
         ps_dim=cat(1,ps_dim,dim(i));
      else
         lev_ind=i;
      end
   end
end
dim(i+1)=netcdf.defDim(ncid,'bnds',2);

var=netcdf.defVar(ncid,var_out.info.out_name,'NC_FLOAT',dim(1:end-1));
attnames=fieldnames(var_out.info);
atts=struct2cell(var_out.info);
for i=1:length(atts)
   if ~isempty(atts{i}) & sum(strcmp(attnames{i},{'frequency';'realm';'dimensions';'type';'out_name'}))==0
      netcdf.putAtt(ncid,var,attnames{i},atts{i});
   end
end

if do_ps==1
   var_ps=netcdf.defVar(ncid,'ps','NC_FLOAT',ps_dim);
   attnames=fieldnames(ps.info);
   atts=struct2cell(ps.info);
   for i=1:length(atts)
      if ~isempty(atts{i}) & sum(strcmp(attnames{i},{'frequency';'realm';'dimensions';'type';'out_name'}))==0
         netcdf.putAtt(ncid,var_ps,attnames{i},atts{i});
      end
   end
   formula_names=fieldnames(var_out.dim{lev_ind}.out.formula);
   formula=struct2cell(var_out.dim{lev_ind}.out.formula);
   for i=1:length(formula_names)
      if sum(size(formula{i}.value)==1)==0
         formula_dims{i}=[dim(end) dim(lev_ind)];
      else
         formula_dims{i}=dim(lev_ind);
      end
      var_formula(i)=netcdf.defVar(ncid,formula_names{i},'NC_FLOAT',formula_dims{i});
      attnames=fieldnames(formula{i}.info{1});
      atts=struct2cell(formula{i}.info{1});
      for j=1:length(atts)
         if ~isempty(atts{j}) & sum(strcmp(attnames{j},{'dimensions';'type';'out_name'}))==0
            netcdf.putAtt(ncid,var_formula(i),attnames{j},atts{j})
         end;
      end
   end
end
for i=1:length(var_out.dim)
   var_dim(i)=netcdf.defVar(ncid,var_out.dim{i}.out.info.out_name,'NC_FLOAT',dim(i));
   var_dim_bnds(i)=netcdf.defVar(ncid,[var_out.dim{i}.out.info.out_name,'_bnds'],'NC_FLOAT',[dim(end) dim(i)]);
   attnames=fieldnames(var_out.dim{i}.out.info);
   atts=struct2cell(var_out.dim{i}.out.info);
   for j=1:length(atts)
       if ~isempty(atts{j}) & sum(strcmp(attnames{j},{'out_name';'must_have_bounds';'stored_direction';'type'}))==0
          netcdf.putAtt(ncid,var_dim(i),attnames{j},atts{j});
          netcdf.putAtt(ncid,var_dim_bnds(i),attnames{j},atts{j}); 
       end
   end   

end

netcdf.endDef(ncid);

start=zeros(length(var_out.dim),1);
ending=start;
for i=1:length(ending)
   ending(i)=length(var_out.dim{i}.out.value);
end
netcdf.putVar(ncid,var,start,ending,var_out.native.value);
if do_ps==1
   netcdf.putVar(ncid,var_ps,ps.value);
   for i=1:length(var_formula)
      netcdf.putVar(ncid,var_formula(i),formula{i}.value);
   end
end

for i=1:length(var_out.dim)
   netcdf.putVar(ncid,var_dim(i),var_out.dim{i}.out.value);
   netcdf.putVar(ncid,var_dim_bnds(i),var_out.dim{i}.out.bnds);
end

netcdf.close(ncid);


end

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

   dim.native.value=a+ps.*b;
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
variables_components={'var';'operation';'eval';'p0';'lat0';'axis'};
variable_out=struct('var',strings(0),'operation',strings(0),'eval',strings(0),'p0',[],'lat0',[],'axis',strings(0));

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
