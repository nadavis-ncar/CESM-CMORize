%Script to CMORize CESM output
%Currently only functioning for the atmosphere component
%Assumes lat/lon grid, single variable timeseries (convert SE/MPAS to gridded output first)
%nadavis@ucar.edu

clear all

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
cesm_globals=json_load(['cesm_global_attributes.json']);
cesm_globals_names=fieldnames(cesm_globals);

%Load/transform data, collect metadata, create output files, for each output spec
for i=5   

   output_specification_file=[output_specification_files(i).folder,'/',output_specification_files(i).name];
   output_file=json_load(output_specification_file);

   output=strrep(output_specification_files(i).name,cmor_specification.cmor_specification_format,'');
   output=strrep(output,'.json','');

   vars=fieldnames(output_file.variable_entry);
   vars_info=struct2cell(output_file.variable_entry); 

   output_var_details=struct2cell(output_file.variable_entry);  
  
   version=['v',datestr(now,'yyyymmdd')];
   tem=0;

   %Output file for each variable
   for v=1:length(vars)
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
      if ~isstruct(variables_list) 
         variables.var{1}=variables_list;
         variables.operation=[];
         variables.eval=[];
         variables.option=[];
      else
         [variables]=translate_cesm_variable(variables_list);
      end
      load_ps=0;
      save_ps=0;
      load_lev=0;

      %Grab global variables for the file
      globals=global_attributes(specification.CV.CV.required_global_attributes,specification.CV.CV,output,cmor_specification,specification,frequency,vars{v},cesm_globals);
      globals=set_index_variables(globals,cmor_specification.cmor_case_name);
   
      %For ambiguous variables
      if ~isempty(variables.option)
         for j=1:length(variables.option)
            if contains(var_info{v}.comment,variables.option(j).string)
               variables.var{1}=variables.option(j).name;
            end
         end
      end
 
      %Set TEM variable(s)
      if isempty(variables.var) & strcmp(variables.operation,'TEM') & tem==0
         disp('TEM variable requested, now calculating all TEM output variables')
         variables.var{1}='Uzm';
         variables.var{2}='Vzm';
         variables.var{3}='Wzm';
         variables.var{4}='VTHzm';
         variables.var{5}='UVzm';
         variables.var{6}='UWzm';
         variables.var{7}='THzm'; 
         tem=1;
      end

      %Either we're not outputting the variable or we've already done it (TEM)
      if ~isempty(variables.var)

	      %Assemble file list(s)
	      for vin=1:length(variables.var)
	      	var_files=dir([dir_input,'*.',variables.var{vin},'.*']);
	      	for i=1:length(var_files)
	            file_list{vin}{i}=[var_files(1).folder,'/',var_files(1).name];
	         end
	      end

      	while ~isempty(file_list{1})

		      %Load variable(s) 
		      for vin=1:length(variables.var)
		         
		         file_name=[file_list{vin}(1).folder,'/',file_list{vin}(1).name];

		         date=ncread(file_name,'date');
		         time=ncread(file_name,'time');
		         dims=parse_string(local_var_spec.dimensions);

		         %Get dimensions of input variable
		         %Deal with useless singular dimension "zlon"
		         remove_list=[];
		         file_structure=ncinfo(file_name);
		         for j=1:length(file_structure.Variables)
		            if strcmp(file_structure.Variables(j).Name,variables.var{vin})
		               for k=1:length(ncinfo(file_name).Variables(j).Dimensions)
		                  if strcmp(ncinfo(file_name).Variables(j).Dimensions(k).Name,'zlon') 
		                     remove_list=k;
		                  else
		                     dim_native_values{k}=ncread(file_name,ncinfo(file_name).Variables(j).Dimensions(k).Name);
		                     dim_native_names{k}=ncinfo(file_name).Variables(j).Dimensions(k).Name;
		                     dim_native(k)=k;
		                  end
		               end       
		            end
		         end
		         dim_native_names(remove_list)=[];
		         dim_native(remove_list)=[];
		         dim_native=dim_native-length(remove_list);
		         dim_native_names_list=dim_native_names;

		         %Load dimension info, connect output dims to input dims
		         for j=1:length(dims)
		            dim{vin}{j}.native.name=translate_cesm(cesm_dictionary,dims{j},'Dimension');   
		            %Single-level/subdomain 
		            if strcmp(dim{vin}{j}.native.name,'no_match')
		               dim{vin}{j}.native.name=translate_cesm(cesm_dictionary,eval(['specification.coordinate.axis_entry.',dims{j},'.out_name']),'Dimension');
		            end
		            %Already interpolated?
		            if strcmp(dim{vin}{j}.native.name,'lev')
		               for k=1:length(dim_native_names)
		                  if strcmp(dim_native_names{k},'ilev')
		                     dim{vin}{j}.native.name='ilev';
		                     load_lev=1;
		                  end
		               end
		            end
		            dim{vin}{j}.native.ind=find(strcmp(dim_native,dim{vin}{j}.native.name),1,'first');
		            dim_native(strcmp(dim{vin}{j}.native.name,dim_native_names))=[];
		            dim_native_names(strcmp(dim{vin}{j}.native.name,dim_native_names))=[]; 
		            dim{vin}{j}.native.value=ncread(file_name,dim{vin}{j}.native.name);
		            dim{vin}{j}.out.name=dims{j};
		            dim{vin}{j}.out.info=dimension_info(specification,dims{j},cesm_dictionary);
		            dim{vin}{j}.interp_special='';
		            dim{vin}{j}.interp=[];
		            
		            %Assimilate MIP lev info for native output
		            if strcmp(dim{vin}{j}.out.info,'lev')
		               dim{vin}{j}.out.info=dimension_info(specification,'standard_hybrid_sigma',cesm_dictionary);
		            end
		            if (strcmp('lev',dim{vin}{j}.native.name) | strcmp('ilev',dim{vin}{j}.native.name)) ...
		                   & sum(strcmp(fieldnames(specification.coordinate.axis_entry),dim{vin}{j}.out.name))>0
		               dim{vin}{j}=get_requested_axis_value(dim{vin}{j},specification.coordinate.axis_entry);
		            end     
		 
		            if strcmp(dim{vin}{j}.native.name,'ilev') | strcmp(dim{vin}{j}.native.name,'lev')
		               dim{vin}{j}.native.value=dim{vin}{j}.native.value*100;
		            end
		             
		            %do we need to interpolate the vertical grid?
		            if strcmp(dim{vin}{j}.native.name,'lev')                
		               lev_dim=j;
		               load_lev=1;
		               load_ps=1;
		               if isempty(dim{vin}{j}.interp)
		                  save_ps=1;
		               end
		               if strcmp(dim{vin}{j}.out.name,'alevel')
		                  dim{vin}{j}.out.name='lev';
		                  dim{vin}{j}.interp='';
		               elseif strcmp(dim{vin}{j}.out.name,'alevhalf')
		                  dim{vin}{j}.out.name='lev';
		                  dim{vin}{j}.interp='';
		               end
		               dim{vin}{j}.interp_special='vertical';
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
		         if strcmp(variables.operation,'TEM') | ~isempty(remove_list) 
		            var{vin}.native.value=squeeze(var{vin}.native.value);
		         end
		         var{vin}.native.name=variables.var{vin};

		         %Load ps
		         if load_ps==1
		            var_files=dir([dir_input,'*.PS.*']);
		            file_name_ps=[var_files(1).folder,'/',var_files(1).name];
		            ps_out.value=ncread(file_name_ps,'PS');
		            for j=1:length(vars)
		               if strcmp(vars{j},'ps')
		                  ps_out.info=vars_info{j};
		               end
		            end
		         end 
      	   end

		      %Inherit dimensions of input var
		      var_out.dim=dim{1};

		      %Note dimensions to be removed (must be empty to save variable)
		      var_out.dim_names_to_remove=dim_native_names;
		      var_out.dim_to_remove=dim_native;

		      %Perform arithmetic operation
		      %Eval converts a string to code and executes
		      if isstruct(variables)
		         if ~isempty(variables.eval) 
		            eval(['var_out.native.value=',variables.eval{1},';'])
		         else
		            var_out.native=var{1}.native;
		            if length(var)>1
		               var_out.native_secondary=var{2}.native;
		            end
		         end 
		      else
		         var_out.native=var{1}.native;
		      end

		      for i=1:length(dim_native_names_list)
		         var_out.native_dim(i).name=dim_native_names_list{i};
		         var_out.native_dim(i).value=dim_native_values{i};
		      end

		      %Special operations - omega to w, age of air, integrate, max value
		      if ~isempty(variables.operation)
		         operation=variables.operation;
		         switch operation
		            case 'age_of_air'
		               disp('calculating age of air')
		               var_out=age_of_air(var_out,variables.p0,variables.lat0);
		            case 'sum'
		               var_out=calculate_sum(var_out,variables.axis,ps_out.value,a,b);
		            case 'ozone_integral'
		               var_out=convert_units(var_out,ps_out.value,a,b,'molar_mixing_ratio_to_DU');
		               if strcmp(fieldnames(var_out),'native_secondary')
		                  var_out=calculate_integral(var_out,ps_out.value,a,b,'tropopause');
		               else
		                  var_out=calculate_integral(var_out,ps_out.value,a,b);
		               end
		            case 'burden'
		               var_out=convert_units(var_out,ps_out.value,a,b,'number_density_to_molar_mixing_ratio');
		               var_out=calculate_integral(var_out,ps_out.value,a,b);
		            case 'max_value'
		               var_out=calculate_maximum(var_out,variables.axis);
		            case 'omega_to_w'
		               var_out.temp=var{2}.native.value;
		               var_out=omega_to_w(var_out,ps_out.value,a,b);
		            case 'TEM'
		               var_out=calculate_tem(var_out);
		         end
		      end

		      %Need to generalize to all cases
		      if ~isempty(dim_native)
		         dim_num_shift=1;
		      else
		         dim_num_shift=0;
		      end

		      %Ensure dimension ordering is correct
		      permute_order=[];
		      for j=1:length(var_out.dim)
		         if strcmp(var_out.dim{j}.native.name,dim_native_names_list{j})
		            permute_order=[permute_order j];
		         else
		            for k=1:length(dim_native_names_list)
		               if strcmp(var_out.dim{j}.native.name,dim_native_names_list{k})
		                  permute_order=[permute_order k];
		               end
		            end
		         end
		      end

		      %Only attempt reordering if necessary - typically for single pressure-level output
		      if length(permute_order) == length(size(var_out.native.value))    
		         var_out.native.value=permute(var_out.native.value,permute_order);
		      end

		      %Do any interpolation
		      for j=1:length(var_out.dim) 
		         if ~isempty(var_out.dim{j}.interp)
		            %Vertical interpolation on native levels, vs. on gridded dimensions
		            if strcmp(var_out.dim{j}.interp_special,'vertical') & load_ps==1
		               var_out.native.value=interpolate_field(var_out.native.value,j+dim_num_shift,var_out.dim{j},NaN,a,b,ps_out.value);
		            else
		               var_out.native.value=interpolate_field(var_out.native.value,j+dim_num_shift,var_out.dim{j},NaN);
		               if load_ps==1
		                  ps_out.value=interpolate_field(ps_out.value,j+dim_num_shift,var_out.dim{j});
		               end
		            end
		            var_out.dim{j}.out.value=var_out.dim{j}.interp;
		         else
		            var_out.dim{j}.out.value=ncread(file_name,var_out.dim{j}.native.name);
		         end
		      end

		      %Averaging
		      for j=1:length(local_var_spec.averaging)
		         if ~strcmp(local_var_spec.averaging{j},'time') & strcmp(local_var_spec.averaging{j},'longitude') & length(remove_list)==0
		            var_out.native.value=squeeze(mean(var_out.native.value,1,'omitnan'));
		            var_out.native.value(isnan(var_out.native.value))=str2num(output_file.Header.missing_value);
		            ps_out.value=mean(ps_out.value,1);
		         end
		      end

		      %Replace all NaNs with missing values  
		      var_out.native.value(isnan(var_out.native.value))=str2num(output_file.Header.missing_value);
		 
		      %Ensure dimension directions are consistent with file specification
		      var_out=stored_direction(var_out);   

		      do_ab=0;
		      %Bnd variables
		      for j=1:length(var_out.dim)
		         if strcmp(var_out.dim{j}.out.name,'lev')
		            var_out.dim{j}.out.bnds = create_bnds(var_out.dim{j}.out,ilev);
		            do_ab=j;
		         elseif strcmp(var_out.dim{j}.out.name,'time')
		            var_out.dim{j}.out.bnds= ncread(file_name,'time_bnds');
		         elseif length(var_out.dim{j}.out.value)>1
		            var_out.dim{j}.out.bnds = create_bnds(var_out.dim{j}.out);
		         end
		      end
	 
		      %Bnd variables (vertical)
		      if do_ab~=0
		         var_out.dim{do_ab}.out.formula.a_bnds.value = create_bnds(a,ai);
		         var_out.dim{do_ab}.out.formula.a_bnds.info = get_formula_variable(specification.formula_terms.formula_entry,'a_bnds');
		         var_out.dim{do_ab}.out.formula.b_bnds.value = create_bnds(b,bi);
		         var_out.dim{do_ab}.out.formula.b_bnds.info = get_formula_variable(specification.formula_terms.formula_entry,'b_bnds');
		         var_out.dim{do_ab}.out.formula.a.value=a;
		         var_out.dim{do_ab}.out.formula.a.info = get_formula_variable(specification.formula_terms.formula_entry,'a');
		         var_out.dim{do_ab}.out.formula.b.value=b;
		         var_out.dim{do_ab}.out.formula.b.info = get_formula_variable(specification.formula_terms.formula_entry,'b');
		         if var_out.fim{do_ab}.flip==1
		            var_out.dim{do_ab}.out.formula.a_bnds.value=flip(var_out.dim{do_ab}.out.formula.a_bnds.value,2);
		            var_out.dim{do_ab}.out.formula.b_bnds.value=flip(var_out.dim{do_ab}.out.formula.b_bnds.value,2);
		            var_out.dim{do_ab}.out.formula.a.value=flip(var_out.dim{do_ab}.out.formula.a.value,1);
		            var_out.dim{do_ab}.out.formula.b.value=flip(var_out.dim{do_ab}.out.formula.b.value,1);
		         end
		      end

		      %Set info for NetCDF
		      var_out.info=output_var_details{v};

		      %MIP-compatible model name
		      for j=1:length(globals)
		         if strcmp(globals(j).name,'source_id')
		            id_index=j;
		         end
		      end

		      %Grid labeling
		      grid_label='gn';
		      if strcmp(output(end),'Z')
		         grid_label=[grid_label,'z'];
		      end

		      %Save to file      
		      if strcmp(variables.operation,'TEM')
		         tem_vars=struct2cell(var_out.tem);
		         for k=1:length(fieldnames(var_out.tem))
		            var_out.native.value=tem_vars{k};
		            var_out.info=output_var_details{find(strcmp(vars,tem_out{k},1,'first'))};
		            dir_output=[cmor_specification.cmor_output_dir,cmor_specification.case_name,'/postprocess/output/',output,'/',tem_out{k},'/',grid_label,'/',version,'/'];
		            if ~exist(dir_output)
		               mkdir(dir_output)
		            end
		            outfile=[dir_output,tem_out{k},'_',output,'_',globals(id_index).value,'_',cmor_specification.cmor_experiment,...
		           '_',cmor_specification.cmor_case_name,'_',grid_label,'_',file_name(end-15:end-3),'.nc'];
		            create_netcdf(var_out,outfile,globals);
		         end
		      else
		         dir_output=[cmor_specification.cmor_output_dir,cmor_specification.case_name,'/postprocess/output/',output,'/',vars{v},'/',grid_label,'/',version,'/'];
		         if ~exist(dir_output)
		            mkdir(dir_output)
		         end
		         outfile=[dir_output,vars{v},'_',output,'_',globals(id_index).value,'_',cmor_specification.cmor_experiment,...
		              '_',cmor_specification.cmor_case_name,'_',grid_label,'_',file_name(end-15:end-3),'.nc'];
		         if save_ps==1
		            create_netcdf(var_out,outfile,globals,ps_out);
		         else
		            create_netcdf(var_out,outfile,globals);
		         end
		      end 
	     
	     		%Remove files from list 
		      for j=1:length(file_list)
		      	file_list{j}(1)=[];
		      end
		      clear variables dim dim_native var_out ps_out
		     
		      toc

         end
      end
   end
end

delete(gcp('nocreate'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Enforces dimension rules
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function var_out=stored_direction(var_out)

for i=1:length(var_out.dim)
   if length(var_out.dim{i}.out.value)>1
   flip_dim=0;
   dim_increasing=(var_out.dim{i}.out.value(2)-var_out.dim{i}.out.value(1))>0;
   if dim_increasing==1
      if strcmp(var_out.dim{i}.out.info.stored_direction,'decreasing')
         flip_dim=1;
      end
   else
      if strcmp(var_out.dim{i}.out.info.stored_direction,'increasing')
         flip_dim=1;
      end
   end
   if flip_dim==1
      var_out.native.value=flip(var_out.native.value,i);
      var_out.dim{i}.out.value=flip(var_out.dim{i}.out.value,i);
      var_out.dim{i}.flip=1;
   else
      var_out.dim{i}.flip=0;
   end
   end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculates TEM output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function var_out=calculate_tem(var_out)

g=9.81;
a=6371e3;
omega=7.292e-5;
H=6800;

lat=var_out.dim{1}.native.value;
lev=var_out.dim{2}.native.value;
time=var_out.dim{3}.native.value;

lev=permute(repmat(lev(:),[1 length(lat) length(time)]),[2 1 3]);
rho=(1.5/1e5)*lev;
z=log(lev/1e5)*-H;
lat=repmat(lat(:),[1 length(lev) length(time)]);
f=2*omega*sind(lat);
cosmat=cosd(lat);
phi=deg2rad(lat);

uzm=variables.var{1};
vzm=variables.var{2};
wzm=variables.var{3};
vthzm=variables.var{4};
uvzm=variables.var{5};
uwzm=variables.var{6};
thzm=variables.var{7};

[dummy,dudz]=gradient(uzm,squeeze(phi(:,1,1)),squeeze(z(1,:,1)));
[dudy,dummy]=gradient(uzm.*cosmat,squeeze(phi(:,1,1)),squeeze(z(1,:,1)));
dudy=dudy.*(1./(a*cosmat));
[dummy,dtheta_dz]=gradient(thzm,squeeze(phi(:,1,1)),squeeze(z(1,:,1)));

var_out.tem.epfy=rho.*a.*cosmat.*(dudz.*vthzm./dtheta_dz-uvzm);
var_out.tem.epfz=rho.*a.*cosmat.*((f-dudy).*vthzm./dtheta_dz-uwzm);

[dFydy,dummy]=gradient(var_out.tem.epfy.*cosmat,squeeze(phi(:,1,1)),squeeze(z(1,:,1)));
dFydy=dFydy.*(1./(a*cosmat));
[dummy,dFzdz]=gradient(var_out.tem.epfz,squeeze(phi(:,1,1)),squeeze(z(1,:,1)));
var_out.tem.utendepfd=(1./(rho.*a.*cosmat)).*(dFydy+dFzdz);

%residual circulation
[dummy,dTdz]=gradient(rho.*vthzm./dtheta_dz,squeeze(phi(:,1,1)),squeeze(z(1,:,1)));
[dTdy,dummy]=gradient(cosmat.*vthzm./dtheta_dz,squeeze(phi(:,1,1)),squeeze(z(1,:,1)));
var_out.tem.vtem=vzm-(dTdz./rho);
var_out.tem.wtem=wzm+(1./(a*cosmat)).*dTdy;

var_out.tem.epfy=var_out.tem.epfy./rho;
var_out.tem.epfz=var_out.tem.epfz./rho;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get requested axis values or single value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dim=get_requested_axis_value(dim,axis_entries)

axes=struct2cell(axis_entries);
axes_entry=find(strcmp(fieldnames(axis_entries),dim.out.name),1,'first');
if ~isempty(axes{axes_entry}.requested)
   requested=axes{axes_entry}.requested;
   for i=1:length(requested)
      dim.interp(i)=str2double(requested{i});
   end
   dim.interp=dim.interp(:);
elseif ~isempty(axes{axes_entry}.value)
   dim.interp=str2double(axes{axes_entry}.value);
else
   error('No value found for axis requesting specified values')
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Expand surface pressure to full pressure field with hybrid coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pres=calculate_pressure(ps,a,b)

ps_mat=permute(repmat(ps,[ones(1,length(size(ps))) length(a)]),[1 2 4 3]);
a=permute(repmat(a(:),[1 size(ps)]),[2 3 1 4]);
b=permute(repmat(b(:),[1 size(ps)]),[2 3 1 4]);
pres=a+b.*ps_mat;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Convert units 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data=convert_units(data,ps,a,b,conversion)

if strcmp(conversion,'number_density_to_molar_mixing_ratio')
   pres=calculate_pressure(ps,a,b);

   rho=pres./(287*data.native_secondary.value);
   n_air=(6.02e23/0.028)*rho;

   data.native.value=(100^3)*data.native.value./n_air;
elseif strcmp(conversion,'molar_mixing_ratio_to_DU')
   data.native.value=data.native.value*(287*273.15/1e5);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Convert omega to w
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data=omega_to_w(data,ps,a,b)

pres=calculate_pressure(ps,a,b);
data.native.value=9.81*(pres./(287*data.temp)).*data.native.value;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate maximum value along arbitrary dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data=calculate_maximum(data,dim_name)

dim_num=data.dim_to_remove(find(strcmp(data.dim_names_to_remove,dim_name),1,'first'));
data.native.value=squeeze(max(data.native.value,[],dim_num));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate sum along arbitrary dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data=calculate_sum(data,dim_name,varargin)

dim_num=data.dim_to_remove(find(strcmp(data.dim_names_to_remove,dim_name),1,'first'));
data.native.value=squeeze(sum(data.native.value,dim_num));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate integral in pressure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data=calculate_integral(data,varargin)

data_out=zeros(size(data.native.value,1),size(data.native.value,2),size(data.native.value,4));

do_trop=0;

ps=varargin{1};
a=varargin{2};
b=varargin{3};

pres=calculate_pressure(ps,a,b);

%Tropospheric integral?
if nargin>5
   if strcmp('tropopause',varargin{4})
      do_trop=1;
      trop=data.native_secondary.value;
   end
end

%Loop over all columns and integrate
for i=1:size(ps,1)
   disp([sprintf('%3.0f',100*i/size(ps,1)),'% complete'])
   for j=1:size(ps,2)
      for t=1:size(data.native.value,4)
         if do_trop==1
            [val,trop_ind]=min(abs(squeeze(trop(i,j,t)) - squeeze(pres(i,j,t,:))));
            if pres(i,j,t,trop_ind)>trop(i,j,t)
               trop_ind=trop_ind-1;
            end
            data_out(i,j,t)=trapz(squeeze(pres(i,j,trop_ind+1:end,t)),data.native.value(i,j,trop_ind+1:end,t));
            data_out(i,j,t)=data_out(i,j,t)+(pres(i,j,trop_ind+1,t)-trop(i,j,t))*data.native.value(i,j,trop_ind,t);
         else
            data_out(i,j,t)=trapz(squeeze(pres(i,j,:,t)),squeeze(data.native.value(i,j,:,t)));
         end
      end
   end
end
   
data.native.value=data_out/9.81; 

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate age of air
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tracer=age_of_air(tracer,p_ref,y_ref)

aoa=zeros(size(tracer.native.value));

%Setup indexing for eval
dim_string='(';
for i=1:length(size(tracer.native.value))
   dim_name=tracer.native_dim(i).name;
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
[val,y_ref]=min(abs(tracer.native_dim(lat_dim).value-y_ref));
[val,p_ref]=min(abs(tracer.native_dim(lev_dim).value-p_ref));
disp(['Reference latitude ',sprintf('%2.1f',tracer.native_dim(lat_dim).value(y_ref))])
disp(['Reference pressure ',sprintf('%2.1f',tracer.native_dim(lev_dim).value(p_ref))])
tracer.native.value=movmean(tracer.native.value,12,time_dim);
reference_timeseries=eval(['squeeze(mean(tracer.native.value',ref_string,',',num2str(lon_dim),'));']);
reference_timeseries=repmat(reference_timeseries(:),[1 length(reference_timeseries)]);

nt=size(tracer.native.value,time_dim);

%Calculate aoa using the matrix method
for y=1:length(tracer.native_dim(lat_dim).value)
   disp([sprintf('%2.2f',100*(y-1)/length(tracer.native_dim(lat_dim).value)),'% complete'])
   for x=1:length(tracer.native_dim(lon_dim).value)
      for p=1:length(tracer.native_dim(lev_dim).value)
         timeseries=eval(['squeeze(tracer.native.value',dim_string,')']);
         timeseries=repmat(timeseries(:),[1 length(timeseries)])';
         [val,time_index]=min(abs(timeseries-reference_timeseries),[],2);
         eval(['aoa',dim_string,'=(time_index''-[1:nt])/12;']);
      end
   end
end

tracer.native.value=aoa;

end
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

end
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

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create NetCDF output file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
   if length(var_out.dim{i}.out.value)>1
      var_dim_bnds(i)=netcdf.defVar(ncid,[var_out.dim{i}.out.info.out_name,'_bnds'],'NC_FLOAT',[dim(end) dim(i)]);
   end
   attnames=fieldnames(var_out.dim{i}.out.info);
   atts=struct2cell(var_out.dim{i}.out.info);
   for j=1:length(atts)
       if ~isempty(atts{j}) & sum(strcmp(attnames{j},{'out_name';'must_have_bounds';'stored_direction';'type';'requested'}))==0
          netcdf.putAtt(ncid,var_dim(i),attnames{j},atts{j});
          if length(var_out.dim{i}.out.value)>1
             netcdf.putAtt(ncid,var_dim_bnds(i),attnames{j},atts{j}); 
          end
       end
   end   
end

netcdf.endDef(ncid);

start=zeros(length(var_out.dim),1);
for i=1:length(var_out.dim)
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
   if length(var_out.dim{i}.out.value)>1
      netcdf.putVar(ncid,var_dim_bnds(i),var_out.dim{i}.out.bnds);
   end
end

netcdf.close(ncid);

end
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
      
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create bounding variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
      for j=2:size(bnds,2)-1
         bnds(1,j)=dim(j)-(dim(j)-dim(j-1))/2;
         bnds(2,j)=dim(j)+(dim(j+1)-dim(j))/2;
      end
      bnds(:,1)=[-(dim(2)-dim(1))/2 (dim(2)-dim(1))/2]+dim(1);
      bnds(:,end)=[-(dim(end)-dim(end-1))/2 (dim(end)-dim(end-1))/2]+dim(end);
   end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Interpolates along specified dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function field_out=interpolate_field(field,interp_dim,dim,missing_value,varargin)

do_native=0;
squeeze_field=0;

%Optional: load in surface pressure, a/b's for vertical interp, generate pressure array
if nargin>4
   a=varargin{1};
   b=varargin{2};
   ps=varargin{3};

   dim.native.value=calculate_pressure(ps,a,b);
   do_native=1;
end
%Set up output field
field_size=size(field);
field_size(interp_dim)=length(dim.interp);
field_out=zeros(field_size);

if do_native==1
   local_field_out=zeros(size(field_out,4),size(field_out,3));
end

%Interpolate
disp('Interpolating')
for i=1:size(field,1)
   disp([sprintf('%0.2f',100*i/size(field,1)),'% complete'])
   for j=1:size(field,2)
      if do_native==1
         local_field=squeeze(field(i,j,:,:))';
         local_dim=squeeze(dim.native.value(i,j,:,:))';
         for k=1:size(local_dim,1)
           local_field_out(k,:)=interp1(local_dim(k,:),local_field(k,:),...
                                  dim.interp,'linear',missing_value);
         end
         field_out(i,j,:,:)=local_field_out';
     else
         fieldsize=length(size(field_out));
         switch fieldsize
            case 2
               field_out(i,j)=interp1(dim.native.value,...
                                  squeeze(field(i,j,:)),dim.interp,...
                                  'linear',missing_value);
            case 3
               field_out(i,j,:)=interp1(dim.native.value,...
                                  squeeze(field(i,j,:)),dim.interp,...
                                  'linear',missing_value);
            case 4
               for k=1:size(field_out,3)
                  field_out(i,j,k,:)=interp1(dim.native.value,...
                                  squeeze(field(i,j,k,:)),dim.interp,...
                                  'linear',missing_value);
               end
         end
      end
   end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Collapse a field to a specified matrix rank
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Expand a field to an arbitrary matrix rank of 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Returns dimension information from specification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Returns the index of the matching field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ind=return_index(structure,field)
 
names=fieldnames(structure);
for i=1:length(names)
   if strcmp(field,names{i})
      ind=i;
   end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parses a string into individual words
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function parsed=parse_string(in)

parsed = strings(0);
while (in ~= "")
   [token,in] = strtok(in);
   parsed = [parsed ; token];
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Translate variable details, make actionable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [variable_out]=translate_cesm_variable(in)

opt_count=1;
in_components=fieldnames(in);
variables_components={'var';'operation';'eval';'p0';'lat0';'axis';'option'};
variable_out=struct('var',strings(0),'operation',strings(0),'eval',strings(0),'p0',[],'lat0',[],'axis',strings(0),'option',[]);

for i=1:length(in_components)
   for j=1:length(variables_components)-1
      if contains(in_components{i}, variables_components{j}) 
         eval(['variable_out.',variables_components{j},'=cat(1,variable_out.',variables_components{j},',in.',in_components{i},');']) 
      end
   end
   if contains(in_components{i}, 'option')
      variable_out.option(opt_count)=struct('var',in_components{i}.var,'string',in_components{i}.string);
      opt_count=opt_count+1;
   end 
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Scan dictionary and return field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=translate_cesm(dictionary,in,field,specification)

out=[];

dictionary_entries=fieldnames(eval(['dictionary.',field]));
for i=1:length(dictionary_entries)
   if strcmp(dictionary_entries{i},in)
      out=eval(['dictionary.',field,'.',dictionary_entries{i}]);
   end
end

%Non-exact match
if isempty(out)
   for i=1:length(dictionary_entries)
      if contains(in,dictionary_entries{i})
         out=eval(['dictionary.',field,'.',dictionary_entries{i}]);
      end
   end
end

if isempty(out)
   out='no_match';
end

end
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

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Remove entries from a list
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load file and decode json
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function decoded=json_load(file)
fid = fopen(file);
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
decoded=jsondecode(str);

end
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

end
