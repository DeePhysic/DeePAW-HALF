#include <hdf5.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int half_is_hdf5_c(const char *filename) { return H5Fis_hdf5(filename) > 0 ? 1 : 0; }

int half_extract_potcar_c(const char *filename, char *output, int output_size) {
  hid_t file=-1,set=-1,type=-1;char *content=NULL,*variable=NULL;size_t length=0,written=0;int fd=-1,status=-1;
  file=H5Fopen(filename,H5F_ACC_RDONLY,H5P_DEFAULT);if(file<0||H5Lexists(file,"input/potcar/content",H5P_DEFAULT)<=0)goto done;
  set=H5Dopen2(file,"input/potcar/content",H5P_DEFAULT);type=H5Dget_type(set);if(set<0||type<0)goto done;
  if(H5Tis_variable_str(type)>0){if(H5Dread(set,type,H5S_ALL,H5S_ALL,H5P_DEFAULT,&variable)<0||!variable)goto done;length=strlen(variable);content=variable;}
  else{length=H5Tget_size(type);content=(char*)malloc(length+1);if(!content||H5Dread(set,type,H5S_ALL,H5S_ALL,H5P_DEFAULT,content)<0)goto done;content[length]='\0';}
  char pattern[]="/tmp/half-potcar-XXXXXX";fd=mkstemp(pattern);if(fd<0||(int)strlen(pattern)>=output_size)goto done;
  while(written<length){ssize_t n=write(fd,content+written,length-written);if(n<=0)goto done;written+=(size_t)n;}
  strcpy(output,pattern);status=0;
done:if(fd>=0)close(fd);if(variable)H5free_memory(variable);else free(content);if(type>=0)H5Tclose(type);if(set>=0)H5Dclose(set);if(file>=0)H5Fclose(file);return status;
}

static hid_t structure_group(hid_t file) {
  const char *paths[] = {"structure/positions", "input/poscar", "results/positions"};
  for (int i=0;i<3;i++) if (H5Lexists(file,paths[i],H5P_DEFAULT)>0) return H5Gopen2(file,paths[i],H5P_DEFAULT);
  return -1;
}

int half_probe_vaspwave_c(const char *filename, int *ntypes, int *nions, int grid[3]) {
  hid_t file=-1, group=-1, set=-1, space=-1; int *counts=NULL, status=-1; hsize_t dims[4];
  file=H5Fopen(filename,H5F_ACC_RDONLY,H5P_DEFAULT);if(file<0)goto done;group=structure_group(file);if(group<0)goto done;
  set=H5Dopen2(group,"number_ion_types",H5P_DEFAULT);if(set<0)goto done;space=H5Dget_space(set);if(space<0||H5Sget_simple_extent_dims(space,dims,NULL)<0)goto done;
  *ntypes=(int)dims[0];counts=(int*)malloc((size_t)*ntypes*sizeof(int));if(!counts)goto done;
  if(H5Dread(set,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,counts)<0)goto done;*nions=0;for(int i=0;i<*ntypes;i++)*nions+=counts[i];
  H5Sclose(space);space=-1;H5Dclose(set);set=H5Dopen2(file,"charge/grid",H5P_DEFAULT);if(set<0)goto done;
  if(H5Dread(set,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,grid)<0)goto done;status=0;
done: free(counts);if(space>=0)H5Sclose(space);if(set>=0)H5Dclose(set);if(group>=0)H5Gclose(group);if(file>=0)H5Fclose(file);return status;
}

int half_read_vaspwave_c(const char *filename, int ntypes, int nions, const int grid[3],
    char *system, int system_size, char *species, int *counts, double *lattice,
    double *positions, int *direct, double *charge) {
  hid_t file=-1,group=-1,set=-1,type=-1,memtype=-1,space=-1,memspace=-1;int status=-1;double scale=1.0;hsize_t dims[4],start[4]={0,0,0,0},count[4];
  file=H5Fopen(filename,H5F_ACC_RDONLY,H5P_DEFAULT);if(file<0)goto done;group=structure_group(file);if(group<0)goto done;
  memset(system,0,(size_t)system_size);strncpy(system,"VASP HDF5 structure",(size_t)system_size-1);
  if(H5Lexists(group,"system",H5P_DEFAULT)>0){set=H5Dopen2(group,"system",H5P_DEFAULT);type=H5Dget_type(set);size_t n=H5Tget_size(type);char *tmp=(char*)calloc(n+1,1);if(!tmp||H5Dread(set,type,H5S_ALL,H5S_ALL,H5P_DEFAULT,tmp)<0){free(tmp);goto done;}strncpy(system,tmp,(size_t)system_size-1);free(tmp);H5Tclose(type);type=-1;H5Dclose(set);set=-1;}
  if(H5Lexists(group,"scale",H5P_DEFAULT)>0){set=H5Dopen2(group,"scale",H5P_DEFAULT);if(H5Dread(set,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,&scale)<0)goto done;H5Dclose(set);set=-1;}
  set=H5Dopen2(group,"lattice_vectors",H5P_DEFAULT);if(set<0||H5Dread(set,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,lattice)<0)goto done;for(int i=0;i<9;i++)lattice[i]*=scale;H5Dclose(set);set=-1;
  set=H5Dopen2(group,"position_ions",H5P_DEFAULT);if(set<0||H5Dread(set,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,positions)<0)goto done;H5Dclose(set);set=-1;
  set=H5Dopen2(group,"number_ion_types",H5P_DEFAULT);if(set<0||H5Dread(set,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,counts)<0)goto done;H5Dclose(set);set=-1;
  memset(species,' ',(size_t)16*ntypes);set=H5Dopen2(group,"ion_types",H5P_DEFAULT);if(set<0)goto done;memtype=H5Tcopy(H5T_C_S1);if(memtype<0||H5Tset_size(memtype,16)<0||H5Dread(set,memtype,H5S_ALL,H5S_ALL,H5P_DEFAULT,species)<0)goto done;H5Tclose(memtype);memtype=-1;H5Dclose(set);set=-1;
  set=H5Dopen2(group,"direct_coordinates",H5P_DEFAULT);if(set<0||H5Dread(set,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,direct)<0)goto done;H5Dclose(set);set=-1;
  H5Gclose(group);group=-1;set=H5Dopen2(file,"charge/charge",H5P_DEFAULT);if(set<0)goto done;space=H5Dget_space(set);int rank=H5Sget_simple_extent_ndims(space);if(rank==4){H5Sget_simple_extent_dims(space,dims,NULL);count[0]=1;count[1]=dims[1];count[2]=dims[2];count[3]=dims[3];if(H5Sselect_hyperslab(space,H5S_SELECT_SET,start,NULL,count,NULL)<0)goto done;memspace=H5Screate_simple(4,count,NULL);if(H5Dread(set,H5T_NATIVE_DOUBLE,memspace,space,H5P_DEFAULT,charge)<0)goto done;}else if(rank==3){if(H5Dread(set,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,charge)<0)goto done;}else goto done;
  status=0;
done:if(memspace>=0)H5Sclose(memspace);if(space>=0)H5Sclose(space);if(memtype>=0)H5Tclose(memtype);if(type>=0)H5Tclose(type);if(set>=0)H5Dclose(set);if(group>=0)H5Gclose(group);if(file>=0)H5Fclose(file);return status;
}

static int scalar(hid_t parent, const char *name, hid_t type, const void *value) {
  hid_t space = H5Screate(H5S_SCALAR);
  hid_t set = space < 0 ? -1 : H5Dcreate2(parent, name, type, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  int status = (set < 0 || H5Dwrite(set, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, value) < 0) ? -1 : 0;
  if (set >= 0) H5Dclose(set); if (space >= 0) H5Sclose(space); return status;
}

static int array(hid_t parent, const char *name, hid_t type, int rank, const hsize_t *dims, const void *value) {
  hid_t space = H5Screate_simple(rank, dims, NULL);
  hid_t set = space < 0 ? -1 : H5Dcreate2(parent, name, type, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  int status = (set < 0 || H5Dwrite(set, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, value) < 0) ? -1 : 0;
  if (set >= 0) H5Dclose(set); if (space >= 0) H5Sclose(space); return status;
}

static int string_scalar(hid_t parent, const char *name, const char *value) {
  hid_t type = H5Tcopy(H5T_C_S1); size_t length = strlen(value); if (!length) length = 1;
  if (type < 0 || H5Tset_size(type, length) < 0) return -1;
  int status = scalar(parent, name, type, value); H5Tclose(type); return status;
}

int half_write_vaspwave_h5_c(const char *filename, const char *system,
    int ntypes, const char *species, const int *counts, int nions,
    const double *lattice, const double *positions, const int *grid,
    const double *charge, double encut, double fermi, int nk, int nb,
    const double *kpoints, const double *eigenvalues, const double *occupations,
    const int *npws, const int64_t *offsets, const float *coefficients) {
  hid_t file=-1, version=-1, structure=-1, pos=-1, charge_group=-1, wave=-1, spin=-1, point=-1, stype=-1;
  int status=-1, major=6, minor=6, patch=0, direct=1;
  double one=1.0, max_npw=0.0, rispin=1.0, rnb=(double)nb, rnk=(double)nk;
  hsize_t d1[1], d2[2], d3[3], d4[4];
  double *celtot=NULL, *fertot=NULL;
  file=H5Fcreate(filename,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT); if(file<0)goto done;
  version=H5Gcreate2(file,"version",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT); if(version<0)goto done;
  if(scalar(version,"major",H5T_NATIVE_INT,&major)||scalar(version,"minor",H5T_NATIVE_INT,&minor)||scalar(version,"patch",H5T_NATIVE_INT,&patch))goto done;
  H5Gclose(version);version=-1;
  structure=H5Gcreate2(file,"structure",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);if(structure<0)goto done;
  pos=H5Gcreate2(structure,"positions",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);if(pos<0)goto done;
  if(string_scalar(pos,"system",system)||scalar(pos,"scale",H5T_NATIVE_DOUBLE,&one))goto done;
  d2[0]=3;d2[1]=3;if(array(pos,"lattice_vectors",H5T_NATIVE_DOUBLE,2,d2,lattice))goto done;
  stype=H5Tcopy(H5T_C_S1);if(stype<0||H5Tset_size(stype,16)<0)goto done;d1[0]=(hsize_t)ntypes;
  if(array(pos,"ion_types",stype,1,d1,species))goto done;H5Tclose(stype);stype=-1;
  if(array(pos,"number_ion_types",H5T_NATIVE_INT,1,d1,counts))goto done;
  d2[0]=(hsize_t)nions;d2[1]=3;if(array(pos,"position_ions",H5T_NATIVE_DOUBLE,2,d2,positions))goto done;
  if(scalar(pos,"direct_coordinates",H5T_NATIVE_INT,&direct))goto done;
  H5Gclose(pos);pos=-1;H5Gclose(structure);structure=-1;
  charge_group=H5Gcreate2(file,"charge",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);if(charge_group<0)goto done;
  d1[0]=3;if(array(charge_group,"grid",H5T_NATIVE_INT,1,d1,grid))goto done;
  d4[0]=1;d4[1]=(hsize_t)grid[2];d4[2]=(hsize_t)grid[1];d4[3]=(hsize_t)grid[0];
  if(array(charge_group,"charge",H5T_NATIVE_DOUBLE,4,d4,charge))goto done;H5Gclose(charge_group);charge_group=-1;
  wave=H5Gcreate2(file,"wave",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);if(wave<0)goto done;
  d2[0]=3;d2[1]=3;if(array(wave,"amat",H5T_NATIVE_DOUBLE,2,d2,lattice))goto done;
  for(int k=0;k<nk;k++)if(npws[k]>max_npw)max_npw=npws[k];max_npw*=8.0;
  if(scalar(wave,"efermi",H5T_NATIVE_DOUBLE,&fermi)||scalar(wave,"enmax",H5T_NATIVE_DOUBLE,&encut)||
     scalar(wave,"rdum",H5T_NATIVE_DOUBLE,&max_npw)||scalar(wave,"rispin",H5T_NATIVE_DOUBLE,&rispin)||
     scalar(wave,"rnb_tot",H5T_NATIVE_DOUBLE,&rnb)||scalar(wave,"rnkpts",H5T_NATIVE_DOUBLE,&rnk))goto done;
  spin=H5Gcreate2(wave,"spin_1",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);if(spin<0)goto done;
  celtot=(double*)calloc((size_t)nb*2,sizeof(double));fertot=(double*)malloc((size_t)nb*sizeof(double));if(!celtot||!fertot)goto done;
  for(int k=0;k<nk;k++){
    char name[64];snprintf(name,sizeof(name),"kpoint_%d",k+1);point=H5Gcreate2(spin,name,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);if(point<0)goto done;
    memset(celtot,0,(size_t)nb*2*sizeof(double));for(int b=0;b<nb;b++){celtot[2*b]=eigenvalues[k*nb+b];fertot[b]=0.5*occupations[k*nb+b];}
    d2[0]=(hsize_t)nb;d2[1]=2;if(array(point,"celtot",H5T_NATIVE_DOUBLE,2,d2,celtot))goto done;
    d1[0]=(hsize_t)nb;if(array(point,"fertot",H5T_NATIVE_DOUBLE,1,d1,fertot))goto done;
    if(scalar(point,"num_planewaves",H5T_NATIVE_INT,&npws[k]))goto done;d1[0]=3;
    if(array(point,"vkpt",H5T_NATIVE_DOUBLE,1,d1,kpoints+3*k))goto done;
    d3[0]=(hsize_t)nb;d3[1]=(hsize_t)npws[k];d3[2]=2;
    if(array(point,"wave",H5T_NATIVE_FLOAT,3,d3,coefficients+offsets[k]))goto done;
    H5Gclose(point);point=-1;
  }
  status=0;
done:
  free(celtot);free(fertot);if(stype>=0)H5Tclose(stype);if(point>=0)H5Gclose(point);if(spin>=0)H5Gclose(spin);
  if(wave>=0)H5Gclose(wave);if(charge_group>=0)H5Gclose(charge_group);if(pos>=0)H5Gclose(pos);
  if(structure>=0)H5Gclose(structure);if(version>=0)H5Gclose(version);if(file>=0)H5Fclose(file);return status;
}
