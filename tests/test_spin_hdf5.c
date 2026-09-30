#include <hdf5.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

int half_write_vaspwave_h5_c(const char *filename, const char *system,
    int ntypes, const char *species, const int *counts, int nions,
    const double *lattice, const double *positions, const int *grid,
    const double *charge, double encut, double fermi, int nk, int nb,
    int nspin, const double *kpoints, const double *eigenvalues,
    const double *occupations, const int *npws, const int64_t *offsets,
    const float *coefficients);

static int read_double(hid_t file, const char *path, double *value) {
  hid_t set=H5Dopen2(file,path,H5P_DEFAULT);
  if(set<0)return -1;
  int status=H5Dread(set,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,value)<0?-1:0;
  H5Dclose(set);return status;
}

int main(void) {
  char path[]="/tmp/half-spin-hdf5-XXXXXX";
  int fd=mkstemp(path);if(fd<0)return 1;close(fd);unlink(path);
  const char species[16]={'X',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '};
  const int counts[1]={1},grid[3]={1,1,1},npws[1]={1};
  const int64_t offsets[1]={0};
  const double lattice[9]={1,0,0,0,1,0,0,0,1},positions[3]={0,0,0};
  const double charge[1]={1},kpoints[3]={0,0,0};
  const double eigenvalues[4]={-1,1,-0.8,1.2};
  const double occupations[4]={1,0,0.5,0};
  const float coefficients[4]={1,0,0,0};
  if(half_write_vaspwave_h5_c(path,"spin test",1,species,counts,1,lattice,positions,grid,
      charge,400,0,1,2,2,kpoints,eigenvalues,occupations,npws,offsets,coefficients))return 2;
  hid_t file=H5Fopen(path,H5F_ACC_RDONLY,H5P_DEFAULT);if(file<0)return 3;
  double rispin=0,up[2]={0,0},down[2]={0,0};int status=0;
  if(read_double(file,"wave/rispin",&rispin)||fabs(rispin-2)>1e-12)status=4;
  if(H5Lexists(file,"wave/spin_2/kpoint_1/wave",H5P_DEFAULT)<=0)status=5;
  if(read_double(file,"wave/spin_1/kpoint_1/fertot",up)||
     read_double(file,"wave/spin_2/kpoint_1/fertot",down))status=6;
  if(fabs(up[0]-1)>1e-12||fabs(up[1])>1e-12||fabs(down[0]-0.5)>1e-12||fabs(down[1])>1e-12)status=7;
  H5Fclose(file);unlink(path);return status;
}
