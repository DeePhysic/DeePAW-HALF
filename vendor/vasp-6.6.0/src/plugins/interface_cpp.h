#ifndef PLUGINS_INTERFACE_H_
#define PLUGINS_INTERFACE_H_

extern "C"
{
    enum Error
    {
        NO_ERROR,
        BUG_IN_CPP,
        IMPORT_ERROR,
        PYTHON_EXCEPTION // this should be 3 to match the Python code
    };

    enum Mode
    {
        FORCE_AND_STRESS,
        MACHINE_LEARNING,
    };

    struct neighbor_list_data
    {
        int size;
        int *neighbors;     // dim: size
        double *distances;  // dim: size
        double *directions; // dim: size, 3
    };

    struct constants_force_and_stress  // also used for machine learning
    {
        double ENCUT;
        double NELECT;
        int shape_grid[3];
        int number_ions;
        int number_ion_types;
        int *ion_types;          // dim: number_ions
        int *atomic_numbers;     // dim: number_ion_types
        double *lattice_vectors; // dim: number_ion_types
        double *positions;       // dim: number_ions, 3
        double *ZVAL;            // dim: number_ion_types
        double *POMASS;          // dim: number_ion_types
        double *forces;          // dim: number_ions, 3
        double *stress;          // dim: 3,3
        double *charge_density;  // dim: shape_grid[1,2,3]
        int mode;  // determines whether force_and_stress or machine_learning interface is called
        neighbor_list_data *neighbor_list; // dim: number_ions
    };

    struct additions_force_and_stress  // also used for machine learning
    {
        double total_energy;
        double *forces; // dim: number_ions, 3
        double *stress; // dim: 3,3
    };

    struct constants_local_potential
    {
        double ENCUT;
        double NELECT;
        int shape_grid[3];
        int number_ions;
        int number_ion_types;
        int *ion_types;            // dim: number_ions
        int *atomic_numbers;       // dim: number_ion_types
        double *lattice_vectors;   // dim: 3,3
        double *positions;         // dim: number_ions, 3
        double *ZVAL;              // dim: number_ion_types
        double *charge_density;    // dim: shape_grid[1,2,3]
		double *hartree_potential; // dim: shape_grid[1,2,3]
		double *ion_potential;     // dim: shape_grid[1,2,3]
		int dipole_moment_size;    // size of the dipole_moment array
		double *dipole_moment;     // dim: 3 or 0 depending on if computed by VASP
        neighbor_list_data *neighbor_list; // dim: number_ions
    };

    struct additions_local_potential
    {
        double total_energy;
        double *total_potential; // dim: shape_grid[1,2,3]
    };

	struct constants_occupancies
	{
		double NELECT;
		double EFERMI;
		double SIGMA;
		int ISMEAR;
		double EMIN;
		double EMAX;
		double NUPDOWN;
	};

	struct additions_occupancies
	{
		double NELECT;
		double EFERMI;
		double SIGMA;
		int ISMEAR;
		double EMIN;
		double EMAX;
		double NUPDOWN;
	};

    struct constants_structure
    {
        int number_ions;
        int number_ion_types;
        int *ion_types;          // dim: number_ions
        int *atomic_numbers;     // dim: number_ion_types
        double *lattice_vectors; // dim: 3,3
        double *positions;       // dim: number_ions, 3
        double *POMASS;          // dim: number_ion_types
		double total_energy;    // dim: 1
        double *forces;          // dim: number_ions, 3
        double *stress;          // dim: 3,3
        int shape_grid[3];
        double *charge_density;  // dim: shape_grid[1,2,3]
        neighbor_list_data *neighbor_list; // dim: number_ions
    };

    struct additions_structure
    {
        double *lattice_vectors; // dim: 3,3
        double *positions;       // dim: number_ions, 3
    };

    int plugins_initialize();
    int plugins_finalize();
    int interface_force_and_stress(const constants_force_and_stress &constants,
                                   additions_force_and_stress &additions);
    int interface_local_potential(const constants_local_potential &constants,
                                  additions_local_potential &additions);
    int interface_update_structure(const constants_structure &constants,
                                   additions_structure &additions);
    int interface_occupancies(const constants_occupancies &constants,
								   additions_occupancies &additions);
}

#endif // PLUGINS_INTERFACE_H_
