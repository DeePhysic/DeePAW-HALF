#include <complex>
#include <iostream>
#include <vector>

#include <pybind11/embed.h>
#include <pybind11/numpy.h>

#include "interface_cpp.h"

namespace py = pybind11;
using namespace py::literals;

namespace vasp
{

typedef std::complex<double> complex_;
typedef py::array_t<int, py::array::f_style> int_array;
typedef py::array_t<double, py::array::f_style> double_array;
typedef py::array_t<complex_> complex_array;

void print_flags(py::array array) {
    py::object flags = array.attr("flags");
	std::cout << std::boolalpha;
    std::cout << "C_CONTIGUOUS : " << flags.attr("c_contiguous").cast<bool>() << std::endl;
    std::cout << "F_CONTIGUOUS : " << flags.attr("f_contiguous").cast<bool>() << std::endl;
    std::cout << "OWNDATA      : " << flags.attr("owndata").cast<bool>() << std::endl;
    std::cout << "WRITEABLE    : " << flags.attr("writeable").cast<bool>() << std::endl;
    std::cout << "ALIGNED      : " << flags.attr("aligned").cast<bool>() << std::endl;
    std::cout << "WRITEBACKIFCOPY : " << flags.attr("writebackifcopy").cast<bool>() << std::endl;
}

template <class Constants>
const py::dict convert_electrons(const Constants &constants)
{
    int size_grid = std::accumulate(
        constants.shape_grid, constants.shape_grid + 3, 1, std::multiplies<int>());
    double_array charge_density(size_grid, constants.charge_density);
    charge_density.resize(constants.shape_grid);
    return py::dict(
        "ENCUT"_a = constants.ENCUT,
        "NELECT"_a = constants.NELECT,
        "ZVAL"_a = double_array(constants.number_ion_types, constants.ZVAL),
        "shape_grid"_a = int_array(3, constants.shape_grid),
        "charge_density"_a = charge_density);
}

template <class Constants>
const py::dict convert_charge_density(const Constants &constants)
{
    int size_grid = std::accumulate(
        constants.shape_grid, constants.shape_grid + 3, 1, std::multiplies<int>());
    double_array charge_density(size_grid, constants.charge_density);
    charge_density.resize(constants.shape_grid);
    return py::dict(
        "shape_grid"_a = int_array(3, constants.shape_grid),
        "charge_density"_a = charge_density);
}

template <class Constants>
const py::dict convert_potentials(const Constants &constants)
{
    int size_grid = std::accumulate(
        constants.shape_grid, constants.shape_grid + 3, 1, std::multiplies<int>());
	double_array hartree_potential(size_grid, constants.hartree_potential);
	hartree_potential.resize(constants.shape_grid);
	double_array ion_potential(size_grid, constants.ion_potential);
	ion_potential.resize(constants.shape_grid);
	return py::dict("hartree_potential"_a = hartree_potential,
					"ion_potential"_a = ion_potential);
}

template <class Constants>
const py::dict convert_dipole(const Constants &constants)
{
	if (constants.dipole_moment_size > 0) {
		double_array dipole_moment(constants.dipole_moment_size, constants.dipole_moment);
		py::dict output("dipole_moment"_a = dipole_moment);
		return output;
	}
	else
	{
		py::dict output;
		return output;
	}
}


template <class Constants>
const py::dict convert_dynamics(const Constants &constants)
{
    int size_forces = constants.number_ions * 3;
    int size_stress = 9;
    double_array forces(size_forces, constants.forces);
    forces.resize({constants.number_ions, 3});
    double_array stress(size_stress, constants.stress);
    stress.resize({3, 3});
    return py::dict(
        "POMASS"_a = double_array(constants.number_ion_types, constants.POMASS),
        "forces"_a = forces,
        "stress"_a = stress);
}

template <class ConstantsOrAdditions>
const py::dict convert_occupancies(const ConstantsOrAdditions &constants_or_additions)
{
	return py::dict(
			"NELECT"_a = constants_or_additions.NELECT,
			"EFERMI"_a = constants_or_additions.EFERMI,
			"SIGMA"_a = constants_or_additions.SIGMA,
			"ISMEAR"_a = constants_or_additions.ISMEAR,
			"EMIN"_a = constants_or_additions.EMIN,
			"EMAX"_a = constants_or_additions.EMAX,
			"NUPDOWN"_a = constants_or_additions.NUPDOWN);
}

template <class Constants>
const py::dict convert_structure(const Constants &constants)
{
    int size_lattice_vectors = 9;
    int size_positions = constants.number_ions * 3;
    double_array lattice_vectors(size_lattice_vectors, constants.lattice_vectors);
    lattice_vectors.resize({3, 3});
    double_array positions(size_positions, constants.positions);
    positions.resize({constants.number_ions, 3});
    return py::dict(
        "number_ions"_a = constants.number_ions,
        "number_ion_types"_a = constants.number_ion_types,
        "lattice_vectors"_a = lattice_vectors,
        "ion_types"_a = int_array(constants.number_ions, constants.ion_types),
        "atomic_numbers"_a = int_array(constants.number_ion_types, constants.atomic_numbers),
        "positions"_a = positions);
}

template <class Constants>
const py::dict convert_neighbor_list(const py::object &module,
                                     const Constants &constants)
{
    if (!constants.neighbor_list)
    {
        return py::dict();
    }
    py::list neighbor_list;
    py::object Neighbors = module.attr("Neighbors");
    for (int ion = 0; ion < constants.number_ions; ion++)
    {
        auto list_ion = constants.neighbor_list[ion];
        int size_directions = list_ion.size * 3;
        double_array directions(size_directions, list_ion.directions);
        directions.resize({list_ion.size, 3});
        py::object neighbors = Neighbors(
            "neighbors"_a = int_array(list_ion.size, list_ion.neighbors),
            "distances"_a = double_array(list_ion.size, list_ion.distances),
            "directions"_a = directions);
        neighbor_list.append(neighbors);
    }
    return py::dict("neighbor_list"_a = neighbor_list);
}

const py::object mutable_double_array(const std::vector<int>& shape, double *data)
{
    int size = std::accumulate(shape.begin(), shape.end(), 1, std::multiplies<int>());
    py::capsule capsule(data, [](void *) {});
    double_array array(size, data, capsule);
    array.resize(shape);
    return array;
}

int force_and_stress(const py::object &module,
                     const constants_force_and_stress &constants,
                     additions_force_and_stress &additions)
{
#ifndef NO_FORCE_AND_STRESS
    py::object ConstantsForceAndStress = module.attr("ConstantsForceAndStress")(
        **convert_electrons(constants),
        **convert_dynamics(constants),
        **convert_structure(constants),
        **convert_neighbor_list(module, constants));
    py::object AdditionForceAndStress = module.attr("AdditionsForceAndStress")(
        "total_energy"_a = additions.total_energy,
        "forces"_a = mutable_double_array({constants.number_ions, 3}, additions.forces),
        "stress"_a = mutable_double_array({3, 3}, additions.stress));
    py::int_ error;
    if (constants.mode == FORCE_AND_STRESS)
    {
        error = module.attr("interface_force_and_stress")(
            ConstantsForceAndStress, AdditionForceAndStress);
    }
#ifndef NO_MACHINE_LEARNING
    else if (constants.mode == MACHINE_LEARNING)
    {
        error = module.attr("interface_machine_learning")(
            ConstantsForceAndStress, AdditionForceAndStress);
    }
#endif
    py::float_ total_energy = AdditionForceAndStress.attr("total_energy");
    additions.total_energy = total_energy.cast<float>();
    return error.cast<int>();
#else
    return NO_ERROR;
#endif
}

int local_potential(const py::object &module,
                    const constants_local_potential &constants,
                    additions_local_potential &additions)
{
#ifndef NO_LOCAL_POTENTIAL
    py::object ConstantsLocalPotential = module.attr("ConstantsLocalPotential")(
        **convert_electrons(constants),
        **convert_structure(constants),
		**convert_potentials(constants),
		**convert_dipole(constants),
        **convert_neighbor_list(module, constants));
    std::vector<int> shape(constants.shape_grid, constants.shape_grid + 3);
    py::object AdditionLocalPotential = module.attr("AdditionsLocalPotential")(
        "total_energy"_a = additions.total_energy,
        "total_potential"_a = mutable_double_array(shape, additions.total_potential));
    py::int_ error = module.attr("interface_local_potential")(
        ConstantsLocalPotential, AdditionLocalPotential);
    py::float_ total_energy = AdditionLocalPotential.attr("total_energy");
    additions.total_energy = total_energy.cast<float>();
    return error.cast<int>();
#else
    return NO_ERROR;
#endif
}

int occupancies(const py::object &module,
				const constants_occupancies &constants,
				additions_occupancies &additions)
{
#ifndef NO_OCCUPANCIES
	py::object ConstantsOccupancies = module.attr("ConstantsOccupancies")(
			**convert_occupancies(constants));
	py::object AdditionsOccupancies = module.attr("AdditionsOccupancies")(
			**convert_occupancies(additions));
    py::int_ error = module.attr("interface_occupancies")(ConstantsOccupancies, AdditionsOccupancies);
    additions.NELECT = AdditionsOccupancies.attr("NELECT").cast<float>();
    additions.EFERMI = AdditionsOccupancies.attr("EFERMI").cast<float>();
    additions.SIGMA = AdditionsOccupancies.attr("SIGMA").cast<float>();
    additions.ISMEAR = AdditionsOccupancies.attr("ISMEAR").cast<int>();
    additions.EMIN = AdditionsOccupancies.attr("EMIN").cast<float>();
    additions.EMAX = AdditionsOccupancies.attr("EMAX").cast<float>();
    additions.NUPDOWN = AdditionsOccupancies.attr("NUPDOWN").cast<float>();
    return error.cast<int>();
#else
    return NO_ERROR;
#endif
}

int update_structure(const py::object &module,
                     const constants_structure &constants,
                     additions_structure &additions)
{
#ifndef NO_STRUCTURE
    py::dict total_energy_constant = py::dict (
		"total_energy"_a = constants.total_energy);
    py::object ConstantsStructure = module.attr("ConstantsStructure")(
        **convert_dynamics(constants),
        **convert_structure(constants),
		**convert_charge_density(constants),
		**total_energy_constant,
        **convert_neighbor_list(module, constants));
    py::object AdditionStructure = module.attr("AdditionsStructure")(
        "lattice_vectors"_a = mutable_double_array({3, 3}, additions.lattice_vectors),
        "positions"_a = mutable_double_array({constants.number_ions, 3}, additions.positions));
    py::int_ error = module.attr("interface_structure")(ConstantsStructure, AdditionStructure);
    return error.cast<int>();
#else
    return NO_ERROR;
#endif
}


template <class Function, class Constants, class Additions>
int try_python_call(Function *function,
                    const Constants &constants,
                    Additions &additions)
{
    try
    {
        py::object module = py::module::import("vasp");
        return function(module, constants, additions);
    }
    catch (py::error_already_set &e)
    {
        if (e.matches(PyExc_ImportError))
        {
            return IMPORT_ERROR;
        }
        else
        {
            std::cerr << e.what() << std::endl;
            return BUG_IN_CPP;
        }
    };
}

} // namespace vasp

extern "C"
{
    int plugins_initialize()
    {
        py::initialize_interpreter();
        return NO_ERROR;
    }

    int plugins_finalize()
    {
        py::finalize_interpreter();
        return NO_ERROR;
    }

    int interface_force_and_stress(const constants_force_and_stress &constants,
                                   additions_force_and_stress &additions)
    {
        return vasp::try_python_call(&vasp::force_and_stress, constants, additions);
    }

    int interface_local_potential(const constants_local_potential &constants,
                                  additions_local_potential &additions)
    {
        return vasp::try_python_call(&vasp::local_potential, constants, additions);
    }

    int interface_update_structure(const constants_structure &constants,
                                   additions_structure &additions)
    {
        return vasp::try_python_call(&vasp::update_structure, constants, additions);
    }

	int interface_occupancies(const constants_occupancies &constants,
						      additions_occupancies &additions)
	{
		return vasp::try_python_call(&vasp::occupancies, constants, additions);
	}
}
