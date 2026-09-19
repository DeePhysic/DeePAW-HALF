from dataclasses import dataclass

import ase
import numpy as np
from ase.calculators.calculator import Calculator
from ase.units import J, Pascal, m
from vasp import AdditionsForceAndStress, ConstantsForceAndStress


@dataclass
class AseForceField:
    calculator: Calculator

    def force_and_stress(
        self,
        constants: ConstantsForceAndStress,
        additions: AdditionsForceAndStress,
    ):
        atoms = self._atoms_from_constants(constants)
        total_energy = atoms.get_potential_energy()
        forces = atoms.get_forces()
        stress = self._stress_in_vasp_units(atoms)
        # convert results to double precision numpy arrays
        additions.total_energy += np.array(total_energy, dtype=np.float64)
        additions.forces += np.array(forces, dtype=np.float64)
        additions.stress += np.array(stress, dtype=np.float64)

    def _atoms_from_constants(self, constants):
        atoms = ase.Atoms(
            cell=constants.lattice_vectors, scaled_positions=constants.positions
        )
        atoms.set_atomic_numbers(constants.atomic_numbers[constants.ion_types])
        atoms.set_pbc(True)
        atoms.calc = self.calculator
        return atoms

    def _stress_in_vasp_units(self, atoms):
        try:
            stress = atoms.get_stress(voigt=False)
            return (-1 * J * atoms.get_volume() / (m**3 * Pascal)) * stress
        except PropertyNotImplementedError:
            return np.zeros((3, 3))
