"""
Check that the CP calculation might not work.
"""
from aiida.plugins.factories import CalculationFactory
from aiida.orm import Code, StructureData, Dict, KpointsData
from aiida.orm.nodes.data.upf import get_pseudos_from_structure
from aiida.engine import run
from ase.lattice.spacegroup import crystal


def main():
    alat = 5.4
    ase_structure = crystal(
        "Si",
        [(0, 0, 0)],
        spacegroup=227,
        cellpar=[alat, alat, alat, 90, 90, 90],
        primitive_cell=True,
    )
    structure = StructureData(ase=ase_structure)
    structure.store()

    calculation = CalculationFactory("quantumespresso.cp")
    print(calculation.__class__.__name__)
    builder = calculation.get_builder()

    # get some code
    code = Code.objects.all()[0]
    builder.code = code

    # Kind of have to import this one first
    structure = StructureData.objects.all()[0]
    builder.structure = structure

    # Parameters that we of course don't know about
    parameters = Dict(
        dict={
            "CONTROL": {
                "calculation": "cp",
                "restart_mode": "from_scratch",
                "wf_collect": False,
                "iprint": 1,
                "isave": 100,
                "dt": 3.0,
                "max_seconds": 25 * 60,
                "nstep": 10,
            },
            "SYSTEM": {
                "ecutwfc": 30.0,
                "ecutrho": 240.0,
                "nr1b": 24,
                "nr2b": 24,
                "nr3b": 24,
            },
            "ELECTRONS": {
                "electron_damping": 1.0e-1,
                "electron_dynamics": "damp",
                "emass": 400.0,
                "emass_cutoff": 3.0,
            },
            "IONS": {"ion_dynamics": "none"},
        }
    )
    builder.parameters = parameters

    pseudos = get_pseudos_from_structure(structure, "SSSP-Efficiency")
    builder.pseudos = pseudos

    kpoints = KpointsData()
    kpoints.set_kpoints_mesh([2, 2, 2])
    kpoints.store()
    # builder.kpoints = kpoints

    builder.metadata.dry_run = True
    builder.metadata.store_provenance = True
    builder.metadata.options = {"resources": {"num_machines": 1}}
    results, node = run.get_node(builder)
    print(results, node, sep="\n")


if __name__ == "__main__":
    main()
