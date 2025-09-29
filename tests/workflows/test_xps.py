"""Tests for the `XpsWorkChain` class."""

from aiida import orm
from aiida.common import LinkType
from plumpy import ProcessState


def test_default(generate_workchain_xps):
    """Test instantiating the WorkChain, then mock its process, by calling methods in the ``spec.outline``."""

    wkchain = generate_workchain_xps()

    assert wkchain.setup() is None
    assert wkchain.should_run_relax() is False

    # generate structure for scf calculation
    wkchain.prepare_structures()
    # run all scf and return the node
    all_scf_inputs = wkchain.run_all_scf()
    for label, node in all_scf_inputs.items():
        # add result to the scf node, used for wkchain.results()
        # "energy" will be used in the get_spectra_by_element calcfunction
        result = orm.Dict({'energy': 0})
        result.store()
        result.base.links.add_incoming(node, link_type=LinkType.RETURN, link_label='output_parameters')
        # set state to finished.
        node.set_process_state(ProcessState.FINISHED)
        node.set_exit_status(0)
        wkchain.ctx[label] = node
    assert wkchain.inspect_all_scf() is None
    # process results
    wkchain.results()
    wkchain.update_outputs()
    assert set(wkchain.node.base.links.get_outgoing().all_link_labels()) == {
        'standardized_structure',
        'supercell_structure',
        'symmetry_analysis_data',
        'output_parameters_ch_scf__site_0',
        'chemical_shifts__Si_cls',
        'final_spectra_cls__Si_cls_spectra',
    }
