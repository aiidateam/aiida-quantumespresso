"""Tests for the `PwCalculationTools` class."""

import numpy as np
import pytest
from aiida import orm
from aiida.common.links import LinkType
from pathlib import Path

from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData


def test_pw_get_scf_accuracy(fixture_localhost, generate_calc_job_node):
    """Test the `PwCalculationTools.get_scf_accuracy` method."""
    entry_point_name = 'quantumespresso.pw'

    # Missing `output_trajectory` node
    node = generate_calc_job_node(entry_point_name, fixture_localhost)
    with pytest.raises(ValueError):
        node.tools.get_scf_accuracy()

    # Missing `scf_accuracy` array
    node = generate_calc_job_node(entry_point_name, fixture_localhost)
    trajectory = orm.ArrayData()
    trajectory.base.links.add_incoming(node, link_type=LinkType.CREATE, link_label='output_trajectory')
    trajectory.store()

    with pytest.raises(ValueError):
        node.tools.get_scf_accuracy()

    # Missing `scf_accuracy_index` array
    node = generate_calc_job_node(entry_point_name, fixture_localhost)
    trajectory = orm.ArrayData()
    trajectory.set_array('scf_accuracy', np.array([1, 1, 1, 2, 2, 2, 2, 2]))
    trajectory.base.links.add_incoming(node, link_type=LinkType.CREATE, link_label='output_trajectory')
    trajectory.store()

    with pytest.raises(ValueError):
        node.tools.get_scf_accuracy()

    node = generate_calc_job_node(entry_point_name, fixture_localhost)
    trajectory = orm.ArrayData()
    trajectory.set_array('scf_accuracy', np.array([1, 1, 1, 2, 2, 2, 2, 2]))
    trajectory.set_array('scf_iterations', np.array([3, 5]))
    trajectory.base.links.add_incoming(node, link_type=LinkType.CREATE, link_label='output_trajectory')
    trajectory.store()

    # Invalid indices, there are only two frames
    with pytest.raises(IndexError):
        node.tools.get_scf_accuracy(index=2)
    with pytest.raises(IndexError):
        node.tools.get_scf_accuracy(index=-3)

    assert np.array_equal(node.tools.get_scf_accuracy(index=0), np.array([1, 1, 1]))
    assert np.array_equal(node.tools.get_scf_accuracy(index=1), np.array([2, 2, 2, 2, 2]))
    assert np.array_equal(node.tools.get_scf_accuracy(index=-1), np.array([2, 2, 2, 2, 2]))
    assert np.array_equal(node.tools.get_scf_accuracy(index=-2), np.array([1, 1, 1]))




def test_pw_get_occupations_not_hubbard(fixture_localhost, generate_calc_job_node):
    """get_occupations should raise ValueError when the node has no structure input (not a HubbardStructureData)."""
    entry_point_name = 'quantumespresso.pw'
    node = generate_calc_job_node(entry_point_name, fixture_localhost)

    # When inputs.structure does not exist the implementation should raise ValueError
    with pytest.raises(ValueError):
        node.tools.get_occupations()



def test_get_occupations(fixture_localhost, generate_calc_job_node):
    """
    Test get_occupations by uploading the local XML file into a temporary retrieved node.
    """
    entry_point_name = 'quantumespresso.pw'
    
    # --- 1. Define Expected Output ---
    expected = [
        {
            'atom_label': 1,
            'atom_specie': 'Fe1',
            'shell': '3d',
            'occupation_matrix': {
                'up': np.array([ 1.43157031e-01, -3.18183278e-03, -3.18183278e-03, -9.86022584e-17,
                                 -6.36366556e-03, -3.18183278e-03,  9.89707556e-01, -2.24844802e-04,
                                 -5.51109604e-03,  2.24844802e-04, -3.18183278e-03, -2.24844802e-04,
                                  9.89707556e-01,  5.51109604e-03,  2.24844802e-04, -1.00036285e-16,
                                 -5.51109604e-03,  5.51109604e-03,  1.43157031e-01, -6.41163792e-15,
                                 -6.36366556e-03,  2.24844802e-04,  2.24844802e-04, -6.41172177e-15,
                                  9.89707556e-01]),
                'down': np.array([ 1.54622465e-01,  2.17031883e-03,  2.17031883e-03, -1.32897775e-16,
                                   4.34063766e-03,  2.17031883e-03,  9.90426589e-01, -1.31551804e-04,
                                   3.75910248e-03,  1.31551804e-04,  2.17031883e-03, -1.31551804e-04,
                                   9.90426589e-01, -3.75910248e-03,  1.31551804e-04, -1.40663373e-16,
                                   3.75910248e-03, -3.75910248e-03,  1.54622465e-01, -6.41762390e-15,
                                   4.34063766e-03,  1.31551804e-04,  1.31551804e-04, -6.41731727e-15,
                                   9.90426589e-01])
            }
        },
        {
            'atom_label': 2,
            'atom_specie': 'Fe2',
            'shell': '3d',
            'occupation_matrix': {
                'up': np.array([ 2.56052690e-01, -1.25945459e-01, -1.25945459e-01,  6.20205994e-16,
                                 -2.51890918e-01, -1.25945459e-01,  9.06416480e-01,  4.25677469e-02,
                                 -2.18143934e-01, -4.25677469e-02, -1.25945459e-01,  4.25677469e-02,
                                  9.06416480e-01,  2.18143934e-01, -4.25677469e-02,  6.23690688e-16,
                                 -2.18143934e-01,  2.18143934e-01,  2.56052690e-01, -5.98747982e-15,
                                 -2.51890918e-01, -4.25677469e-02, -4.25677469e-02, -5.98658387e-15,
                                  9.06416480e-01]),
                'down': np.array([ 4.34421823e-01,  1.75448148e-01,  1.75448148e-01, -1.46695264e-15,
                                   3.50896296e-01,  1.75448148e-01,  7.68564262e-01,  1.10597256e-01,
                                   3.03885106e-01, -1.10597256e-01,  1.75448148e-01,  1.10597256e-01,
                                   7.68564262e-01, -3.03885106e-01, -1.10597256e-01, -1.47976561e-15,
                                   3.03885106e-01, -3.03885106e-01,  4.34421823e-01, -5.58088706e-15,
                                   3.50896296e-01, -1.10597256e-01, -1.10597256e-01, -5.58370980e-15,
                                   7.68564262e-01])
            }
        },
        {
            'atom_label': 3,
            'atom_specie': 'O1',
            'shell': '2p',
            'occupation_matrix': {
                'up': np.array([ 0.79784573,  0.00425208,  0.00425208,  0.        ,  0.        ,
                                 0.00425208,  0.79784573, -0.00425208,  0.        ,  0.        ,
                                 0.00425208, -0.00425208,  0.79784573,  0.        ,  0.        ,
                                 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                 0.        ,  0.        ,  0.        ,  0.        ,  0.        ]),
                'down': np.array([ 0.80131533, -0.00665188, -0.00665188,  0.        ,  0.        ,
                                  -0.00665188,  0.80131533,  0.00665188,  0.        ,  0.        ,
                                  -0.00665188,  0.00665188,  0.80131533,  0.        ,  0.        ,
                                   0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                   0.        ,  0.        ,  0.        ,  0.        ,  0.        ])
            }
        },
        {
            'atom_label': 4,
            'atom_specie': 'O1',
            'shell': '2p',
            'occupation_matrix': {
                'up': np.array([ 0.79784573,  0.00425208,  0.00425208,  0.        ,  0.        ,
                                 0.00425208,  0.79784573, -0.00425208,  0.        ,  0.        ,
                                 0.00425208, -0.00425208,  0.79784573,  0.        ,  0.        ,
                                 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                 0.        ,  0.        ,  0.        ,  0.        ,  0.        ]),
                'down': np.array([ 0.80131533, -0.00665188, -0.00665188,  0.        ,  0.        ,
                                  -0.00665188,  0.80131533,  0.00665188,  0.        ,  0.        ,
                                  -0.00665188,  0.00665188,  0.80131533,  0.        ,  0.        ,
                                   0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                   0.        ,  0.        ,  0.        ,  0.        ,  0.        ])
            }
        },
    ]

    # --- 2. Setup Structure  ---
    base_structure = orm.StructureData(cell=[[10, 0, 0], [0, 10, 0], [0, 0, 10]])
    base_structure.append_atom(position=(0, 0, 0), symbols='Fe')
    
    structure = HubbardStructureData.from_structure(base_structure)
    structure.store()

    # --- 3. Generate Node with Inputs ---
    inputs = {'structure': structure}
    node = generate_calc_job_node(entry_point_name, fixture_localhost, inputs=inputs)

    # --- 4. Upload Local XML to the Node ---
    local_xml_path = Path(__file__).parent / 'data-file-schema.xml'
    
    if not local_xml_path.exists():
        pytest.fail(f"Could not find 'data-file-schema.xml' at {local_xml_path}")

    # Create the Retrieved folder
    retrieved = orm.FolderData()
    retrieved.put_object_from_file(local_xml_path, 'data-file-schema.xml')
    
    # Link retrieved to the calculation (Output links allow linking to stored nodes)
    retrieved.base.links.add_incoming(node, link_type=LinkType.CREATE, link_label='retrieved')
    retrieved.store()

    # --- 5. Run Test ---
    result = node.tools.get_occupations(reshape=False)

    assert result is not None
    assert_nested_almost_equal(result, expected)
    

# Helper for Numpy comparison 
def assert_nested_almost_equal(actual, expected):
    if isinstance(expected, dict):
        assert actual.keys() == expected.keys()
        for key in expected:
            assert_nested_almost_equal(actual[key], expected[key])
    elif isinstance(expected, list):
        assert len(actual) == len(expected)
        for i in range(len(expected)):
            assert_nested_almost_equal(actual[i], expected[i])
    elif isinstance(expected, np.ndarray):
        np.testing.assert_allclose(actual, expected, atol=1e-8)
    else:
        assert actual == expected