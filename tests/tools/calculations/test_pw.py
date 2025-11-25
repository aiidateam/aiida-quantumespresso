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
    """get_occupations should raise ValueError when the node has a structure input which is not a HubbardStructureData."""
    entry_point_name = 'quantumespresso.pw'

    base_structure = orm.StructureData(cell=[[10, 0, 0], [0, 10, 0], [0, 0, 10]])
    base_structure.append_atom(position=(0, 0, 0), symbols='Fe')
    base_structure.store()
    

    node = generate_calc_job_node(entry_point_name, fixture_localhost, inputs={'structure': base_structure})

    with pytest.raises(ValueError):
        node.tools.get_occupations()

@pytest.mark.parametrize(
    'xml_filename',
    [
        'unpolarized-data-file-schema.xml',
        'collinear-data-file-schema.xml',
        'noncollinear-data-file-schema.xml',
    ]
)
def test_get_occupations(fixture_localhost, generate_calc_job_node, data_regression, num_regression, xml_filename, filepath_tests):
    """
    Test get_occupations by uploading local XML files into a temporary retrieved node.
    Tests unpolarized, collinear polarized, and non-collinear calculations.
    """
    entry_point_name = 'quantumespresso.pw'
    
    # --- 1. Setup Structure  ---
    base_structure = orm.StructureData(cell=[[10, 0, 0], [0, 10, 0], [0, 0, 10]])
    base_structure.append_atom(position=(0, 0, 0), symbols='Fe')
    
    structure = HubbardStructureData.from_structure(base_structure)
    structure.store()

    # --- 2. Generate Node with Inputs ---
    inputs = {'structure': structure}
    node = generate_calc_job_node(entry_point_name, fixture_localhost, inputs=inputs)

    # --- 3. Upload Local XML to the Node ---
    
    local_xml_path = Path(filepath_tests) / 'tools' / 'calculations' /  xml_filename 
    
    if not local_xml_path.exists():
        pytest.fail(f"Could not find '{xml_filename}' at {local_xml_path}")

    # Create the Retrieved folder
    retrieved = orm.FolderData()
    retrieved.put_object_from_file(local_xml_path, 'data-file-schema.xml')
    
    # Link retrieved to the calculation
    retrieved.base.links.add_incoming(node, link_type=LinkType.CREATE, link_label='retrieved')
    retrieved.store()

    # --- 4. Run Test ---
    result = node.tools.get_occupations(reshape=False)

    # Prepare data for regressions
    metadata = []
    numeric_data = {}
    
    for i, item in enumerate(result):
        metadata.append({
            'atom_index': item['atom_index'],
            'kind_name': item['kind_name'],
            'manifold': item['manifold'],
        })
        
        # Handle different occupation_matrix structures
        occ_matrix = item['occupations']
        for key, val in occ_matrix.items():
            numeric_data[f'atom_{i}_{key}'] = val
    
    # Check with pytest-regressions
    data_regression.check({'occupations_metadata': metadata})
    num_regression.check(numeric_data, default_tolerance={'atol': 1e-8, 'rtol': 0})