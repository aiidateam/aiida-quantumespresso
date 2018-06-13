# -*- coding: utf-8 -*-
from aiida.cmdline.utils.decorators import with_dbenv
from aiida.common.exceptions import NotExistent


@with_dbenv()
def validate_hubbard_parameters(structure, parameters, hubbard_u=None, hubbard_v=None, hubbard_file_pk=None):
    """
    Validate Hubbard input parameters and update the parameters input node accordingly. If a valid hubbard_file_pk
    is provided, the node will be loaded and returned

    :param structure: the StructureData node that will be used in the inputs
    :param parameters: the ParameterData node that will be used in the inputs
    :param hubbard_u: the Hubbard U inputs values from the cli
    :param hubbard_v: the Hubbard V inputs values from the cli
    :param hubbard_file_pk: a pk referencing a SinglefileData with Hubbard parameters
    :returns: the loaded SinglefileData node with Hubbard parameters if valid pk was defined, None otherwise
    """
    from aiida.orm import load_node
    from aiida.orm.data.singlefile import SinglefileData

    if [v is None for v in [hubbard_u, hubbard_v, hubbard_file_pk]].count(True) > 1:
        raise ValueError('the hubbard_u, hubbard_v and hubbard_file_pk options are mutually exclusive')

    hubbard_file = None

    if hubbard_file_pk:

        try:
            hubbard_file = load_node(pk=hubbard_file_pk)
        except NotExistent:
            ValueError('{} is not a valid pk'.format(hubbard_file_pk))
        else:
            if not isinstance(hubbard_file, SinglefileData):
                ValueError('Node<{}> is not a SinglefileData but {}'.format())

        parameters['SYSTEM']['lda_plus_u'] = True
        parameters['SYSTEM']['lda_plus_u_kind'] = 2
        parameters['SYSTEM']['hubbard_parameters'] = 'file'

    elif hubbard_v:

        parameters['SYSTEM']['lda_plus_u'] = True
        parameters['SYSTEM']['lda_plus_u_kind'] = 2
        parameters['SYSTEM']['hubbard_parameters'] = 'input'
        parameters['SYSTEM']['hubbard_v'] = []

        for value in hubbard_v:
            parameters['SYSTEM']['hubbard_v'].append(value)

    elif hubbard_u:

        structure_kinds = structure.get_kind_names()
        hubbard_kinds = [value[0] for value in hubbard_u]

        if not set(hubbard_kinds).issubset(structure_kinds):
            raise ValueError('the kinds in the specified starting Hubbard values is not a strict subset of the kinds in the structure')

        parameters['SYSTEM']['lda_plus_u'] = True
        parameters['SYSTEM']['lda_plus_u_kind'] = 0
        parameters['SYSTEM']['hubbard_u'] = {}

        for kind, value in hubbard_u:
            parameters['SYSTEM']['hubbard_u'][kind] = value

    return hubbard_file
