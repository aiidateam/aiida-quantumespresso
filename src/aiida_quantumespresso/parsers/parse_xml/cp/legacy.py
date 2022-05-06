# -*- coding: utf-8 -*-
"""Code that was written to parse the legacy XML format of Quantum ESPRESSO, which was deprecated in version 6.4."""
from xml.dom.minidom import parseString

from aiida_quantumespresso.parsers import QEOutputParsingError
from aiida_quantumespresso.parsers.parse_xml.legacy import (
    parse_xml_child_attribute_str,
    parse_xml_child_bool,
    parse_xml_child_float,
    parse_xml_child_integer,
    parse_xml_child_str,
    read_xml_card,
    xml_card_cell,
    xml_card_exchangecorrelation,
    xml_card_header,
    xml_card_ions,
    xml_card_planewaves,
    xml_card_spin,
)

units_suffix = '_units'
default_energy_units = 'eV'
default_k_points_units = '1 / angstrom'
default_length_units = 'Angstrom'


# TODO: the xml has a lot in common with pw, maybe I should avoid duplication of code
# or maybe should I wait for the new version of data-file.xml ?
def parse_cp_xml_output(data):
    """Parse xml data data must be a single string, as returned by file.read() (notice the difference with
    parse_text_output!) On output, a dictionary with parsed values.

    Democratically, we have decided to use picoseconds as units of time, eV for energies, Angstrom for lengths.
    """
    import copy

    dom = parseString(data)

    parsed_data = {}

    #CARD HEADER
    parsed_data = copy.deepcopy(xml_card_header(parsed_data, dom))

    # CARD CONTROL

    cardname = 'CONTROL'
    target_tags = read_xml_card(dom, cardname)

    tagname = 'PP_CHECK_FLAG'
    parsed_data[tagname.lower()] = parse_xml_child_bool(tagname, target_tags)

    # CARD STATUS

    cardname = 'STATUS'
    target_tags = read_xml_card(dom, cardname)

    tagname = 'STEP'
    attrname = 'ITERATION'
    parsed_data[(tagname + '_' + attrname).lower()] = int(parse_xml_child_attribute_str(tagname, attrname, target_tags))

    tagname = 'TIME'
    attrname = 'UNITS'
    value = parse_xml_child_float(tagname, target_tags)

    units = parse_xml_child_attribute_str(tagname, attrname, target_tags)
    if units not in ['pico-seconds']:
        raise QEOutputParsingError(f'Units {units} are not supported by parser')
    parsed_data[tagname.lower()] = value

    tagname = 'TITLE'
    parsed_data[tagname.lower()] = parse_xml_child_str(tagname, target_tags)

    # CARD CELL
    parsed_data, lattice_vectors, volume = copy.deepcopy(xml_card_cell(parsed_data, dom))

    # CARD IONS
    parsed_data = copy.deepcopy(xml_card_ions(parsed_data, dom, lattice_vectors, volume))

    # CARD PLANE WAVES

    parsed_data = copy.deepcopy(xml_card_planewaves(parsed_data, dom, 'cp'))

    # CARD SPIN
    parsed_data = copy.deepcopy(xml_card_spin(parsed_data, dom))

    # CARD EXCHANGE_CORRELATION
    parsed_data = copy.deepcopy(xml_card_exchangecorrelation(parsed_data, dom))

    # TODO CARD OCCUPATIONS

    # CARD BRILLOUIN ZONE
    # TODO: k points are saved for CP... Why?

    cardname = 'BRILLOUIN_ZONE'
    target_tags = read_xml_card(dom, cardname)

    tagname = 'NUMBER_OF_K-POINTS'
    parsed_data[tagname.replace('-', '_').lower()] = parse_xml_child_integer(tagname, target_tags)

    tagname = 'UNITS_FOR_K-POINTS'
    attrname = 'UNITS'
    metric = parse_xml_child_attribute_str(tagname, attrname, target_tags)
    if metric not in ['2 pi / a']:
        raise QEOutputParsingError(
            f'Error parsing attribute {attrname}, tag {tagname} inside {target_tags.tagName}, units unknown'
        )
    parsed_data[tagname.replace('-', '_').lower()] = metric

    # TODO: check what happens if one does not use the monkhorst pack in the code
    tagname = 'MONKHORST_PACK_GRID'
    try:
        a = target_tags.getElementsByTagName(tagname)[0]
        value = [int(a.getAttribute('nk' + str(i + 1))) for i in range(3)]
        parsed_data[tagname.replace('-', '_').lower()] = value
    except:
        raise QEOutputParsingError(f'Error parsing tag {tagname} inside {target_tags.tagName}.')

    tagname = 'MONKHORST_PACK_OFFSET'
    try:
        a = target_tags.getElementsByTagName(tagname)[0]
        value = [int(a.getAttribute('k' + str(i + 1))) for i in range(3)]
        parsed_data[tagname.replace('-', '_').lower()] = value
    except:
        raise QEOutputParsingError(f'Error parsing tag {tagname} inside {target_tags.tagName}.')

    try:
        kpoints = []
        for i in range(parsed_data['number_of_k_points']):
            tagname = 'K-POINT.' + str(i + 1)
            a = target_tags.getElementsByTagName(tagname)[0]
            b = a.getAttribute('XYZ').replace('\n', '').rsplit()
            value = [float(s) for s in b]

            metric = parsed_data['units_for_k_points']
            if metric == '2 pi / a':
                value = [float(s) / parsed_data['lattice_parameter'] for s in value]

                weight = float(a.getAttribute('WEIGHT'))

                kpoints.append([value, weight])

        parsed_data['k_point'] = kpoints
    except:
        raise QEOutputParsingError(f'Error parsing tag K-POINT.# inside {target_tags.tagName}.')

    tagname = 'NORM-OF-Q'
    # TODO decide if save this parameter
    parsed_data[tagname.replace('-', '_').lower()] = parse_xml_child_float(tagname, target_tags)

    # CARD PARALLELISM
    # can be optional

    try:
        cardname = 'PARALLELISM'
        target_tags = read_xml_card(dom, cardname)

        tagname = 'GRANULARITY_OF_K-POINTS_DISTRIBUTION'
        parsed_data[tagname.lower().replace('-', '_')] = parse_xml_child_integer(tagname, target_tags)

        tagname = 'NUMBER_OF_PROCESSORS'
        parsed_data[tagname.lower().replace('-', '_')] = parse_xml_child_integer(tagname, target_tags)

        tagname = 'NUMBER_OF_PROCESSORS_PER_POOL'
        parsed_data[tagname.lower().replace('-', '_')] = parse_xml_child_integer(tagname, target_tags)

        tagname = 'NUMBER_OF_PROCESSORS_PER_IMAGE'
        parsed_data[tagname.lower().replace('-', '_')] = parse_xml_child_integer(tagname, target_tags)

        tagname = 'NUMBER_OF_PROCESSORS_PER_TASKGROUP'
        parsed_data[tagname.lower().replace('-', '_')] = parse_xml_child_integer(tagname, target_tags)

        tagname = 'NUMBER_OF_PROCESSORS_PER_POT'
        parsed_data[tagname.lower().replace('-', '_')] = parse_xml_child_integer(tagname, target_tags)

        tagname = 'NUMBER_OF_PROCESSORS_PER_BAND_GROUP'
        parsed_data[tagname.lower().replace('-', '_')] = parse_xml_child_integer(tagname, target_tags)

        tagname = 'NUMBER_OF_PROCESSORS_PER_DIAGONALIZATION'
        parsed_data[tagname.lower().replace('-', '_')] = parse_xml_child_integer(tagname, target_tags)
    except QEOutputParsingError:
        pass

    # CARD TIMESTEPS

    cardname = 'TIMESTEPS'
    target_tags = read_xml_card(dom, cardname)

    for tagname in ['STEP0', 'STEPM']:
        try:
            tag = target_tags.getElementsByTagName(tagname)[0]

            try:
                second_tagname = 'ACCUMULATORS'
                second_tag = tag.getElementsByTagName(second_tagname)[0]
                data = second_tag.childNodes[0].data.rstrip().split()  # list of floats
                parsed_data[second_tagname.replace('-', '_').lower()] = [float(i) for i in data]
            except:
                pass

            second_tagname = 'IONS_POSITIONS'
            second_tag = tag.getElementsByTagName(second_tagname)[0]
            third_tagname = 'stau'
            third_tag = second_tag.getElementsByTagName(third_tagname)[0]
            list_data = third_tag.childNodes[0].data.rstrip().split()
            list_data = [float(i) for i in list_data]
            # convert to matrix
            val = []
            mat = []
            for i, data in enumerate(list_data):
                val.append(data)
                if (i + 1) % 3 == 0:
                    mat.append(val)
                    val = []
            parsed_data[(second_tagname + '_' + third_tagname).replace('-', '_').lower()] = mat
            third_tagname = 'svel'
            third_tag = second_tag.getElementsByTagName(third_tagname)[0]
            list_data = third_tag.childNodes[0].data.rstrip().split()
            list_data = [float(i) for i in list_data]
            # convert to matrix
            val = []
            mat = []
            for i, data in enumerate(list_data):
                val.append(data)
                if (i + 1) % 3 == 0:
                    mat.append(val)
                    val = []
            parsed_data[(second_tagname + '_' + third_tagname).replace('-', '_').lower()] = mat
            try:
                third_tagname = 'taui'
                third_tag = second_tag.getElementsByTagName(third_tagname)[0]
                list_data = third_tag.childNodes[0].data.rstrip().split()
                list_data = [float(i) for i in list_data]
                # convert to matrix
                val = []
                mat = []
                for i, data in enumerate(list_data):
                    val.append(data)
                    if (i + 1) % 3 == 0:
                        mat.append(val)
                        val = []
                parsed_data[(second_tagname + '_' + third_tagname).replace('-', '_').lower()] = mat
            except:
                pass

            try:
                third_tagname = 'cdmi'
                third_tag = second_tag.getElementsByTagName(third_tagname)[0]
                list_data = third_tag.childNodes[0].data.rstrip().split()
                parsed_data[(second_tagname + '_' + third_tagname).replace('-', '_').lower()
                            ] = [float(i) for i in list_data]
            except:
                pass

            try:
                third_tagname = 'force'
                third_tag = second_tag.getElementsByTagName(third_tagname)[0]
                list_data = third_tag.childNodes[0].data.rstrip().split()
                list_data = [float(i) for i in list_data]
                # convert to matrix
                val = []
                mat = []
                for i, data in enumerate(list_data):
                    val.append(data)
                    if (i + 1) % 3 == 0:
                        mat.append(val)
                        val = []
                parsed_data[(second_tagname + '_' + third_tagname).replace('-', '_').lower()] = mat
            except:
                pass

            second_tagname = 'IONS_NOSE'
            second_tag = tag.getElementsByTagName(second_tagname)[0]
            third_tagname = 'nhpcl'
            third_tag = second_tag.getElementsByTagName(third_tagname)[0]
            parsed_data[(second_tagname + '_' + third_tagname).replace('-', '_').lower()
                        ] = float(third_tag.childNodes[0].data)
            third_tagname = 'nhpdim'
            third_tag = second_tag.getElementsByTagName(third_tagname)[0]
            parsed_data[(second_tagname + '_' + third_tagname).replace('-', '_').lower()
                        ] = float(third_tag.childNodes[0].data)
            third_tagname = 'xnhp'
            third_tag = second_tag.getElementsByTagName(third_tagname)[0]
            parsed_data[(second_tagname + '_' + third_tagname).replace('-', '_').lower()
                        ] = float(third_tag.childNodes[0].data)
            try:
                third_tagname = 'vnhp'
                third_tag = second_tag.getElementsByTagName(third_tagname)[0]
                parsed_data[(second_tagname + '_' + third_tagname).replace('-', '_').lower()
                            ] = float(third_tag.childNodes[0].data)
            except:
                pass

            try:
                second_tagname = 'ekincm'
                second_tag = tag.getElementsByTagName(second_tagname)[0]
                parsed_data[second_tagname.replace('-', '_').lower()] = float(second_tag.childNodes[0].data)
            except:
                pass

            second_tagname = 'ELECTRONS_NOSE'
            second_tag = tag.getElementsByTagName(second_tagname)[0]
            try:
                third_tagname = 'xnhe'
                third_tag = second_tag.getElementsByTagName(third_tagname)[0]
                parsed_data[(second_tagname + '_' + third_tagname).replace('-', '_').lower()
                            ] = float(third_tag.childNodes[0].data)
            except:
                pass
            try:
                third_tagname = 'vnhe'
                third_tag = second_tag.getElementsByTagName(third_tagname)[0]
                parsed_data[(second_tagname + '_' + third_tagname).replace('-', '_').lower()
                            ] = float(third_tag.childNodes[0].data)
            except:
                pass

            second_tagname = 'CELL_PARAMETERS'
            second_tag = tag.getElementsByTagName(second_tagname)[0]
            try:
                third_tagname = 'ht'
                third_tag = second_tag.getElementsByTagName(third_tagname)[0]
                list_data = third_tag.childNodes[0].data.rstrip().split()
                list_data = [float(i) for i in list_data]
                # convert to matrix
                val = []
                mat = []
                for i, data in enumerate(list_data):
                    val.append(data)
                    if (i + 1) % 3 == 0:
                        mat.append(val)
                        val = []
                parsed_data[(second_tagname + '_' + third_tagname).replace('-', '_').lower()] = mat
            except:
                pass
            try:
                third_tagname = 'htvel'
                third_tag = second_tag.getElementsByTagName(third_tagname)[0]
                list_data = third_tag.childNodes[0].data.rstrip().split()
                list_data = [float(i) for i in list_data]
                # convert to matrix
                val = []
                mat = []
                for i, data in enumerate(list_data):
                    val.append(data)
                    if (i + 1) % 3 == 0:
                        mat.append(val)
                        val = []
                parsed_data[(second_tagname + '_' + third_tagname).replace('-', '_').lower()] = mat
            except:
                pass
            try:
                third_tagname = 'gvel'
                third_tag = second_tag.getElementsByTagName(third_tagname)[0]
                list_data = third_tag.childNodes[0].data.rstrip().split()
                list_data = [float(i) for i in list_data]
                # convert to matrix
                val = []
                mat = []
                for i, data in enumerate(list_data):
                    val.append(data)
                    if (i + 1) % 3 == 0:
                        mat.append(val)
                        val = []
                parsed_data[(second_tagname + '_' + third_tagname).replace('-', '_').lower()] = mat
            except:
                pass

            second_tagname = 'CELL_NOSE'
            second_tag = tag.getElementsByTagName(second_tagname)[0]
            try:
                third_tagname = 'xnhh'
                third_tag = second_tag.getElementsByTagName(third_tagname)[0]
                list_data = third_tag.childNodes[0].data.rstrip().split()
                list_data = [float(i) for i in list_data]

                # convert to matrix
                val = []
                mat = []
                for i, data in enumerate(list_data):
                    val.append(data)
                    if (i + 1) % 3 == 0:
                        mat.append(val)
                        val = []
                parsed_data[(second_tagname + '_' + third_tagname).replace('-', '_').lower()] = mat
            except:
                pass
            try:
                third_tagname = 'vnhh'
                third_tag = second_tag.getElementsByTagName(third_tagname)[0]
                list_data = third_tag.childNodes[0].data.rstrip().split()
                list_data = [float(i) for i in list_data]
                # convert to matrix
                val = []
                mat = []
                for i, data in enumerate(list_data):
                    val.append(data)
                    if (i + 1) % 3 == 0:
                        mat.append(val)
                        val = []
                parsed_data[(second_tagname + '_' + third_tagname).replace('-', '_').lower()] = mat
            except:
                pass
        except Exception as e:
            raise QEOutputParsingError(f'Error parsing CARD {cardname}')

    # CARD BAND_STRUCTURE_INFO

    cardname = 'BAND_STRUCTURE_INFO'
    target_tags = read_xml_card(dom, cardname)

    tagname = 'NUMBER_OF_ATOMIC_WFC'
    parsed_data[tagname.lower().replace('-', '_')] = parse_xml_child_integer(tagname, target_tags)

    tagname = 'NUMBER_OF_ELECTRONS'
    parsed_data[tagname.lower().replace('-', '_')] = int(parse_xml_child_float(tagname, target_tags))

    tagname = 'NUMBER_OF_BANDS'
    parsed_data[tagname.lower().replace('-', '_')] = parse_xml_child_integer(tagname, target_tags)

    tagname = 'NUMBER_OF_SPIN_COMPONENTS'
    parsed_data[tagname.lower().replace('-', '_')] = parse_xml_child_integer(tagname, target_tags)

    # TODO
    # - EIGENVALUES (that actually just contains occupations)
    #   Why should I be interested in that, if CP works for insulators only?
    # - EIGENVECTORS
    # - others TODO are written in the function

    return parsed_data
