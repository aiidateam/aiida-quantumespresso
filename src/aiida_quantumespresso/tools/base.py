# -*- coding: utf-8 -*-
"""Base class for parser of Quantum ESPRESSO pw.x and cp.x input files based on generic parser of `qe-tools` package."""
import re

from aiida.common import InputValidationError


class StructureParseMixin:
    """Mixin that extends ``~qe_tools.parsers.qeinputparser.QeInputFile`` to parse a ``StructureData``."""

    # pylint: disable=too-few-public-methods

    def get_structuredata(self):
        """Return a StructureData object based on the data in the input file.

        All of the names corresponding of the ``Kind`` objects composing the ``StructureData`` object will match those
        found in the ``ATOMIC_SPECIES`` block, so the pseudo potentials can be linked to the calculation using the kind
        name for each specific type of atom (in the event that you wish to use different pseudo's for two or more of the
        same atom).

        :return: structure data node of the structure defined in the input file.
        :rtype: :class:`~aiida.orm.nodes.data.structure.StructureData`
        """
        from aiida.orm.nodes.data.structure import Kind, Site, StructureData

        valid_elements_regex = re.compile(
            """
            (?P<ele>
    H  | He |
    Li | Be | B  | C  | N  | O  | F  | Ne |
    Na | Mg | Al | Si | P  | S  | Cl | Ar |
    K  | Ca | Sc | Ti | V  | Cr | Mn | Fe | Co | Ni | Cu | Zn | Ga | Ge | As | Se | Br | Kr |
    Rb | Sr | Y  | Zr | Nb | Mo | Tc | Ru | Rh | Pd | Ag | Cd | In | Sn | Sb | Te | I  | Xe |
    Cs | Ba | Hf | Ta | W  | Re | Os | Ir | Pt | Au | Hg | Tl | Pb | Bi | Po | At | Rn |
    Fr | Ra | Rf | Db | Sg | Bh | Hs | Mt |

    La | Ce | Pr | Nd | Pm | Sm | Eu | Gd | Tb | Dy | Ho | Er | Tm | Yb | Lu | # Lanthanides
    Ac | Th | Pa | U  | Np | Pu | Am | Cm | Bk | Cf | Es | Fm | Md | No | Lr | # Actinides
            )
            [^a-z]  # Any specification of an element is followed by some number
                    # or capital letter or special character.
        """, re.X | re.I
        )

        data = self.structure
        species = self.atomic_species

        structure = StructureData()
        structure.base.attributes.set('cell', data['cell'])

        for mass, name, pseudo in zip(species['masses'], species['names'], species['pseudo_file_names']):
            try:
                symbols = valid_elements_regex.search(pseudo).group('ele').capitalize()
            except Exception as exception:
                raise InputValidationError(
                    f'could not determine element name from pseudo name: {pseudo}'
                ) from exception
            structure.append_kind(Kind(name=name, symbols=symbols, mass=mass))

        for symbol, position in zip(data['atom_names'], data['positions']):
            structure.append_site(Site(kind_name=symbol, position=position))

        return structure
