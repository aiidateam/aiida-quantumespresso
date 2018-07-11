# -*- coding: utf-8 -*-
import click
from aiida.utils.cli import command
from aiida.utils.cli import options


@command()
@options.code(callback_kwargs={'entry_point': 'quantumespresso.pw'})
@options.structure()
#@options.pseudo_family()
@options.daemon()
@click.option(
    '-z', '--protocol', type=click.Choice(['standard']), default='standard', show_default=True,
    help='the protocol to use for the workflow'
)
def launch(
    code, structure, pseudo_family, daemon, protocol):
    """
    Run the PwBandStructureWorkChain for a given input structure 
    to compute the band structure for the relaxed structure
    """
    from aiida.orm.data.base import Str
    from aiida.orm import DataFactory
    from aiida.orm.utils import WorkflowFactory
    from aiida.work.launch import run, submit

    PwBandStructureWorkChain = WorkflowFactory('quantumespresso.pw.band_structure')
    ParameterData = DataFactory('parameter')

    inputs = {
        'code': code,
        'structure': structure,
        ## Pseudo family is not anymore a parameter of this workflow. Instead,
        ## you should already have the pseudos in your DB, or pass a pseudo_data
        ## modifier below, with the MD5 of the pseudos you want to use.

        ## If you don't have the SSSP pseudopotentials
        ## (that you can download from here for SSSP v.1.0:
        ## https://www.materialscloud.org/archive/2018.0001/v2
        ##  
        ## you can get the dictionary with the md5 with the
        ## following code snippet, but you still need then to specify "cutoff" and "dual"
        ## for all relevant pseudos!
        ##
        ## CODE SNIPPET:
        ##
        ## def get_md5_dict(family_name):
        ##     UpfData = DataFactory('upf')
        ##     family = UpfData.get_upf_group(family_name)
        ##     return {node.element: {'md5': node.md5sum} for node in family.nodes}
        'protocol': ParameterData(dict={
            'name': 'theos-ht-1.0',
        }),
        ## or (to apply modifiers):
        # 'protocol': ParameterData(dict={
        #     'name: 'theos-ht-1.0',
        #     'modifiers': {
        #         'parameters': 'fast',
        #         'pseudo': 'SSSP-efficiency-1.0'
        #     }
        # })
        ## or (for custom-specified pseudos and cutoffs):
        # 'protocol': ParameterData(dict={
        #     'name: 'theos-ht-1.0',
        #     'modifiers': {
        #         'parameters': 'fast',
        #         'pseudo': 'custom',
        #         'pseudo_data': {
        #             "Ag": {
        #                 "cutoff": "50", 
        #                 "dual": "4", 
        #                 "filename": "Ag_ONCV_PBE-1.0.upf", 
        #                 "md5": "96da9acec54ba82f98e06bace1ebc465", 
        #                 "pseudopotential": "SG15"
        #             }, 
        #             "Al": {
        #                 "cutoff": "30", 
        #                 "dual": "8", 
        #                 "filename": "Al.pbe-n-kjpaw_psl.1.0.0.UPF", 
        #                 "md5": "4d58055f5a69695be6f94701d50bfe3f", 
        #                 "pseudopotential": "100PAW"
        #             }, 
        #             # ...
        #     }
        # })
    }

    if daemon:
        workchain = submit(PwBandStructureWorkChain, **inputs)
        click.echo('Submitted {}<{}> to the daemon'.format(PwBandStructureWorkChain.__name__, workchain.pk))
    else:
        run(PwBandStructureWorkChain, **inputs)
