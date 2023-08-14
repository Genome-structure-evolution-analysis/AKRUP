import os
import click

import AKRUP
from AKRUP.load_blockinfo import BlockInfo
from AKRUP.dotplot_trajectory import DotPlotTrajectory
from AKRUP.event_Correspondence import PolyploidyEvent
from AKRUP.inferred_ancestral_karyotype import ConservedAncestralRegions
from AKRUP.dotplot_CSR_event import DotplotCsrKs
from AKRUP.blast_block_dotplot import DotPlot
from AKRUP.geteventblock import EventBlock
from AKRUP.ancestral_genome_seqs import ancestral_karyotype
from AKRUP.ColinearScan import RunColinearScan
from AKRUP.Blast import RunBlast
from AKRUP.ks_fig import KsDistrubute
from AKRUP.Ks import RunKs
from AKRUP.ancestral_karyotype import KaryotypeFig


from AKRUP.funcbase import load_conf

path = AKRUP.__path__[0]

sections = {'runblast': 'blast',
            'runcolinearscan': 'colinearscan',
           'runks': 'ks',
           'loadblock': 'loadblock',
           'dotplot': 'dotplot',
           'blockdotplot': 'blockdotplot',
           'eventblock': 'eventblock',
           'eventdotplot': 'eventdotplot',
           'ksfigure': 'ksdistribute',
           'csrdotplot': 'csrdotplot',
           'event_correspondence': 'Polyploidy_CSR',
           'trajectorydotplot': 'karyotype',
           'inferranckaryotype': 'ancestral',
           'anckaryotypefig': 'ancestralfig',
           'ancgenomeseqs': 'ancestralseq',
           }

optionfiles = {
            'rb': 'run_blast.conf',
            'rc': 'run_ColinearScan.conf',
            'rk': 'run_ks.conf',
            'lk': 'Loadblock.conf',
            'd': 'blast_dotplot.conf',
            'bd': 'block_dotplot.conf',
            'eb': 'event_block.conf',
            'ed': 'event_dotplot.conf',
            'kf': 'ksdistribute.conf',
            'cd': 'CSR_dotplot.conf',
            'ec': 'Polyploidy_CSR.conf',
            'td': 'dotplot_trajectory.conf',
            'iak': 'infer_anckaryotype.conf',
            'akf': 'ancestral_plotfig.conf',
            'ags': 'ancestralseq.conf',
           }

def runmodule(argument):
    switch = {
        'runblast': RunBlast,
        'runcolinearscan': RunColinearScan,
        'runks': RunKs,
        'loadblock': BlockInfo,
        'dotplot': DotPlot,
        'blockdotplot': DotPlot,
        'eventblock': EventBlock,
        'eventdotplot': DotplotCsrKs,
        'csrdotplot': DotplotCsrKs,
        'ksfigure': KsDistrubute,
        'event_correspondence': PolyploidyEvent,
        'trajectorydotplot': DotPlotTrajectory,
        'inferranckaryotype': ConservedAncestralRegions,
        'anckaryotypefig': KaryotypeFig,
        'ancgenomeseqs': ancestral_karyotype,
    }
    return switch.get(argument)

@click.command()
@click.help_option("-h", "--help", help='Show this help message and exit')
@click.version_option(version='1.0.4',message='version %(version)s', package_name='AKRUP')
@click.option('-rb', '--runblast', help='Search for potential homologous gene pairs',type=click.Path(exists=True))
@click.option('-rc', '--runcolinearscan', help='Infer genomic collinearity information',type=click.Path(exists=True))
@click.option('-rk', '--runks', help='Calculate Ka/Ks for homologous gene pairs',type=click.Path(exists=True))
@click.option('-d', '--dotplot', help='Show homologous gene dotplot',type=click.Path(exists=True))
@click.option('-bd', '--blockdotplot', help='Show synteny block dotplot',type=click.Path(exists=True))
@click.option('-eb', '--eventblock', help='Obtain event-related syntenic region',type=click.Path(exists=True))
@click.option('-ed', '--eventdotplot', help='Show event-related syntenic region dotplot',type=click.Path(exists=True))
@click.option('-kf', '--ksfigure', help='Draw Ks distribution',type=click.Path(exists=True))
@click.option('-lk', '--loadblock', help='Load collinearity information',type=click.Path(exists=True))
@click.option('-ec', '--event-correspondence', help='Extract event-related syntenic region',type=click.Path(exists=True))
@click.option('-cd', '--csrdotplot', help='Show continuous syntenic regions dotplot',type=click.Path(exists=True))
@click.option('-iak', '--inferranckaryotype', help='Inferring ancestral karyotypes',type=click.Path(exists=True))
@click.option('-akf', '--anckaryotypefig', help='Draw karyotypes figure',type=click.Path(exists=True))
@click.option('-ags', '--ancgenomeseqs', help='Extraction of ancestral genome sequence',type=click.Path(exists=True))
@click.option('-td', '--trajectorydotplot', help='Show ancestal karyotype trajectory',type=click.Path(exists=True))
@click.option('-e', '--example', help='Displays the configured parameters', 
  type=click.Choice(['rb', 'rc', 'rk','d','bd', 'eb', 'ed', 'kf', 'lk', 'ec', 'cd', 'iak', 'akf', 'ags', 'td']))
def main(**args):
    '''
    AKRUP: Ancestral Karyotype Reconstruction Universal Pipeline
    
    \b
    Include "bottom-up" inferences of ancestral karyotypes and
    "top-down" inferences of ancient chromosome evolutionary trajectories
    '''
    try:
        for k, v in args.items():
            if k == 'example' and v != None:
                return click.echo(open(os.path.join(path, f'example/{optionfiles.get(v)}')).read()+'\n')
            elif v != None:
                with open(v, 'r', encoding='utf-8-sig') as hanf:
                    click.echo(hanf.read()+'\n')
                options = load_conf(v, sections.get(k))

                if k in ['dotplot', 'blockdotplot', 'csrdotplot', 'eventdotplot']:
                    p = runmodule(k)(options, k)
                    p.run()
                else:
                    p = runmodule(k)(options)
                    p.run()

    except PermissionError:
        click.echo(click.style('Permission denied', fg='red'))
    except Exception as e:
        click.echo(click.style(f'Error: {e}', fg='red'))
    except:
        click.echo(click.style('Please check the configuration file or parameters', fg='red'))

if __name__=='__main__':
    main()