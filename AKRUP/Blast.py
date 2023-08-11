import subprocess
from multiprocessing import cpu_count

import AKRUP
from AKRUP.funcbase import *

cpu_num = cpu_count()
current_path = AKRUP.__path__[0]
platform_name = platform.platform().split('-')[0]

class RunBlast:
    def __init__(self, options):

        self.num_thread = '8' # auto
        self.evalue = '1e-5'
        self.outfmt = '6'
        self.max_target_seqs = '10'
        self.querypep = 'query pep file'
        self.subjectpep = 'subject pep file'
        self.outblast = 'save blast file (spec_spec.blast)'

        for k, v in options:
            setattr(self, str(k), v)

        if self.num_thread == 'auto':
            self.num_thread = int(cpu_num/2)
        if platform_name == 'Windows':
            
            self.makeblastdb = os.path.join(current_path, 'ini/makeblastdb.exe')
            self.blastp = os.path.join(current_path, 'ini/blastp.exe')

        elif platform_name == 'Linux':
            self.makeblastdb = os.path.join(current_path, 'ini/makeblastdb')
            self.blastp = os.path.join(current_path, 'ini/blastp')

    def run(self):
        path = get_temdir()
        if not os.path.isdir(path):
            os.mkdir(path)

        try:
            subprocess.call(f'{self.makeblastdb} -in {self.subjectpep} -dbtype prot -out {path}/bt.db', shell=True)

            p = subprocess.call(f'{self.blastp} -query {self.querypep} -out {self.outblast} -db {path}/bt.db'
                                f' -outfmt {self.outfmt} -evalue {self.evalue} -max_target_seqs {self.max_target_seqs} -num_threads '
                                f'{self.num_thread}', shell=True)
            if p != 0:
                print(f'Warning!!! {p}')
        except Exception as e:
            print(f'Error: {e}')
        finally:
            shutil.rmtree(f'{path}')
