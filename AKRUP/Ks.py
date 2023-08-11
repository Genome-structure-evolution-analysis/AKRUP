import subprocess
from multiprocessing import cpu_count
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import as_completed

from Bio.Align.Applications import ClustalwCommandline
from Bio import SeqIO

import AKRUP
from AKRUP.funcbase import *

cpu_num = cpu_count()
best_thread_num = min(int(cpu_num/2), 5)
current_path = AKRUP.__path__[0]
platform_name = platform.platform().split('-')[0]


class CalculateExceptin(Exception):
    "This is to check for exceptions in calculating Ks"
    def __init__(self):
        pass


class RunKs:
    def __init__(self, options):

        self.species_cds1 = 'cds1 file'
        self.species_cds2 = 'cds2 file'
        self.block_file = 'block file'
        self.save_ks_file = 'save ks file  (*.ks.txt)'

        for k, v in options:
            setattr(self, str(k), v)

        self.ksscript = os.path.join(current_path, 'ini/calculate.Ks.pl')
        if platform_name == 'Linux':
            self.clustalw = os.path.join(current_path, 'ini/clustalw2')
        elif platform_name == 'Windows':
            self.clustalw = os.path.join(current_path, 'ini/clustalw2.exe')
        elif platform_name == 'macOS':
            pass

    @staticmethod
    def check_file(file_path):
        n = os.path.isfile(file_path)
        fn_name = os.path.basename(file_path)
        assert n, f"{fn_name}, File not exist!!!"   

    def get_seq(self):
        peps, cdss = {}, {}
        for cds in [self.species_cds1, self.species_cds2]:
            self.check_file(cds)
            for seq_record in SeqIO.parse(cds, 'fasta'):
                seq_pep = seq_record.seq.translate()
                seq_pep = str(seq_pep).replace('.', '').replace('*', '')
                cdss[seq_record.id] = seq_record.seq
                peps[seq_record.id] = seq_pep

        return peps, cdss

    def get_block_pair(self, path):
        for li in open(path):
            li = li.strip()
            if not li: continue
            if 'overlap' in li or 'self' in li or \
                    'LOCALE' in li or 'MAXIMUM' in li or \
                    'p-value' in li or '>' in li or 'the' in li:
                continue
            else:
                lis = re.split('\\s', li)
                yield [lis[0], lis[2]]

    def run_ks(self, path, input_file, out_file, cdsfile, ks_sf, cond, id1, id2, ksrow=''):
        try:
            clustalw_cline = ClustalwCommandline(self.clustalw, infile=f'{path}/{input_file}', outfile=f'{path}/{out_file}', tree=False)

            clustalw_cline()

            output = subprocess.run(f'perl {self.ksscript} {path} {out_file} {cdsfile} {id1} {id2}', encoding='utf-8', shell=True, capture_output=True)
            ksrow = output.stdout
            ksrow = re.findall('\<(.*?)\\s+\>', ksrow.strip())
        except Exception as e:
            raise f'Error: {e}'
        except:
            ksrow = f'{id1}\t{id2}\t-2\t-2\n'
        else:
            if ksrow:
                ksrow = ksrow[0]
            else:
                ksrow = f'{id1}\t{id2}\t-2\t-2\n'
        finally:

            saverow = '\t'.join(re.split('\\s', ksrow.strip()))
            if cond.acquire():
                ks_sf.write(saverow+'\n')
                ks_sf.flush()
                cond.release()

            os.remove(f'{path}/{input_file}')
            os.remove(f'{path}/{out_file}')
            os.remove(f'{path}/{cdsfile}')
            dnd_file = out_file.replace('aln', 'dnd')
            os.remove(f'{path}/{dnd_file}')
            return output.returncode, output.stderr

    def run(self):
        flag = True
        try:
            middle_path = 'run_middle_ks'+get_temdir()
            if not os.path.isdir(middle_path):
                os.mkdir(middle_path)
            cond = threading.Condition()
            thread = []
            # thead_pool = ThreadPoolExecutor(best_thread_num)
            with ThreadPoolExecutor(best_thread_num) as thead_pool:
                ks_sf = open(self.save_ks_file, 'w')

                peps, cdss = self.get_seq()
                pairs = self.get_block_pair(self.block_file)
                num = 0
                for pair in pairs:
                    if pair[0] == pair[1]: continue
                    pair_pep = f'{num}.pep'
                    pair_cds = f'{num}.cds'
                    pair_aln = f'{num}.aln'
                    try:
                        savepep = open(f'{middle_path}/{pair_pep}', 'w')
                        savecds = open(f'{middle_path}/{pair_cds}', 'w')
                        if pair[0] not in peps or pair[1] not in peps or\
                            pair[0] not in cdss or pair[1] not in cdss:
                            continue
                        savepep.write(f'>{pair[0]}\n{peps[pair[0]]}\n')
                        savepep.write(f'>{pair[1]}\n{peps[pair[1]]}\n')
                        savecds.write(f'>{pair[0]}\n{cdss[pair[0]]}\n')
                        savecds.write(f'>{pair[1]}\n{cdss[pair[1]]}\n')
                    finally:
                        savepep.close()
                        savecds.close()

                    handle = thead_pool.submit(self.run_ks, middle_path, pair_pep, pair_aln, pair_cds, ks_sf, cond, pair[0], pair[1])
                    thread.append(handle)
                    num += 1

                for th in as_completed(thread):
                    code, stderr = th.result()
                    if code != 0:
                        print(stderr)
                        flag = False
                        raise CalculateExceptin

        except CalculateExceptin:
            print("Calculate Ks Error: Check the dependent environment!!!!")
        else:
            if flag:
                print('Calculation Complete!!!')
        finally:
            ks_sf.close()
            shutil.rmtree(f'{middle_path}')
