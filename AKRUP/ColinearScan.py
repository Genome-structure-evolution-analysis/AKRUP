import subprocess
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import as_completed

import AKRUP
from AKRUP.funcbase import *

cpu_num = cpu_count()
current_path = AKRUP.__path__[0]

platform_name = platform.platform().split('-')[0]


class RunColinearScan:
    def __init__(self, options):
        # self.num_process = 'num_process'  # num/auto
        self.num_thread = 'num/auto'  # num/auto
        self.lens_file1 = 'lens1 file'
        self.lens_file2 = 'lens2 file'
        self.gff_file1 = 'gff1 file'
        self.gff_file2 = 'gff2 file'
        self.blast_file = 'blast file'
        self.save_block_file = 'save block file (*.block.rr.txt)'
        

        for k, v in options:
            setattr(self, str(k), v)

        if self.num_thread == 'auto':
            self.num_thread = int(cpu_num/2)

        if platform_name == 'Linux':
            self.blockscan = os.path.join(current_path, 'ini/blockscan')
        elif platform_name == 'Windows':
            self.blockscan = os.path.join(current_path, 'ini/blockscan.exe')
        elif platform_name == 'macOS':
            pass

    @staticmethod
    def check_file(file_path):
        n = os.path.isfile(file_path)
        file_name = os.path.basename(file_path)
        if not n:
            print(f'{file_name}---File not exist!!!')

    @staticmethod
    def check_del(file_path):
        n = os.path.isfile(file_path)
        if n:
            os.remove(file_path)

    @staticmethod
    def get_chr_num(gene):
        spec_na = re.findall('^\\D+', gene)[0]
        chr_num = re.search('\\d+', gene)[0]
        return [spec_na, int(chr_num)]

    @staticmethod
    def isself(pos):
        start1, start2, end1, end2 = pos
        if abs(min(start1, end1) - min(start2, end2)) < 30 and \
           abs(max(start2, end2) - max(start1, end1)) < 50:
            return 1
        else:
            return 0

    @staticmethod
    def isoverlap(pos):
        s1, e1 = min(pos[0], pos[2]), max(pos[0], pos[2])
        s2, e2 = min(pos[1], pos[3]), max(pos[1], pos[3])
        bs1, be1 = min(pos[4], pos[6]), max(pos[4], pos[6])
        bs2, be2 = min(pos[5], pos[7]), max(pos[5], pos[7])
        if s1 >= bs1-30 and e1 <= be1+30 and s2 >= bs2-50 and e2 <= be2+50:
            return 1
        else:
            return 0

    def get_lens(self, lens_path):
        lens_dit = {}
        self.check_file(lens_path)
        lens_f = open(lens_path, 'r')
        name = os.path.basename(lens_path).split('.')[0]
        for li in lens_f:
            lis = li.strip().split()
            lis[0] = str(int(lis[0]))
            lens_dit[lis[0]] = lis[1]
        chr_num = len(lens_dit.keys())

        return lens_dit, chr_num

    def get_gff(self, gff_path):
        self.check_file(gff_path)
        gff_dit = {}
        ch = {'+': '1', '-': '-1'}
        gff_f = open(gff_path, 'r')
        for li in gff_f:
            lis = li.strip().split()
            gff_dit[lis[5]] = ch[lis[3]]+' '+lis[6]
        return gff_dit

    def remove_redundancy(self, block_path, sf_path):
        block, old_flag, reason, isadd = [], [], '', 1
        block_f = open(block_path, 'r')
        sf = open(sf_path, 'a+')
        for li in block_f:
            if re.search('^\\s+', li):
                sf.write(li)
            elif any(([x in li for x in ['the', '+']])):
                sf.write(li)
            elif '>' in li:
                start = re.split('\\s+', block[0])
                end = re.split('\\s+', block[-1])
                start1, start2 = int(float(start[1])), int(float(start[3]))
                end1, end2 = int(float(end[1])), int(float(end[3]))

                for i in range(len(old_flag)):
                    pos = [start1, start2, end1, end2]
                    if self.get_chr_num(start[0])[0] == self.get_chr_num(start[2])[0]:
                        flag1 = self.isself(pos)
                        if flag1:
                            isadd, reason = 0, "self homologous"
                    lis = old_flag[i].split('\t')
                    pos.extend([int(x) for x in lis])
                    flag2 = self.isoverlap(pos)
                    if flag2:
                        isadd, reason = 0, f'overlap with block {i}+1th'
                        break
                if isadd:
                    old_flag.append(
                        '\t'.join([str(start1), str(start2), str(end1), str(end2)]))
                    for row in block:
                        sf.write(row)
                    sf.write(li)
                else:
                    sf.write(li)
                    sf.write(reason+'\n')
                    isadd = 1
                block = []
            else:
                block.append(li)

    def blast_get_pair(self, blast_path, path):
        if not os.path.isdir(path):
            os.mkdir(path)
        na1, na2 = os.path.basename(blast_path).split('.')[0].split('_')
        sf_path = f'{path}/{na1}_all_{na2}_all.pairs'
        sf = open(sf_path, 'w')
        try:
            self.check_file(blast_path)
            gff1 = self.get_gff(self.gff_file1)
            gff2 = self.get_gff(self.gff_file2)
            df = pd.read_csv(blast_path, sep='\t', names=list(
                range(1, 13)), index_col=False)
            df2 = df[(df[12].astype('float') >= 100) & (df[2] != df[1])]
            df3 = df2.loc[:, [1, 2]]
            for li in df3.values.tolist():
                flag = '\t'.join(li)
                if any(([x in flag for x in ['Un', 'Scaffold', 'random']])):
                    continue
                if na1 in li[0] and na2 in li[1]:
                    str1 = gff1.get(li[0], 0)
                    str2 = gff2.get(li[1], 0)
                    if not all((str1, str2)):
                        continue
                    sf.write(f'{li[0]} {str1} {li[1]} {str2}\n')
        finally:
            sf.close()

        try:
            purged_path = f'{path}/{na1}_all_{na2}_all.purged'
            purged_save = open(purged_path, 'w')

            pair_path = open(sf_path)
            pair_lines = pair_path.readlines()
            count_dit = defaultdict(int)
            for pair in pair_lines:
                pairs = re.split('\\s+', pair)
                count_dit[pairs[0]] += 1
                count_dit[pairs[3]] += 1
            for pair in pair_lines:
                pairs = re.split('\\s+', pair.strip())
                if pairs[0] == pairs[1]:continue
                if count_dit[pairs[0]] >= 50 or count_dit[pairs[3]] >=50:continue
                purged_save.write(pair)
        finally:
            purged_save.close()
            pair_path.close()

        lens1, chr_num1 = self.get_lens(self.lens_file1)
        lens2, chr_num2 = self.get_lens(self.lens_file2)

        if not os.path.isdir(f'{path}/pairs'):
            os.mkdir(f'{path}/pairs')
        else:
            shutil.rmtree(f'{path}/pairs')
            # os.system(f'rm -rf {path}/pairs')
            os.mkdir(f'{path}/pairs')
        purged_f = open(purged_path, 'r')
        try:
            for li in purged_f:
                lis = re.split('\\s+', li)
                spe1, chr1 = self.get_chr_num(lis[0])
                spe2, chr2 = self.get_chr_num(lis[3])
                sf = open(f'{path}/pairs/{spe1}.{chr1}.{spe2}.{chr2}.pair', 'a+')
                sf.write(li)
        finally:
            purged_f.close()
            sf.close()

        return lens1, lens2, chr_num1

    def run_colinearscan(self, len1, len2, pair_path, block_path):
        try:
            p = subprocess.call(f'{self.blockscan} -chr1len {len1} -chr2len '
                                f'{len2} -mg1 50 -mg2 50 {pair_path} > {block_path}',
                                shell=True)
            if p != 0:
                print(f'Warning!!! {p}')
        except Exception as e:
            print(e)

    def run(self):
        try:
            path = os.path.basename(self.blast_file).split('.')[0]
            # if not os.path.isdir('blk'):
            #     os.mkdir('blk')
            lens1, lens2, chr_num = self.blast_get_pair(self.blast_file, path)

            thead_pool = ThreadPoolExecutor(self.num_thread)
            files = os.listdir(f'./{path}/pairs/')
            blk_path = f'./{path}/block'
            if not os.path.isdir(blk_path):
                os.mkdir(blk_path)
            thread = []
            for file in files:
                fl = file.split('.')
                if fl[1] not in lens1 or fl[3] not in lens2: continue
                len1, len2 = lens1[fl[1]], lens2[fl[3]]
                blk_name = '.'.join(fl[:-1])+'.blk'
                pair_path, block_path = f'./{path}/pairs/{file}', f'./{path}/block/{blk_name}'
                handle = thead_pool.submit(self.run_colinearscan, len1, len2, pair_path, block_path)
                thread.append(handle)

            for th in as_completed(thread):
                th.result()
            thead_pool.shutdown()
            # sf = f'./blk/{path}.block.rr.txt'

            self.check_del(self.save_block_file)
            block_path = f'./{path}/block/'
            files = os.listdir(block_path)
            for file in files:
                try:
                    file_path = os.path.join(block_path, file)
                    self.remove_redundancy(file_path, self.save_block_file)
                except Exception as e:
                    print(f'error: {e}')

            partiallys = []
            spec1, spec2 = path.split('_')
            blockpath = os.path.dirname(self.save_block_file)
            if spec1 == spec2:
                sf1 = os.path.join(blockpath, f'{path}.partially.block.rr.txt')
                # sf1 = f'{path}.partially.block.rr.txt'
                self.check_del(sf1)
                for i in range(1, chr_num+1):
                    for j in range(i, chr_num+1):
                        partiallys.append(f'{spec1}.{i}.{spec2}.{j}.blk')
                for file in partiallys:
                    file_path = os.path.join(block_path, file)
                    if not os.path.isfile(file_path):
                        continue
                    try:
                        self.remove_redundancy(file_path, sf1)
                    except Exception as e:
                        print(f'error: {e}')
        except PermissionError:
            raise 'Accidental termination!!!'
        except Exception as e:
            raise f'Error: {e}'
        finally:
            shutil.rmtree(f'{path}')
