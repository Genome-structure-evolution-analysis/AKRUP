from AKRUP.funcbase import *


class PolyploidyEvent:
    def __init__(self, options):
        self.corr_file = '*.top.correspondence.txt'
        self.blockinfo = '*.EventRelate_block.information.csv'
        self.save_file = '*.Polyploidy-block.information.csv'
        
        for k, v in options:
            setattr(self, str(k), v)
        self.sf_handle = open(self.save_file, 'w')

    @staticmethod
    def get_chr_num(gene):
        spec_na = re.findall('^\\D+', gene)[0]
        chr_num = re.search('\\d+', gene)[0]
        return str(int(chr_num))

    def classific_event_module(self, blockinfos, corr_file):
        middle_moduleinfo = []
        handle_corr = [x.strip() for x in open(corr_file)]
        middle_moduleinfo.append(blockinfos[0])

        module_num = 1
        for li in handle_corr:
            rows = li.strip().split()
            for block in rows:
                pos = block.split(':')
                chr1, chr2 = self.get_chr_num(pos[1]), self.get_chr_num(pos[3])
                s1, e1, s2, e2 = [int(x) for x in pos[2].split('-')] + [int(x) for x in pos[4].split('-')]
                for b_line in blockinfos:
                    b_lines = b_line.strip().split(',')
                    if not all((b_lines[1] == chr1, b_lines[2] == chr2)):
                        continue
                    s11, e11, s22, e22 = float(b_lines[3]), float(b_lines[4]), float(b_lines[5]), float(b_lines[6])
                    CSR_num = self.get_chr_num(pos[0])
                    if s11 >= s1 and e11 <= e1 and s22 >= s2 and e22 <= e2:
                        b_lines[24] = CSR_num
                        b_lines[25] = str(module_num)
                        if len(pos) == 6:
                            b_lines[26] = pos[5]
                        middle_moduleinfo.append(','.join(b_lines))
                module_num += 1
        return middle_moduleinfo

    def run(self):
        block_infos = [line.strip() for line in open(self.blockinfo)]
        block_infos = self.classific_event_module(block_infos, self.corr_file)
        for i in range(len(block_infos)):
            if i == 0:
                self.sf_handle.write(block_infos[i]+'\n')
            else:
                lis = block_infos[i].split(',')
                lis[0] = str(i)
                self.sf_handle.write(','.join(lis)+'\n')
