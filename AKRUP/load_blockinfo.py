from AKRUP.funcbase import *


class BlockInfo:
    def __init__(self, options):

        self.lens_file1 = 'lens1 file'
        self.lens_file2 = 'lens2 file'
        self.gff_file1 = 'gff1 file'
        self.gff_file2 = 'gff2 file'
        self.blast_file = 'blast file'
        self.block_file = 'block file'
        self.ks_file = 'ks file'
        self.save_file = 'block.information.csv'
        self.score = 100
        self.evalue = 1e-5
        self.repnum = 20
 
        for k, v in options:
            if str(k) in ['score', 'evalue', 'repnum']:
                setattr(self, str(k), float(v))
            else:
                setattr(self, str(k), v)

        self.save_file = open(self.save_file, 'w')

    @staticmethod
    def get_chr_num(gene):
        # spec_na = re.findall('^\\D+', gene)[0]
        chr_num = re.search('\\d+', gene)[0]
        return int(chr_num)

    def get_block(self, path):
        da = open(path, 'r')
        blocks = []
        block = []
        for li in da:
            li = li.strip()
            if 'the' in li:
                block.append(li)
            elif '>' in li:
                block.append(li)
                blocks.append(block)
                block = []
            elif 'overlap' in li or 'self' in li or \
                    '#' in li or 'MAXIMUM' in li:
                continue
            elif li:
                lis = re.split('\\s', li)
                block.append(lis)

        return blocks

    def get_gff(self, gff_path, lens_path):
        DotplotBase.check_file(gff_path)
        DotplotBase.check_file(lens_path)
        new_gff, lens_dit, gen_chr, gen_order = [], {}, {}, {}
        lens_f = open(lens_path, 'r')
        for li in lens_f:
            lis = li.strip().split()
            lis[0] = str(self.get_chr_num(lis[0]))
            lens_dit[lis[0]] = lis[1]

        for li in open(gff_path, 'r'):
            lis = li.strip().split()
            gen_chr[lis[5]] = str(self.get_chr_num(lis[0]))
            gen_order[lis[5]] = int(lis[6])
            if str(self.get_chr_num(lis[0])) not in lens_dit:
                continue
            new_line = [str(lis[0]), int(lis[1]), int(lis[2]), str(lis[3]), str(lis[5]), int(lis[6])]
            new_gff.append(new_line)

        return gen_chr, new_gff

    def blast_hocv(self, new_blast):
        score = {}
        for key in new_blast.keys():
            gene_list = new_blast[key]
            pa = {key + "_" + ge: [0] * 8 for ge in gene_list}
            score.update(pa)

        for i in range(1, 9):
            for key in new_blast.keys():
                gene_list = new_blast[key]
                bluenum = 4 + i
                for ge in gene_list[:i]:
                    score[key + "_" + ge][i - 1] = 1
                for ge in gene_list[i:bluenum]:
                    score[key + "_" + ge][i - 1] = 0
                for ge in gene_list[bluenum:]:
                    score[key + "_" + ge][i - 1] = -1
        return score

    def caculate_ks(self, block, pair_ks):
        ks = []
        for lines in block:
            pair = lines[0] + '_' + lines[2]
            pair1 = lines[2] + '_' + lines[0]
            if pair in pair_ks:
                ks.append(pair_ks[pair])
            elif pair1 in pair_ks:
                ks.append(pair_ks[pair1])
            else:
                ks.append(-1)
        ks_va = '_'.join([str(x) for x in ks])
        ks = sorted(ks)
        new_ks = list(filter(lambda x: float(x) > 0, ks))
        if len(new_ks) == 0:
            return ks_va, 0, 0
        ks_median = np.median(new_ks)
        ks_mean = np.mean(new_ks)

        return ks_va, ks_median, ks_mean

    def get_ks(self, ks_file):
        ks_f = open(ks_file, 'r')
        pair_ks = {}
        for row in ks_f:
            rows = row.strip().split()
            pair = f"{rows[0]}_{rows[1]}"
            pair_ks[pair] = float(rows[3])

        return pair_ks

    def extract_blockinfo(self, blocks, blast_score, new_gff1, ge_chr1, new_gff2, ge_chr2, pair_ks, num_or=1):
        self.save_file.write(','.join(['num', 'chr1', 'chr2', 'start1', 'end1', 'start2', 'end2', 'pvalue', 'length',
                             'ks_median', 'ks_average', 'hocv1', 'hocv2', 'hocv3', 'hocv4', 'hocv5', 'hocv6', 'hocv7', 
                             'hocv8', 'block1', 'block2', 'ks', 'density1', 'density2', 'CSR', 'module', 'Ts\n']))

        for block in blocks:
            gff1_ges = {x[4]: 1 for x in new_gff1}
            gff2_ges = {x[4]: 1 for x in new_gff2}
            if block[1][0] not in gff1_ges or block[1][2] not in gff2_ges:
                continue
            chr1, chr2 = ge_chr1[block[1][0]], ge_chr2[block[1][2]]
            length = len(block[1:-1])
            p_value = re.findall('p-value : (\\d+\\.\\d+)', block[-1])[0]
            ks_va, ks_median, ks_mean = self.caculate_ks(block[1:-1], pair_ks)
            block_order1 = [str(int(float(pos[1]))) for pos in block[1:-1]]
            block_order2 = [str(int(float(pos[3]))) for pos in block[1:-1]]
            start1, end1 = block_order1[0], block_order1[-1]
            start2, end2 = block_order2[0], block_order2[-1]

            bk_pos1 = '_'.join(block_order1)
            bk_pos2 = '_'.join(block_order2)
            dedsity1 = length / (abs(int(end1) - int(start1)) + 1)
            dedsity2 = length / (abs(int(end2) - int(start2)) + 1)
            pairs = [f'{x[0]}_{x[2]}' for x in block[1:-1]]
            hocvs = [-1, -1, -1, -1, -1, -1, -1, -1]
            for i in range(1, 9):
                all_sco = [blast_score[x][i - 1] for x in pairs if x in blast_score]
                hocvs[i - 1] = np.mean(all_sco)
            hocvs = [str(v) for v in hocvs]
            rows = [str(v) for v in [num_or, chr1, chr2, start1, end1, start2, end2, p_value, length, ks_median, 
                                     ks_mean, hocvs[0], hocvs[1], hocvs[2], hocvs[3], hocvs[4], hocvs[5], hocvs[6], 
                                     hocvs[7], bk_pos1, bk_pos2, ks_va, dedsity1, dedsity2, '0', '0', '0']]
            self.save_file.write(','.join(rows) + '\n')
            num_or += 1

    def run(self):
        blocks = self.get_block(self.block_file)
        gen_chr1, new_gff1 = self.get_gff(self.gff_file1, self.lens_file1)
        gen_chr2, new_gff2 = self.get_gff(self.gff_file2, self.lens_file2)
        gff1_ges = {x[4]: 1 for x in new_gff1}
        gff2_ges = {x[4]: 1 for x in new_gff2}
        newblast = DotplotBase.getnewblast(self.blast_file, self.score, self.evalue, self.repnum, gff1_ges, gff2_ges)
        blast_score = self.blast_hocv(newblast)
        pair_ks = self.get_ks(self.ks_file)
        self.extract_blockinfo(blocks, blast_score, new_gff1, gen_chr1, new_gff2, gen_chr2, pair_ks)
