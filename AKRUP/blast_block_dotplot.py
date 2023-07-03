from AKRUP.funcbase import *


class DotPlot(DotplotBase):
    def __init__(self, options, section):
        DotplotBase.__init__(self)
        self.section = section
        # self.draw_type = 'order'
        self.blast_file = 'blast file'
        self.block_file = 'block file'
        self.gff_file1 = 'gff1 file'
        self.gff_file2 = 'gff2 file'
        self.score = 100
        self.evalue = 1e-5
        self.repnum = 20
        self.hitnum = 5
        self.multiple = 1

        self.block_num = 0

        for k, v in options:
            if str(k) in ['score', 'evalue', 'repnum', 'hitnum', 'multiple', 'dpi', 'block_num']:
                setattr(self, str(k), float(v))
            else:
                setattr(self, str(k), v)


    def select_blast(self, blast_path, block_path):
        gene_pair_list, pair_pos, last_na1, last_na2, result_pair, num = [], {}, '', '', [], 1
        gene_dit = self.get_block(block_path)
        for li in open(blast_path):
            lis = li.strip().split()
            if 'Scaffold' in li or lis[0] == lis[1]:
                continue
            if lis[0] in gene_dit and lis[1] in gene_dit[lis[0]]:
                gene_pair_list.append(f'{lis[0]}\t{lis[1]}')
            if last_na1 == lis[0] and last_na2 != lis[1]:
                num += 1
            else:
                num = 1
            pair_pos[f'{lis[0]}\t{lis[1]}'] = str(num)
            last_na1, last_na2 = lis[0], lis[1]

        for pair in gene_pair_list:
            if pair in pair_pos:
                row = pair.split('\t')
                row.append(pair_pos[pair])
                result_pair.append(row)

        return result_pair

    def get_block(self, path):
        gene_dit = {}
        block, blocks = [], []
        for li in open(path):
            li = li.strip()
            if '>' in li:
                if len(block) >= self.block_num:
                    for row in block:
                        rows = re.split('\\s', row)
                        if rows[0] in gene_dit:
                            gene_dit[rows[0]].append(rows[2])
                        else:
                            gene_dit[rows[0]] = [rows[2]]
                block = []
            elif any(([x in li for x in ['overlap', 'LOCALE',
                                         'self', 'the', '>']])):
                continue
            elif li:
                block.append(li)

        return gene_dit

    def getnewblast1(self, blast, loc_1, loc_2):
        newblast = []
        for lis in blast:
            if not all((lis[0] in loc_1, lis[1] in loc_2)):
                continue
            newblast.append(lis)
        return newblast

    def pair_positon1(self, data, loc1, loc2, colors):
        pos1, pos2, newcolor = [], [], []
        gl_start1, gl_start2 = 11 / 12, 1 / 12
        for row in data:
            x = gl_start1 - loc1[row[0]]
            y = gl_start2 + loc2[row[1]]

            if row[2] == '1':
                color = colors[0]
                pos1.append(y)
                pos2.append(x)
                newcolor.append(color)
            elif row[2] == '2':
                color = colors[1]
                pos1.append(y)
                pos2.append(x)
                newcolor.append(color)

        return pos1, pos2, newcolor

    def run(self):
        plt.figure(figsize=(8, 8), dpi=self.dpi)
        root = plt.axes([0, 0, 1, 1])
        gl1, gl2 = 5/6, 5/6
        gene_loc_1, step1, lens_1, spec1_chr_dict = self.gene_location(gl1, self.lens_file1, self.gff_file1)
        gene_loc_2, step2, lens_2, spec2_chr_dict = self.gene_location(gl2, self.lens_file2, self.gff_file2)

        align1 = dict(color='black', fontsize=18, rotation=90,style='italic')
        align2 = dict(color='black', fontsize=18, rotation=0,style='italic')
        chr_alian = dict(color='black', fontsize=14, rotation=0, style='normal')
        self.plot_chr1(lens_1, gl1, gl2, step1, '', self.genome1_name, align1, chr_alian)
        self.plot_chr2(lens_2, gl1, gl2, step2, '', self.genome2_name, align2, chr_alian)

        if self.section == 'dotplot':
            blast = self.getnewblast(self.blast_file, self.score, self.evalue, self.repnum, gene_loc_1, gene_loc_2)
            x, y, colors = self.pair_positon(blast, gene_loc_1, gene_loc_2)
            plt.scatter(x, y, s=0.5, c=colors, alpha=0.5, edgecolors=None, linewidths=0, marker='o')

        elif self.section == 'blockdotplot':
            blast = self.select_blast(self.blast_file, self.block_file)
            blast = self.getnewblast1(blast, gene_loc_1, gene_loc_2)
            x, y, colors = self.pair_positon1(blast, gene_loc_1, gene_loc_2, self.colors)
            plt.scatter(x, y, s=1, c=colors, alpha=1, edgecolors=None, linewidths=0)

        root.set_xlim(0, 1)
        root.set_ylim(0, 1)
        root.set_axis_off()
        plt.savefig(self.savefile)
