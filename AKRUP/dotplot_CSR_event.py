from scipy import stats
from bisect import bisect
from AKRUP.funcbase import *


class DotplotCsrKs(DotplotBase):
    def __init__(self, options, section):
        DotplotBase.__init__(self)
        
        self.block_info = 'block info file'
        self.top_ancestor_file = ''  # A.*.top.color.pos.txt
        self.left_ancestor_file = ''  # A.*.left.color.pos.txt
        self.hocv_depth = 1
        self.hocv = -1
        self.block_num = 5
        self.ks_s, self.ks_e = 0, 3
        self.range_k = 0.3
        self.peaks = ''  # a,u,sigma
        self.pkcolor = 'orange'
        self.peakflag = 'TRUE'

        for k, v in options:
            if str(k) in ['hocv_depth', 'hocv', 'dpi', 'block_num', 'range_k']:
                setattr(self, str(k), float(v))
            else:
                setattr(self, str(k), v)


        if self.peaks:
            self.peaks = [float(x) for x in self.peaks.split(',')]
        self.hocv_depth = int(self.hocv_depth)
        self.section = section

    def get_chr_lens(self, lens_file, gl):
        len_pos, chr_dict, chr_lens, n = 1, {}, {}, 0
        self.check_file(lens_file)
        for li in open(lens_file):
            lis = li.strip().split()
            spec_chr = self.get_spec_chr(lis[0])
            chr_lens[spec_chr] = float(lis[len_pos])
            chr_dict[spec_chr] = float(n)
            n += float(lis[len_pos])
        total_lens = reduce(lambda x, y: int(x)+int(y), chr_lens.values())
        step = gl/total_lens

        return step, chr_lens, chr_dict

    def read_blockinfo(self, block_info_file, bk_num, hocv, ks_s, ks_e, hocv_depth):
        new_block, hocv_pos = [], 10
        hocv_pos += hocv_depth
        block_list = [x.strip().split(',') for x in open(block_info_file)]
        for bk in block_list[1:]:
            if not int(bk[8]) >= bk_num:
                continue
            if float(bk[9]) < ks_s or float(bk[9]) > ks_e:
                continue
            if float(bk[hocv_pos]) >= hocv:
                new_block.append(bk)

        return new_block

    def pair_positon(self, location):
        pos1, pos2, newcolor = [], [], []
        gl_start1, gl_start2 = 11 / 12, 1 / 12
        for row in location:
            x = gl_start1 - row[0]
            y = gl_start2 + row[1]
            pos1.append(y)
            pos2.append(x)
            newcolor.append(row[2])

        return pos1, pos2, newcolor

    def pair_positon1(self, location, color, peaks, rang):
        pos1, pos2, newcolor = [], [], []
        gl_start1, gl_start2 = 11 / 12, 1 / 12
        for row in location:
            x = gl_start1 - row[0]
            y = gl_start2 + row[1]
            pos1.append(y)
            pos2.append(x)
            if peaks-rang <= row[2] <= peaks+rang:
                newcolor.append(color)
            else:
                newcolor.append('gray')

        return pos1, pos2, newcolor

    def plot_curve(self, data, ax2, range_k, color):

        # t=np.arange(data[1]-0.8,data[1]+0.8,0.01)
        t=np.arange(data[1]-1,data[1]+1,0.01)
        # for k in data:
        y=[0 for k in t]
            # for i in range(0,int((len(k)-1)/3)):    
        y1=stats.norm.pdf(t,float(data[1]),float(data[2]))*float(data[0])
        y = [i + j for i, j in zip(y, y1)]
        style='-'
        ax2.plot(t,y,linestyle=style,color=color,label='0',linewidth=1.5)
        ax2.vlines([data[1]], 0, max(y), linestyles='dashed', colors='red')
        maxy1, maxy2 = y[bisect(t, data[1]+range_k)], y[bisect(t, data[1]-range_k)]
        ax2.fill(t, y, facecolor=color, color=color)
        plt.text(data[1], max(y)+0.06, f'μ={round(data[1], 4)}', color='black', fontsize=10, family='Times New Roman', horizontalalignment="center")
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        plt.xticks([])
        plt.yticks([])

    def run(self):
        plt.figure(figsize=(8, 8), dpi=self.dpi)
        root = plt.axes([0, 0, 1, 1])
        gl1, gl2 = 5/6, 5/6
        bk_info_list = self.read_blockinfo(self.block_info, self.block_num, self.hocv, self.ks_s, self.ks_e, self.hocv_depth)

        step1, lens_1, spec1_chr_dict = self.get_chr_lens(self.lens_file1, gl1)
        step2, lens_2, spec2_chr_dict = self.get_chr_lens(self.lens_file2, gl2)

        align1 = dict(color='black', fontsize=18, rotation=90,style='italic')
        align2 = dict(color='black', fontsize=18, rotation=0,style='italic')
        chr_alian = dict(color='black', fontsize=14, rotation=0, style='normal')
        self.plot_chr1(lens_1, gl1, gl2, step1, '', self.genome1_name, align1, chr_alian)
        self.plot_chr2(lens_2, gl1, gl2, step2, '', self.genome2_name, align2, chr_alian)

        # genes_loc = self.get_gene_location(bk_info_list, spec1_chr_dict, spec2_chr_dict, step1, step2)
        if self.top_ancestor_file:
            self.get_color_location(root, spec1_chr_dict, spec2_chr_dict, step1, step2, self.top_ancestor_file, 'top1')
        if self.left_ancestor_file:
            self.get_color_location(root, spec1_chr_dict, spec2_chr_dict, step1, step2, self.left_ancestor_file, 'left')
        
        if self.section == 'csrdotplot':
            genes_loc = self.get_gene_location(bk_info_list, spec1_chr_dict, spec2_chr_dict, step1, step2, 24)
            x, y, colors = self.pair_positon(genes_loc)  # 根据CSR定义颜色
            cm = plt.cm.get_cmap('gist_rainbow')  # gist_rainbow winter Paired tab10 jet
            sc = plt.scatter(x, y, s=float(1), c=colors,
                             alpha=1, edgecolors=None, linewidths=0, marker='o', vmin=0, vmax=max(colors), cmap=cm)

        elif self.section == 'eventdotplot':
            genes_loc = self.get_gene_location(bk_info_list, spec1_chr_dict, spec2_chr_dict, step1, step2, 9)
            x, y, colors = self.pair_positon1(genes_loc, self.pkcolor, self.peaks[1], self.range_k)  # 根据median
            plt.scatter(x, y, s=1.5, c=colors, alpha=0.5, edgecolors=None, linewidths=0, marker='o')
            
            if self.peakflag.upper() == 'TRUE':
                plt.text(0.05, 0.044, 'Speciation', color='black', fontsize=9, family='Times New Roman', verticalalignment="center", rotation=90)
                ax2 = plt.axes([1/12, 0.016, 1.5/12, 0.05])
                self.plot_curve(self.peaks, ax2, self.range_k, self.pkcolor)

        root.set_xlim(0, 1)
        root.set_ylim(0, 1)
        root.set_axis_off()
        plt.savefig(self.savefile)
