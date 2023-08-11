from AKRUP.funcbase import *


class DotPlotTrajectory(DotplotBase):
    def __init__(self, options):
        DotplotBase.__init__(self)
        # self.draw_type = 'order'
        self.blast_file = 'blast file'
        self.gff_file1 = 'gff1 file'
        self.gff_file2 = 'gff2 file'
        self.process_name = 'WGD name/num fusion'
        self.score = 100
        self.evalue = 1e-5
        self.repnum = 20
        self.hitnum = 5
        self.multiple = 1
        self.top_trajectory_ancestor_file = 'trajectory color (*.ancestor_trajectory_conf.txt)'
        self.top_redefining_ancestor_file = 'redefining color file (*.ancestor_genome_conf.txt)'
        self.left_ancestor_file = 'ancestor color file (*.ancestor_genome_conf.txt)'

        for k, v in options:
            if str(k) in ['score', 'evalue', 'repnum', 'hitnum', 'multiple', 'dpi']:
                setattr(self, str(k), float(v))
            else:
                setattr(self, str(k), v)

    def plot_arrows(self, ax, name):
        ax.plot([0, 0.87], [0.45, 0.45], '-', markersize=1, lw=2, color='red')
        ax.plot([0.87, 0.87], [0.45, 0.45], '>', markersize=5, lw=0.2, color='red')
        ax.text(0.45, 0.465, name, horizontalalignment="center", verticalalignment="center",fontsize=8, color = 'black')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])

    def plot_karyotype(self, left_ancestor_list, top_ancestor_list1, top_ancestor_list2):

        l1 = max(len(set([x[0] for x in left_ancestor_list])), 0)
        l2 = max(len(set([x[0] for x in top_ancestor_list1])), 0)
        l3 = max(len(set([x[0] for x in top_ancestor_list2])), 0)
        weight = 84*7/39/(l1+l2+l3)
        l1 = l1*weight*2
        l2 = l2*weight*2
        l3 = l3*weight*2
        
        if top_ancestor_list1:
            ax4 = plt.axes([(1+1.5/7*l1)/12, 0.004, 1/12, 0.05])
            self.plot_arrows(ax4, self.process_name)
        if top_ancestor_list2:
            ax5 = plt.axes([(1+1.5/7*l1+3/12*l2+1)/12, 0.004, 1/12, 0.05])
            self.plot_arrows(ax5, 'Redefining')

        ax3 = plt.axes([1/12, 0.004, (1.5/7*l1)/12, 0.05])  #位置[左,下,右,上]
        self.plot_chromosome_fig(ax3, 0.5, 1, left_ancestor_list)
        ax2 = plt.axes([(1+1.5/7*l1+1)/12, 0.004, (3/12*l2)/12, 0.05])
        self.plot_chromosome_fig(ax2, 0.5, 1, top_ancestor_list1)
        ax1 = plt.axes([(1+1.5/7*l1+3/12*l2+2)/12, 0.004, (3/12*l3)/12, 0.05])
        self.plot_chromosome_fig(ax1, 0.5, 1, top_ancestor_list2)

    def run(self):
        plt.figure(figsize=(8, 8), dpi=self.dpi)
        root = plt.axes([0, 0, 1, 1])
        gl1, gl2 = 5/6, 5/6
        gene_loc_1, step1, chr1_lens, spec1_chr_dict = self.gene_location(gl1, self.lens_file1, self.gff_file1)
        gene_loc_2, step2, chr2_lens, spec2_chr_dict = self.gene_location(gl2, self.lens_file2, self.gff_file2)

        align1 = dict(color='black', fontsize=18, rotation=90,style='italic')
        align2 = dict(color='black', fontsize=18, rotation=0,style='italic')
        chr_alian = dict(color='black', fontsize=14, rotation=0, style='normal')
        self.plot_chr1(chr1_lens, gl1, gl2, step1, '', self.genome1_name, align1, chr_alian, 2240/2400, 1/12, 17/240)
        self.plot_chr2(chr2_lens, gl1, gl2, step2, '', self.genome2_name, align2, chr_alian, 1/12, 224/240, 225/240)

        blast = self.getnewblast(self.blast_file, self.score, self.evalue, self.repnum, gene_loc_1, gene_loc_2)
        x, y, colors = self.pair_positon(blast, gene_loc_1, gene_loc_2, 224/240, 1/12)
        plt.scatter(x, y, s=0.5, c=colors, alpha=0.5, edgecolors=None, linewidths=0, marker='o')

        left_ancestor_list, top_ancestor_list1, top_ancestor_list2 = [], [], []
        if self.top_trajectory_ancestor_file:
            top_ancestor_list1 = self.get_color_location(root, spec1_chr_dict, spec2_chr_dict, step1, step2, 
                                        self.top_trajectory_ancestor_file, 'top1', 224/240, 1/12, 205/2400)
        if self.top_redefining_ancestor_file:
            top_ancestor_list2 = self.get_color_location(root, spec1_chr_dict, spec2_chr_dict, step1, step2, 
                                        self.top_redefining_ancestor_file, 'top2', 224/240, 1/12, 205/2400)
        if self.left_ancestor_file:
            left_ancestor_list = self.get_color_location(root, spec1_chr_dict, spec2_chr_dict, step1, step2, 
                                        self.left_ancestor_file, 'left', 224/240, 1/12, 205/2400)

        self.plot_karyotype(left_ancestor_list, top_ancestor_list1, top_ancestor_list2)

        root.set_xlim(0, 1)
        root.set_ylim(0, 1)
        root.set_axis_off()
        plt.savefig(self.savefile)
