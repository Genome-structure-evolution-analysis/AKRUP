from AKRUP.funcbase import *


class DotplotBlock(DotplotBase):
    def __init__(self, name, lens_dit_ilst, block_list, top_color_list, left_color_list, name_change, ancestral_color, save_path):
        DotplotBase.__init__(self)
        self.name = name
        self.save_path = save_path
        self.bk_info_list = block_list
        self.name_change = name_change
        self.lens_dit_ilst = lens_dit_ilst
        self.top_color_list = top_color_list
        self.left_color_list = left_color_list
        self.ancestral_color = ancestral_color

        self.align1 = dict(color='black', fontsize=18, rotation=90,style='italic')
        self.align2 = dict(color='black', fontsize=18, rotation=0,style='italic')
        self.chr_alian = dict(color='black', fontsize=14, rotation=0, style='normal')

    def get_chr_lens(self, lens_list, gl):
        len_pos, chr_dict, chr_lens, n = 1, {}, {}, 0
        for lis in lens_list:
            chr_lens[lis[0]] = float(lis[len_pos])
            chr_dict[lis[0]] = float(n)
            n += float(lis[len_pos])
        total_lens = reduce(lambda x, y: int(x)+int(y), chr_lens.values())
        step = gl/total_lens

        return step, chr_lens, chr_dict

    def pair_positon(self, location):
        pos1, pos2, newcolor = [], [], []
        gl_start1, gl_start2 = 11 / 12, 1 / 12
        for row in location:
            x = gl_start1 - row[0]
            y = gl_start2 + row[1]
            pos1.append(y)
            pos2.append(x)
            ty_name = 'CSR'+str(int(row[2]))
            color = self.ancestral_color[ty_name]
            newcolor.append(color)

        return pos1, pos2, newcolor

    def plotfig(self, root, spec1, spec2):
        gl1, gl2 = 5/6, 5/6
        step1, lens_1, spec1_chr_dict = self.get_chr_lens(self.lens_dit_ilst[spec1], gl1)
        step2, lens_2, spec2_chr_dict = self.get_chr_lens(self.lens_dit_ilst[spec2], gl2)
        self.plot_chr1(lens_1, gl1, gl2, step1, '', self.name_change[spec1], self.align1, self.chr_alian)
        self.plot_chr2(lens_2, gl1, gl2, step2, '', self.name_change[spec2], self.align2, self.chr_alian)
        genes_loc = self.get_gene_location(self.bk_info_list, spec1_chr_dict, spec2_chr_dict, step1, step2, 24)
        if self.top_color_list[:1]:
            self.get_color_location(root, spec1_chr_dict, spec2_chr_dict, step1, step2, self.top_color_list, 'top1')
        if self.left_color_list[:1]:
            self.get_color_location(root, spec1_chr_dict, spec2_chr_dict, step1, step2, self.left_color_list, 'left')
        
        x, y, colors = self.pair_positon(genes_loc)
        plt.scatter(x, y, s=0.9, c=colors, alpha=1, edgecolors=None, linewidths=0)

    def plotfig_dot(self, spec1, spec2):
        plt.figure(figsize=(8, 8), dpi=300)
        root = plt.axes([0, 0, 1, 1])
        gl1, gl2 = 5/6, 5/6
        step1, lens_1, spec1_chr_dict = self.get_chr_lens(self.lens_dit_ilst[spec1], gl1)
        step2, lens_2, spec2_chr_dict = self.get_chr_lens(self.lens_dit_ilst[spec2], gl2)
        self.plot_chr1(lens_1, gl1, gl2, step1, '', self.name_change[spec1], self.align1, self.chr_alian)
        self.plot_chr2(lens_2, gl1, gl2, step2, '', self.name_change[spec2], self.align2, self.chr_alian)
        genes_loc = self.get_gene_location(self.bk_info_list, spec1_chr_dict, spec2_chr_dict, step1, step2, 24)
        if self.top_color_list[:1]:
            self.get_color_location(root, spec1_chr_dict, spec2_chr_dict, step1, step2, self.top_color_list, 'top1')
        if self.left_color_list[:1]:
            self.get_color_location(root, spec1_chr_dict, spec2_chr_dict, step1, step2, self.left_color_list, 'left')
        
        x, y, colors = self.pair_positon(genes_loc)
        plt.scatter(x, y, s=0.9, c=colors, alpha=1, edgecolors=None, linewidths=0)
        root.set_xlim(0, 1)
        root.set_ylim(0, 1)
        root.set_axis_off()
        plt.savefig(f'{self.save_path}/A-{spec1}_{spec2}.ancestral_dotplot.png')
        plt.savefig(f'{self.save_path}/A-{spec1}_{spec2}.ancestral_dotplot.pdf')

    def plot_ancestral(self, ancestor_list, width, max_col, spec1, spec2):
        fig, ax = plt.subplots(figsize=(8, 8))
        self.plot_chromosome_fig(ax, width, max_col, ancestor_list, False)
        ax.set_xticks([])
        ax.set_yticks([])
        plt.savefig(f'{self.save_path}/A-{spec1}_{spec2}.karyotype_chromosome.png', dpi=500)
        plt.savefig(f'{self.save_path}/A-{spec1}_{spec2}.karyotype_chromosome.pdf', dpi=500)

    def plot_ancestral_before(self, ancestor_list, width, max_col, spec1, spec2, i):
        fig, ax = plt.subplots(figsize=(8, 8))
        self.plot_chromosome_fig(ax, width, max_col, ancestor_list, False)
        ax.set_xticks([])
        ax.set_yticks([])
        plt.savefig(f'{self.save_path}/A-{spec1}_{spec2}.before-{i+1}-WGD-karyotype_chromosome.png', dpi=500)
        plt.savefig(f'{self.save_path}/A-{spec1}_{spec2}.before-{i+1}-WGD-karyotype_chromosome.pdf', dpi=500)

    def main(self, ancestor_list, wgds_chr_colors):
        if platform_name == 'Linux':
            font = {'family' : 'DejaVu Sans',
                'weight' : 'normal'
                        }
        elif platform_name == 'Windows':
            font = {'family' : 'Times New Roman',
                'weight' : 'normal'
                        }

        spec1, spec2 = self.name.split('_')
        max_col = len(wgds_chr_colors)+1
        self.plotfig_dot(spec1, spec2)
        self.plot_ancestral(ancestor_list, 0.5, max_col, spec1, spec2)

        for i in range(len(wgds_chr_colors)):
            conf_list = wgds_chr_colors[i]
            # ax_wgd = fig.add_subplot(grid[2:3, 1+i:2+i])
            self.plot_ancestral_before(conf_list, 0.5, max_col, spec1, spec2, i)

        fig = plt.figure(figsize=(8, 8), dpi=300)

        grid = plt.GridSpec(3, max_col, hspace=0.05, wspace=0.1)
        ax_main = fig.add_subplot(grid[0:2, 0:max_col])
        
        self.plotfig(ax_main, spec1, spec2)
        ax_ancestor = fig.add_subplot(grid[2:3, 0:1])
        
        self.plot_chromosome_fig(ax_ancestor, 0.5, max_col, ancestor_list, False)
        ax_ancestor.set_xlabel(f'{self.name_change[spec2]} ancestral chromosome', font, size=min(15, int(15*1/(max_col*0.5))))

        for i in range(len(wgds_chr_colors)):
            conf_list = wgds_chr_colors[i]
            ax_wgd = fig.add_subplot(grid[2:3, 1+i:2+i])
            self.plot_chromosome_fig(ax_wgd, 0.5, max_col, conf_list, False)
            ax_wgd.set_xlabel(f'{self.name_change[spec2]} before-{i+1} WGD chromosome', font, size=min(15, int(15*1/(max_col*0.5))))

        plt.subplots_adjust(left=1/24, right=1, top=23/24, bottom=1/12)

        ax_main.set_xlim(0.05, 0.95)
        ax_main.set_axis_off()

        plt.savefig(f'{self.save_path}/A-{spec1}_{spec2}.ancestral.png')
        plt.savefig(f'{self.save_path}/A-{spec1}_{spec2}.ancestral.pdf')