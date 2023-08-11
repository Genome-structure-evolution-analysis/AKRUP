from AKRUP.funcbase import *

class KaryotypeFig(DotplotBase):
    def __init__(self, options):
        DotplotBase.__init__(self)
        self.name = ''
        self.ancestor_file = ''
        self.anc_reverse = 'TRUE'
        self.frame_flag = 'TRUE'
        self.savefile = 'ancestral-result.png'  # png pdf svg
        self.chr_name = 'TRUE'
        
        for k, v in options:
            setattr(self, k, v)

        if not self.ancestor_file:
            raise AssertionError("Please check the ancestor file")

        self.align = dict(family='Times New Roman', horizontalalignment="center", verticalalignment="center")  # Arial

    def draw_frame(self, ancestor_list, step, column_gap, left_curb, top_curb, max_height, max_chrlens):
        ancestor_list = sorted(ancestor_list, key=lambda x: x[0])
        max_chr = max([int(x[0]) for x in ancestor_list])
        frame = []
        if self.anc_reverse.upper() == 'TRUE':
            start_height = max_height+top_curb-2*top_curb-max_chrlens*step

        for anc in ancestor_list:
            num_chr = int(anc[0])
            middlens = anc[2]-anc[1]
            if self.anc_reverse.upper() == 'TRUE':
                y = start_height+step*(anc[1]-1)
            else:
                y = max_height+top_curb-2*top_curb-step*anc[2]
            x = left_curb+column_gap*num_chr+column_gap*2*(num_chr-1)
            frame.append([x, y, column_gap*2, step*middlens, anc[3]])
            if self.chr_name.upper() == 'TRUE' and (anc[0] == '1' or anc[0] == str(max_chr)):
                plt.text(x+column_gap, max_height, str(anc[0]), color='black', rotation=0, **self.align, fontstyle='normal', fontsize=15)
        
        return frame

    def run(self):
        ancestor_list = []
        with open(self.ancestor_file) as ancfile:
            for li in ancfile:
                lis = li.strip().split()
                lis[1] = int(lis[1])
                lis[2] = int(lis[2])
                ancestor_list.append(lis)

        max_chrlens = max([x[2] for x in ancestor_list])
        chr_num = len(set([x[0] for x in ancestor_list]))
        chr_lens = {x[0]: x[2] for x in ancestor_list}

        left_curb, top_curb, max_height, max_width = 0.05, 0.03, 0.9, 0.9
        step, column_gap = (max_height-top_curb*3)/max_chrlens, max_width/(chr_num*3+1)
        fig, ax = plt.subplots(figsize=(10/15*chr_num, 8), dpi=500)

        if self.frame_flag.upper() == 'TRUE':
            self.Rectangle(ax, (left_curb, top_curb), max_width, max_height, color=None, edgecolor='black', alpha=1, linew=1, flag=False)
            plt.text(max_width*0.5+left_curb, max_height+top_curb*2, self.name, color='black',
                 fontsize=20 if chr_num > 2 else 17, **self.align, fontstyle='italic')
        else:
            plt.text(max_width*0.5+left_curb, max_height+top_curb, self.name, color='black',
                     fontsize=20 if chr_num > 2 else 17, **self.align, fontstyle='italic')

        frame_loc = self.draw_frame(ancestor_list, step, column_gap, left_curb, top_curb, max_height, max_chrlens)
        for x in frame_loc:
            self.Rectangle(ax, x[:2], x[2], x[3], x[4], 1)

        if self.chr_name.upper() == 'TRUE':
            x1 = left_curb+column_gap*(2-1)+column_gap*2*(2-1)
            x2 = left_curb+column_gap*(chr_num)+column_gap*2*(chr_num-1)
            ax.plot([x1, x2], [max_height, max_height], '-', markersize=14, lw=1, color='black')
            ax.plot([x2, x2], [max_height, max_height], '>', markersize=5, lw=0.15, color='black')

        plt.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0)
        ax.set_axis_off()
        plt.savefig(self.savefile)

