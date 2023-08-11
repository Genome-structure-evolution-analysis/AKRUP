import os
import re
import time
import random
import shutil
import platform
import threading
import configparser
from math import sqrt, pi
from functools import reduce
from itertools import groupby

from tempfile import NamedTemporaryFile
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

import AKRUP


if os.name == "nt":
    os.system("")

platform_name = platform.platform().split('-')[0]


current_path = AKRUP.__path__[0]
ini_path = os.path.join(current_path, 'ini')
if platform_name == 'Linux':
    os.system(f'chmod -R 775 {ini_path}')
    # os.chmod(ini_path, 0o775)


def load_conf(file, section):
    conf = configparser.ConfigParser()
    conf.read(file, 'utf-8')
    return conf.items(section)

def get_temname(mode):
    fp = NamedTemporaryFile(mode=mode, delete=False)
    return fp

def get_temdir():
    a = time.localtime()
    time_list = [str(i) for i in a[2:6]]
    out = "".join(time_list) + str(random.randint(0, 100))

    return out

def await_run(content, event):
    list_circle = ["\\", "|", "/", "—"]
    i = 1
    while not event.isSet():
        time.sleep(0.25)
        # print("\r{} {}".format(content, list_circle[i % 4]), end="")
        print("\r\033[0;32m{}....{}\033[0m".format(content, list_circle[i % 4]), end="")
        i += 1

class DotplotBase:
    def __init__(self):
        self.genome1_name = 'left name'
        self.genome2_name = 'top name'
        self.lens_file1 = 'lens1 file'
        self.lens_file2 = 'lens2 file'
        self.savefile = 'savefile'  # png pdf svg

        self.colors = ['red', 'blue', 'white']
        self.align = dict(family='Times New Roman', horizontalalignment="center", 
                    verticalalignment="center", weight='semibold')

    @staticmethod
    def get_spec_chr(spec_chr):
        spec_chr = re.sub('^\\D+', '', spec_chr)
        spec_chr = re.sub('^0', '', spec_chr)
        return spec_chr

    @classmethod
    def check_file(self, file_path):
        n = os.path.isfile(file_path)
        file_name = os.path.basename(file_path)
        if not n:
            print(f'{file_name}---File not exist!!!')

    def gene_location(self, gl, lens_file, gff_file):
        len_pos, gff_pos = 1, 6
        chr_dict, chr_lens, loc_gene, n = {}, {}, {}, 0
        self.check_file(lens_file)
        for li in open(lens_file):
            lis = li.strip().split()
            spec_chr = self.get_spec_chr(lis[0])
            chr_lens[spec_chr] = float(lis[len_pos])
            chr_dict[spec_chr] = float(n)
            n += float(lis[len_pos])
        total_lens = reduce(lambda x, y: int(x)+int(y), chr_lens.values())
        step = gl/total_lens
        self.check_file(gff_file)
        for li in open(gff_file):
            lis = li.strip().split()
            spec_chr = self.get_spec_chr(lis[0])
            if spec_chr not in chr_dict:
                continue
            loc = (chr_dict[spec_chr] + float(lis[gff_pos]))*step
            loc_gene[lis[5]] = loc

        return loc_gene, step, chr_lens, chr_dict

    @classmethod
    def getnewblast(self, blast, score, evalue, repnum, loc_1, loc_2):
        newblast = {}
        for li in open(blast):
            lis = li.strip().split()
            if not all((float(lis[11]) >= score, float(lis[10]) < evalue, lis[0] != lis[1])):
                continue
            if not all((lis[0] in loc_1, lis[1] in loc_2)):
                continue
            if lis[0] in newblast and lis[1] in {x[0]: 1 for x in newblast[lis[0]]}:
                continue
            if lis[0] in newblast and len(newblast[lis[0]]) < repnum:
                newblast[lis[0]].append([lis[1], lis[11]])
            else:
                newblast[lis[0]] = [[lis[1], lis[11]]]
        for key in newblast.keys():
            gene_list = newblast[key]
            new_gene_list = sorted(gene_list, key=lambda x: [1])
            new_gene_list = [ge[0] for ge in new_gene_list]
            newblast[key] = new_gene_list

        return newblast

    def pair_positon(self, blast, loc1, loc2, gl_start1=11/12, gl_start2=1/12):
        pos1, pos2, newcolor = [], [], []
        # gl_start1, gl_start2 = 11/12, 1/12
        for k, v in blast.items():
            for i in range(len(v)):
                if i < self.multiple:
                    color = self.colors[0]
                elif i <= self.hitnum+self.multiple:
                    color = self.colors[1]
                else:
                    color = self.colors[2]
                pos1.append(gl_start2 + loc2[v[i]])
                pos2.append(gl_start1 - loc1[k])
                newcolor.append(color)
        return pos1, pos2, newcolor

    @staticmethod
    def plot_line(x, y):
        plt.plot(x, y, linestyle='-', color='black', linewidth=0.25)
        plt.plot(x, y, linestyle='-', color='black', linewidth=0.75, alpha=0.5)

    def plot_chr1(self, lens, gl, gl2, step, mark, name, name_align, chrnum_align, 
                    gl_start=11/12, start_x=1/12, mark_y=17/240):
        n = 0
        for k in lens.keys():
            n += lens[k]
            mark_new = str(mark) + str(k)
            x = gl_start - float(n) * step
            mark_x = x + 0.5 * lens[k] * step
            self.plot_line([start_x, start_x + gl2], [x, x])
            plt.text(mark_y-0.01, mark_x, mark_new, **self.align, **chrnum_align)
        self.plot_line([start_x, start_x + gl2], [gl_start, gl_start])
        plt.text(mark_y-0.04, 0.5 * (2 * gl_start - gl), name, **self.align, **name_align)

    def plot_chr2(self, lens, gl, gl2, step, mark, name, name_align, chrnum_align, 
                    gl_start=1/12, start_x=11/12, mark_y=223/240):
        n = 0
        for k in lens.keys():
            n += lens[k]
            mark_new = str(mark) + str(k)
            x = gl_start + float(n) * step
            mark_x = x - 0.5 * lens[k] * step
            self.plot_line([x, x], [start_x, start_x - gl2])
            plt.text(mark_x, mark_y+0.005, mark_new,**self.align, **chrnum_align)
        self.plot_line([gl_start, gl_start], [start_x, start_x - gl2])
        plt.text(0.5 * (2 * gl_start + gl), mark_y + 0.04, name, **self.align, **name_align)

    @staticmethod
    def Rectangle(ax, loc, width, height, color, alpha, linew=None, style='-', flag=True, edgecolor=None):
        p = mpatches.Rectangle(loc, width, height, edgecolor=None, 
            facecolor=color, alpha=alpha, lw=linew, ls=style, fill=flag)
        ax.add_patch(p)

    def get_color_location(self, ax, loc1, loc2, step1, step2, color_pos_file, class_type, 
                            gl_start1=11/12, gl_start2=1/12, top_start=165/2400):
        new_lines, ancestor_chr = [], []
        if isinstance(color_pos_file, list):
            color_list = color_pos_file
        else:
            color_list = [x.strip().split() for x in open(color_pos_file)]
        # for li in open(color_pos_file):
        for lis in color_list:
            # lis = li.strip().split()
            lis[1] = int(lis[1])
            lis[2] = int(lis[2])
            ancestor_chr.append(lis)
            length = abs(int(lis[2])-int(lis[1]))+1
            lis.append(length)
            new_lines.append(lis)

        new_lines = sorted(new_lines, key=lambda x: int(x[0]))
        temp_list = groupby(new_lines, lambda x: x[0])
        for name, group in temp_list:
            group_list = list(group)
            max_num = max([x[2] for x in group_list])
            min_num = min([x[1] for x in group_list])
            if class_type == 'top1':
                for lis in group_list:
                    x1 = gl_start2 + (loc2[lis[0]]+int(lis[1]))*step2
                    x2 = gl_start2 + (loc2[lis[0]]+int(lis[2]))*step2
                    if lis[2] == max_num and lis[1] == min_num:
                        self.Rectangle(ax, (x1+2/2400, top_start), lis[4]*step2-2/2400, 30/2400, lis[3], 1)
                    elif lis[2] == max_num:
                        self.Rectangle(ax, (x1, top_start), lis[4]*step2, 30/2400, lis[3], 1)
                    elif lis[1] == min_num:
                        self.Rectangle(ax, (x1+2/2400, top_start), lis[4]*step2, 30/2400, lis[3], 1)
                    else:
                        self.Rectangle(ax, (x1, top_start), lis[4]*step2, 30/2400, lis[3], 1)

            if class_type == 'top2':
                for lis in group_list:
                    x1 = gl_start2 + (loc2[lis[0]]+int(lis[1]))*step2
                    x2 = gl_start2 + (loc2[lis[0]]+int(lis[2]))*step2
                    if lis[2] == max_num and lis[1] == min_num:
                        self.Rectangle(ax, (x1+2/2400, top_start-35/2400), lis[4]*step2-2/2400, 30/2400, lis[3], 1)
                    elif lis[2] == max_num:
                        self.Rectangle(ax, (x1, top_start-35/2400), lis[4]*step2, 30/2400, lis[3], 1)
                    elif lis[1] == min_num:
                        self.Rectangle(ax, (x1+2/2400, top_start-35/2400), lis[4]*step2, 30/2400, lis[3], 1)
                    else:
                        self.Rectangle(ax, (x1, top_start-35/2400), lis[4]*step2, 30/2400, lis[3], 1)

            elif class_type == 'left':
                for lis in group_list:
                    x1 = gl_start1 - (loc1[lis[0]]+int(lis[1]))*step1
                    x2 = gl_start1 - (loc1[lis[0]]+int(lis[2]))*step1
                    if lis[2] == max_num and lis[1] == min_num:
                        self.Rectangle(ax, (2206/2400, x2), 30/2400, lis[4]*step1-1/2400, lis[3], 1)
                    elif lis[2] == max_num:
                        self.Rectangle(ax, (2206/2400, x2), 30/2400, lis[4]*step1, lis[3], 1)
                    elif lis[1] == min_num:
                        self.Rectangle(ax, (2206/2400, x2), 30/2400, lis[4]*step1-1/2400, lis[3], 1)
                    else:
                        self.Rectangle(ax, (2206/2400, x2), 30/2400, lis[4]*step1, lis[3], 1)

        return ancestor_chr

    def plot_chromosome_fig(self, ax, width, max_col, ancestor_list=None, flag=True):
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
        if not ancestor_list:
            return 
        ancestor_lens = pd.DataFrame(ancestor_list)
        ancestor_lens[0] = ancestor_lens[0].astype(str)
        ancestor_lens[3] = ancestor_lens[3].astype(str)
        # ancestor_lens[4] = ancestor_lens[4].astype(int)
        # ancestor_lens[4] = ancestor_lens[4] / ancestor_lens[4].max()
        chrs = ancestor_lens[0].drop_duplicates().to_list()
        ax.bar(chrs, 1, color='white', alpha=0)
        step = 1/ancestor_lens[2].max()
        for index, row in ancestor_lens.iterrows():
            self.Rectangle(ax, [chrs.index(row[0])-width*0.5,
                            row[1]*step], width, (row[2]-row[1])*step, row[3], 1)  # row[4]

        x1 = chrs.index(chrs[0])-width*0.5
        x2 = chrs.index(chrs[-1])+width*0.5
        if flag:
            ax.text(chrs.index(chrs[0]), ancestor_lens[2].max()*step+0.2, str(min([int(x) for x in chrs])), horizontalalignment="center", verticalalignment="center", fontsize = min(8, int(12*1/(max_col*0.5))), color = 'black')
            ax.text(chrs.index(chrs[-1]), ancestor_lens[2].max()*step+0.2, str(max([int(x) for x in chrs])), horizontalalignment="center", verticalalignment="center", fontsize = min(8, int(12*1/(max_col*0.5))), color = 'black')
        else:
            ax.text(chrs.index(chrs[0]), ancestor_lens[2].max()*step+0.1, chrs[0], horizontalalignment="center", verticalalignment="center", fontsize = min(12, int(12*1/(max_col*0.5))), color = 'black')
            ax.text(chrs.index(chrs[-1]), ancestor_lens[2].max()*step+0.1, chrs[-1], horizontalalignment="center", verticalalignment="center", fontsize = min(12, int(12*1/(max_col*0.5))), color = 'black')


        ax.plot([x1, x2], [ancestor_lens[2].max()*step+0.05, ancestor_lens[2].max()*step+0.05], '-', markersize=min(14, int(14*1/(max_col*0.5))), lw=1, color='black')
        ax.plot([x2, x2], [ancestor_lens[2].max()*step+0.05, ancestor_lens[2].max()*step+0.05], '>', markersize=min(5, int(5*1/(max_col*0.5))), lw=0.15, color='black')
        ax.tick_params(labelsize=15)

    def get_gene_location(self, bk_info, len1, len2, step1, step2, pos):
        gene_location, class_num = [], 0
        new_bks = [bk for bk in bk_info if bk[1] in len1 and bk[2] in len2]
        for bk in new_bks:
            pos1 = [int(po) for po in bk[19].split('_')]
            pos2 = [int(po) for po in bk[20].split('_')]
            class_num = float(bk[pos])
            for i, j in zip(pos1, pos2):
                loc1 = (len1[bk[1]] + float(i)) * step1
                loc2 = (len2[bk[2]] + float(j)) * step2
                gene_location.append([loc1, loc2, class_num])

        return gene_location