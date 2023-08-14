# -*- coding: utf-8 -*-
# @Author: wjq
# @Date:   2023-02-04 10:53:22
# @Last Modified by:   wjq
# @Last Modified time: 2023-08-13 23:56:39

import copy
from itertools import islice, combinations, product, chain
from collections import Counter
from AKRUP.ancestral_dotplot import *

class ConservedAncestralRegions:
    def __init__(self, options):

        self.species = 'refname_spec1_spec2'
        self.bk_files = 'refspec_spec1:refspec_spec1.Polyploidy-block.information.csv,refspec_spec2:refspec_spec2.Polyploidy-block.information.csv'
        self.len_files ='refspec:refspec.lens,spec1:spec1.lens,spec2:spec2.lens'
        self.recentwgdchr = 'spec1:spec1_recent_wgd_chr.txt,spec2:spec2_recent_wgd_chr.txt'
        self.wgds = 'spec1:wgdnum,spec2:wgdnum (spec2:2_2 WGD occurred twice after disagreement with refspec)'
        self.latin_name = ' refspec:refspec genome name,spec1:spec1 genome name,spec2:spec2 genome name'
        self.hocv_depths='refspec_spec1:num,refspec_spec2:num'
        self.infer_name = ''  # refspec_spec1,refspec_spec2
        self.intergenomicratio = 1 # orthologous synteny ratio (default: 1)
        self.block_num = 5
        self.hocv = -1 # (-1< hocv <1)
        self.select_ref_ancestor_spec = ''  # sepc name

        self.common_wgd = False
        self.infer_wgd_flag = True
        self.conserved_spec = ''  # sepc name
        self.save_path = '.'
        
        self.color_pos = 0
        for k, v in options:
            if not v:
                continue
            if str(k) in ['hocv']:
                setattr(self, str(k), float(v))
            elif str(k) in ['intergenomicratio', 'block_num']:
                setattr(self, str(k), int(v))
            else:
                setattr(self, str(k), v)

        self.wgds = {y[0]:y[1] for y in [x.split(':') for x in self.wgds.split(',')]}
        self.hocv_depths = {y[0]:int(y[1]) for y in [x.split(':') for x in self.hocv_depths.split(',')]}
        self.latin_name = {y[0]:y[1] for y in [x.split(':') for x in self.latin_name.split(',')]}
        self.bk_files = {y[0]:y[1] for y in [x.split(':') for x in self.bk_files.split(',')]}
        self.len_files = {y[0]:y[1] for y in [x.split(':') for x in self.len_files.split(',')]}
        self.recentwgdchr = {y[0]:y[1] for y in [x.split(':') for x in self.recentwgdchr.split(',')]}
        self.infer_name = self.infer_name.strip().split(',')

        if self.infer_wgd_flag.upper() == 'TRUE' and not self.infer_name[0]:
            raise AssertionError("infer_name and infer_wgd_flag must be consistent")

        self.colors = {'1':'#a6c952', '2':'#2f7ec2','3':'#009c75', "4": '#00b5ef',
                        "5": '#eb2690', "6": '#ffd00a', "7": '#eb2d2d', "8": "#fc8d62", 
                        "9": "orange", "10": "#CCCC99", "11": "#CC99FF", "12": "c",
                        "13": "#CC9900", "14": "#1E90FF", "15": "#32CD32", "16": "#FFA07A",
                        "17": "#BDB76B", "18": "#D2691E", "19": "#7B68EE", "20": "#CC66FF",
                        "21": "#FF1493", "22": "#FFD700", "23": "#FF69B4", "24": "#CCFF33", 
                        '25':'#1b9e77','26':'#decbe4'}
    
    def read_lens(self, lens_files):
        spec_lens, lens_list = {}, {}
        for name, lenfile in lens_files.items():
            midd_lens, midd_len_list = {}, []
            for li in open(lenfile):
                lis = li.strip().split()
                schr = int(DotplotBase.get_spec_chr(lis[0]))
                midd_lens[schr] = int(lis[1])
                lis[0] = str(schr)
                midd_len_list.append(lis)
            spec_lens[name] = midd_lens
            lens_list[name] = midd_len_list
        return spec_lens, lens_list

    def count_CSR_lens(self, new_mdinfos):
        # 统计每个物种所有染色体中 相同CSR的总长度  用于推断祖先染色体的比例（取出最长的作为标准 1）
        # 统计每个物种中每个染色体中相同CSR的比例  用于推断加倍前后染色体的比例
        CSR_top_all_lens, CSR_top_each_len = {}, {}

        for key in new_mdinfos.keys():
            spec_mdinfo = new_mdinfos[key]
            midd_top_all_len, midd_top_each_len = {}, {}
            for lis in spec_mdinfo:
                midd_top_all_len[lis[14]] = midd_top_all_len.get(lis[14], 0) + int(lis[7])
            sort_spec_mdinfo = sorted(spec_mdinfo, key=lambda x: x[1])
            temp_list = groupby(sort_spec_mdinfo, lambda x: x[1])
            for name, group in temp_list:
                group_list = list(group)
                chr_CSR_lens = {}
                for lis in group_list:
                    chr_CSR_lens[lis[14]] = chr_CSR_lens.get(lis[14], 0) + int(lis[7])
                midd_top_each_len[name] = chr_CSR_lens

            CSR_top_all_lens[key] = midd_top_all_len
            CSR_top_each_len[key] = midd_top_each_len
        return CSR_top_all_lens, CSR_top_each_len

    def product_ancestral_chromosome(self, CSR_color_lens, color_order):
        # color_order = []
        # if len(color_order)
        # order_color = [x[3] for x in CSR_color_lens]
        # order_color = sorted(list(set(order_color)),key=order_color.index)
        new_color_order = sorted(color_order.items(), key=lambda x: x[1])
        order_color = [x[0] for x in new_color_order]

        sort_CSR_color_lens = sorted(CSR_color_lens, key=lambda x: x[3])
        temp_CSR_color_lens = groupby(sort_CSR_color_lens, lambda x: x[3])
        all_color_lens = {}
        for color, group in temp_CSR_color_lens:
            color_list = list(group)
            total_len = 0
            for lis in color_list:
                total_len += abs(lis[2]-lis[1]+1)
            all_color_lens[color] = total_len

        # 统计祖先染色体的颜色和长度大小
        # specs = species.split('_')
        # key = specs[0]+'_'+specs[1]
        # spec_CSR_lens = CSR_lens[key]
        ancestral_chromosome_conf = []
        ancestral_color_lens, ancestral_order = {}, 1
        # sort_CSR_color = sorted(list(define_color.items()), key=lambda x: x[1])
        # temp_list = groupby(sort_CSR_color, lambda x: x[1])
        # for color, group in temp_list:
        for color in order_color:
            # ancestrals = list(group)
            # for anc in ancestrals:
                # ancestral_color_lens[ancestral_order] = ancestral_color_lens.get(ancestral_order, 0)+all_color_lens[anc[1]]
            # ancestral_color_lens[ancestral_order] = [1, ancestral_color_lens[ancestral_order], color, alpla]
            # ancestral_color_lens[ancestral_order] = [1, ancestral_color_lens[ancestral_order], color]

            ancestral_chromosome_conf.append([ancestral_order, 1, all_color_lens[color], color])
            ancestral_order += 1
        # max_lens = max([x[1] for x in ancestral_color_lens.values()]) 
        # for c in ancestral_color_lens:
            # new_row = [c]+ancestral_color_lens[c]
            # ancestral_chromosome_conf.append(new_row)

        # print(ancestral_chromosome_conf)
        return ancestral_chromosome_conf

    def product_ancestral_chromosome1(self, species, define_color, CSR_lens, alpla):
        # 统计祖先染色体的颜色和长度大小
        specs = species.split('_')
        key = specs[0]+'_'+specs[1]
        spec_CSR_lens = CSR_lens[key]
        ancestral_chromosome_conf = []
        ancestral_color_lens, ancestral_order = {}, 1
        sort_CSR_color = sorted(list(define_color.items()), key=lambda x: x[1])
        temp_list = groupby(sort_CSR_color, lambda x: x[1])
        for color, group in temp_list:
            ancestrals = list(group)
            for anc in ancestrals:
                ancestral_color_lens[ancestral_order] = ancestral_color_lens.get(ancestral_order, 0)+spec_CSR_lens[anc[0]]
            # ancestral_color_lens[ancestral_order] = [1, ancestral_color_lens[ancestral_order], color, alpla]
            ancestral_color_lens[ancestral_order] = [1, ancestral_color_lens[ancestral_order], color]
            ancestral_order += 1
        # max_lens = max([x[1] for x in ancestral_color_lens.values()]) 
        for c in ancestral_color_lens:
            new_row = [c]+ancestral_color_lens[c]
            ancestral_chromosome_conf.append(new_row)

        # print(ancestral_chromosome_conf)
        return ancestral_chromosome_conf

    def product_wgd_chromosome(self, species, target_spec, define_color, CSR_each_lens, wgds_chromosome, wgds, alpla):
        
        wgds_names = wgds[target_spec].split('_')
        wgds_chrs_color, wgds_max_lens = [], []
        specs = species.split('_')
        key = specs[0]+'_'+target_spec
        spec_each_lens = CSR_each_lens[key]
        for wgd_chrs in wgds_chromosome:
            midd_chromosome_conf = []
            midd_chr_color, chr_order = {c:[] for c in range(1, len(wgd_chrs)+1)}, 1
            for chr_li in wgd_chrs:
                for CSR in chr_li[4]:
                    CSR_lens = spec_each_lens[chr_li[1]][CSR]
                    color = define_color[CSR]
                    midd_chr_color[chr_order].append([CSR_lens, color])
                chr_order += 1

            for c in midd_chr_color.keys():
                order = 1
                md_rows = midd_chr_color[c]
                # chr_lens = sum([num[0] for num in md_rows])
                break_poins = '_'.join([str(num[0]) for num in md_rows])
                for row in md_rows:
                    # conf_line = [c, order, order+row[0], row[1], alpla]
                    conf_line = [c, order, order+row[0], row[1]]
                    midd_chromosome_conf.append(conf_line)
                    order += row[0]

            wgds_chrs_color.append(midd_chromosome_conf)


        return wgds_names, wgds_chrs_color
        # print(spec_each_lens)

    def map_spec_modeule(self, spec_names, mdinfos):
        # 功能：
        """
        相对与参考物种的一条染色体选择颜色, 最大覆盖度相对应的module 进行相同命名, 之后对应到compair物种的每一条染色体
        """

        # 进一步优化要实现：
        """
        然后统计每一条染色体内的所有module(可以统计类型，或者module的范围), 之后取哪两条染色体交集最大
        物种间的交集是推断祖先(会有合并的情况，也就是一个物种的一条染色体，对应另外物种的多条染色体)
        物种内的是推断加倍前后
        """

        # 取出一个少多倍化的为参考, 并为module定义类型名称 
        modules_colorinfo, modules_infos = {}, {}
        specs = spec_names.split('_')
        ref_spec = specs[0]

        # 为每个物种内的module 映射名称 以及染色(根据参考物种染色体定义), 计算module占染色体百分比
        for spec in specs[1:]:
            key = ref_spec+'_'+spec
            md_info = mdinfos[key]
            md_color_pos, new_md_info = [], []
            # sf = open(f'{key}.top.color.pos.txt', 'w')
            for lis in md_info:
                module_chr_ratio = 0
                color = self.colors[str(lis[0])]
                md_ty = 'CSR'+lis[11]
                top_row = lis[:2]+lis[4:6]+[color, md_ty, lis[12]]
                # sf.write('\t'.join([str(x) for x in top_row])+'\n')
                new_row = lis[:6]+[color, md_ty, lis[12]]
                md_color_pos.append(new_row)

                new_lis = lis+[color, md_ty]
                new_md_info.append(new_lis)
            modules_infos[key] = new_md_info
            modules_colorinfo[key] = md_color_pos

        return modules_infos, modules_colorinfo

    def infer_ancestral_karyotype(self, species, intergenomicratio, chr_types, color_pos):
        sf_car_chromosome = open(f'{self.save_path}/A-{self.species}-ancestral_chromosome_CAR.txt', 'w')

        specs = species.split('_')
        define_color_module, ancestral_chr = self.get_adjacent_between_module(species, intergenomicratio, copy.deepcopy(chr_types), color_pos)
        ancestral_chr = [[str(x[0]), x[1][0], x[2]] for x in ancestral_chr]
        
        # ancestral_chr
        all_CSR, new_ancestral_chr = [], []
        ancestral_chr = sorted(ancestral_chr, key=lambda x: len(x[2].split('_')), reverse=True)
        for lis in ancestral_chr:
            csrs = lis[2].split('_')
            if len(list(set(csrs) & set(all_CSR))) == 0:
                new_ancestral_chr.append(lis)
            all_CSR.extend(csrs)

        new_ancestral_chr = sorted(new_ancestral_chr, key=lambda x: int(x[0]))
        # group_ratio = groupby(new_ratio, lambda x:int(x[0]))


        print(f'chromosome\t{specs[1]}_chr\t{specs[2]}_chr\t\033[1;32mCAR\033[0m')
        sf_car_chromosome.write(f'chromosome\t{specs[1]}_chr\t{specs[2]}_chr\tCAR\n')
        print()
        chrnum = 1
        color_order = {}
        for lis in new_ancestral_chr:
            CSR = lis[2].split('_')[0]
            color = define_color_module[CSR]
            color_order[color] = chrnum
            print(f'chromosome:{chrnum}\t{lis[0]}\t{lis[1]}\t\033[1;32m{lis[2]}\033[0m')
            sf_car_chromosome.write(f'chromosome:{chrnum}\t{lis[0]}\t{lis[1]}\t{lis[2]}\n')
            chrnum += 1

        sf_car_chromosome.close()

        return define_color_module, color_order

    def infer_ancestral_karyotype_common(self, species, conserved_spec, wgds, chr_types, recent_wgdchr, flag=False):
        sf_car_chromosome = open(f'{self.save_path}/A-{self.species}-ancestral_chromosome_CAR.txt', 'w')
        specs = species.split('_')
        all_wgds_chromosome = self.get_adjacent_within_module(species, conserved_spec, wgds, copy.deepcopy(chr_types), recent_wgdchr, flag)
        ancestral_chr = all_wgds_chromosome[-1]
        # print("{:<10}{:<15}{:<15}{:<15}".format("chromosome",f"{specs[1]}_chr", f"{specs[1]}_chr", "\033[1;32mCAR\033[0m"))

        all_CSR, new_ancestral_chr = [], []
        ancestral_chr = sorted(ancestral_chr, key=lambda x: len(x[4]), reverse=True)
        for lis in ancestral_chr:
            csrs = lis[4]
            if len(list(set(csrs) & set(all_CSR))) == 0:
                new_ancestral_chr.append(lis)
            all_CSR.extend(csrs)

        new_ancestral_chr = sorted(new_ancestral_chr, key=lambda x: int(x[0]))

        print(f'chromosome\t{conserved_spec}_chr\t{conserved_spec}_chr\t\033[1;32mCAR\033[0m')
        sf_car_chromosome.write(f'chromosome\t{specs[1]}_chr\t{specs[1]}_chr\tCAR\n')
        print()
        chrnum = 1
        CSR_order = {}
        for lis in new_ancestral_chr:
            CSR = lis[4][0]
            CSR_order[CSR] = chrnum
            # color = define_color_module[CSR]
            # color_order[color] = chrnum
            link_model = '_'.join(lis[4])
            # print("chromosome:{}{:^15}{:^15}{:<15}".format(chrnum, lis[0], lis[1], f"\033[1;32m{link_model}\033[0m"))
            print(f'chromosome:{chrnum}\t{lis[0]}\t{lis[1]}\t\033[1;32m{link_model}\033[0m')
            sf_car_chromosome.write(f'chromosome:{chrnum}\t{lis[0]}\t{lis[1]}\t{link_model}\n')
            chrnum += 1

        sf_car_chromosome.close()

        return all_wgds_chromosome, CSR_order

    def print_ancestral_before_after(self, species, wgd_num, all_wgds_chromosome):
        sf_car_chromosome = open(f'{self.save_path}/A-{species}-chromosome_CAR.WGD-before.txt', 'w')
        specs = species.split('_')
        wgds = wgd_num[specs[1]].split('_')

        print('\n')
        print(f'chromosome\t{specs[1]}_chr\t{specs[1]}_chr\t\033[1;32mCAR\033[0m')
        sf_car_chromosome.write(f'chromosome\t{specs[1]}_chr\t{specs[1]}_chr\tCAR\n')
        print()
        # CSR_order = {}
        for ancestral_chr in all_wgds_chromosome:
            name = wgds.pop()
            sf_car_chromosome.write('\n')
            print(f'WGD {name} Before duplication: {len(ancestral_chr)} chromosome')
            sf_car_chromosome.write(f'WGD {name} Before duplication: {len(ancestral_chr)} chromosome\n')
            chrnum = 1
            for lis in ancestral_chr:
                link_model = '_'.join(lis[4])
                print(f'chromosome:{chrnum}\t{lis[0]}\t{lis[1]}\t\033[1;32m{link_model}\033[0m')
                sf_car_chromosome.write(f'chromosome:{chrnum}\t{lis[0]}\t{lis[1]}\t{link_model}\n')
                chrnum += 1

        sf_car_chromosome.close()

        return all_wgds_chromosome

    def get_adjacent_between_module(self, species, sub_num, chr_all_types, color_pos):
        # 判断
        """
        # 判断哪些是同一个CSR  之后判断哪些module是相邻的
        # 判断完之后 判断物种module的连接顺序
        """

        """
        先物种间推断，之后发生多次加倍的 再推断物种内

        """
        
        define_type_module, colors = {}, list(self.colors.values())

        specs = species.split('_')
        ref_spec = specs[0]
        init_key = ref_spec+'_'+specs[1]
        init_types = chr_all_types[init_key]  # 每个物种内 每条染色体包含的CSR

        iteration_num, ancestral_chr = 0, []
        # Complete ancestral chromosome
        while iteration_num < 10:  # 10000
            midd = {int(i):[] for i in init_types.keys()}
            ke = ref_spec+'_'+specs[2]
            com_types = chr_all_types[ke]
            combin_chrs = list(product(init_types, com_types))
            ratio = []
            for chr_combin in combin_chrs:
                ref_chr, com_chr = chr_combin
                ref = init_types[ref_chr]
                com = com_types[com_chr]
                ref_mds = set([x[3] for x in ref])
                com_mds = set([x[3] for x in com])
                comm_ty = list(ref_mds & com_mds)
                remain_ref = ref_mds-set(comm_ty)
                remain_com = com_mds-set(comm_ty)
                ref_comm_ratio = len(comm_ty)/len(ref_mds)
                com_comm_ratio = len(comm_ty)/len(com_mds)
                ratio.append([ref_chr, com_chr, ref_comm_ratio, com_comm_ratio, comm_ty, list(remain_ref), list(remain_com)])
            new_ratio = sorted(ratio, key=lambda x: int(x[0]))
            group_ratio = groupby(new_ratio, lambda x:int(x[0]))

            for name, group in group_ratio:
                new_group = sorted(list(group), key=lambda x:(float(x[2]), float(x[3])), reverse=True)
                pair_chr = new_group[:sub_num]
                result_gr = [x for x in pair_chr if x[2] > 0 and x[3] > 0]
                # print([x[4] for x in pair_chr if x[2] > 0 and x[3] > 0])
                midd[name] = result_gr

            # print(midd)
            # 判断是不是所有物种里面相同的结构  # 一般是两个物种 以及一个参考物种
            for chr_key in midd:
                midd_max_ty, link_dit, CSRs = '', {}, []
                lines = midd[chr_key]
                for x in lines:
                    link_k = '_'.join(sorted(x[4]))
                    if link_k in link_dit:
                        link_dit[link_k].append(x[1])
                    else:
                        link_dit[link_k] = [x[1]]
                    CSRs.append(link_k)
                if CSRs:
                    CSRs_count = Counter(CSRs).items()
                    a, b = [], []
                    for k, v in CSRs_count:
                        a.append(k)
                        b.append(v)
                    sum_b = sum(b)
                    b = [x/sum_b for x in b]

                    max_ty = Counter(CSRs).most_common(1)[0]
                else:
                    continue
                # 一般至少有一半保留相同的结构
                if max_ty[1] >= sub_num/2:
                    midd_max_ty = max_ty[0]

                    # 属于祖先结构
                    ancestral_chr.append([chr_key, link_dit[midd_max_ty], midd_max_ty])

                    # 删除第一轮最大交集的 并进行颜色定义
                    tyss = midd_max_ty.split('_')
                    color, flag = colors[color_pos], 1
                    if any(([(t in define_type_module) for t in tyss])):
                        # 如果某个ty之前定义了颜色, 那么这个最大的交集的所有ty都使用定义的颜色
                        ke = list(set(tyss) & set(define_type_module.keys()))[0]
                        color = define_type_module[ke]
                        flag = 0
                    for tty in tyss:
                        define_type_module[tty] = color
                    if flag:
                        color_pos += 1

                    for spec in specs[1:]:
                        key = specs[0]+'_'+spec
                        del_chr = []
                        # 将所选物种中 每条染色体上包含祖先节点的ty 都删掉
                        for k in chr_all_types[key].keys():
                            new_tt = []
                            for t in chr_all_types[key][k]:
                                if not t[3] in tyss:
                                    new_tt.append(t)
                            if new_tt:
                                chr_all_types[key][k] = new_tt
                            else:
                                del_chr.append(k)
                        for c in del_chr:
                            chr_all_types[key].pop(c)
            iteration_num += 1

        # 打印
        # print(f'chromosome order\tref_chr\tcom_chr\tConnection pattern')
        # chrnum = 1
        # for lis in ancestral_chr:
        #     print(f'chromosome:{chrnum}\t{lis[0]}\t{lis[1][0]}\t{lis[2]}')
        #     chrnum += 1

        return define_type_module, ancestral_chr

    def get_adjacent_within_module(self, sepcies, target_spec, wgd_num, chr_all_types, define_recent_wgd_file, flagpr=True):
        """
        推断物种内的话,可以自己输入祖先核型, 然后根据参考物种的核型来进一步
        推断物种内推断加倍前后, 参考图的那种情况
        """
        # 读取定义的最近加倍文件
        define_recent_chr = {}
        for li in open(define_recent_wgd_file):
            lis = li.strip().split()
            if lis[0] in define_recent_chr:
                define_recent_chr[lis[0]].append(lis[1])
            else:
                define_recent_chr[lis[0]] = [lis[1]]


        all_wgds_chromosome = []
        sub_num = {'2': 1, '3': 2, '5': 4, '4':3, '8':7}  #发生一次乘几的多倍化
        wgds = wgd_num[target_spec].split('_')
        specs = sepcies.split('_')
        ref_spec = specs[0]
        key = ref_spec+"_"+target_spec
        com_types = chr_all_types[key]
        old_types = copy.deepcopy(com_types)

        # 最近加倍产生
        na = wgds.pop(-1)
        wgd_recent_order = []
        for k, values in define_recent_chr.items():
            ref = com_types[k]
            for v in values:
                com = com_types[v]
                ref_mds = set([x[3] for x in ref])
                com_mds = set([x[3] for x in com])
                comm_ty = list(ref_mds & com_mds)
                remain_ref = ref_mds-set(comm_ty)
                remain_com = com_mds-set(comm_ty)
                ref_comm_ratio = len(comm_ty)/len(ref_mds)
                com_comm_ratio = len(comm_ty)/len(com_mds)
                if ref_comm_ratio < 0.25 or com_comm_ratio < 0.25:
                    continue
                wgd_recent_order.append([k, v, ref_comm_ratio, com_comm_ratio, comm_ty])
        all_wgds_chromosome.append(wgd_recent_order)
        # print(f'chromosome order\tref_chr\tcom_chr\tConnection pattern')
        # chrnum = 1
        # for lis in wgd_recent_order:
            # link_model = '_'.join(lis[4])
            # print(f'chromosome:{chrnum}\t{lis[0]}\t{lis[1]}\t{link_model}')
            # chrnum += 1
        # if flagpr:
            # print(f'WGD {na} Before duplication: {len(wgd_recent_order)} chromosome')
        deal_chr = sorted(list(set([x[0] for x in wgd_recent_order])), key=lambda x: int(x))

        for dc in deal_chr:
            if dc in com_types:
                com_types.pop(dc)
                
        while wgds:
            """
            注意 每次加倍循环到底的时候 都要深拷贝上一次加倍剩下的染色体
            # 先把最大交集的取完,然后把里面最大交集的删除 再进行下一轮交集  直到这次加倍所有的块取完
            # 取完之后 或者有一些易位之类的（设置条件忽略  大部分都删除了，只有这个没有 那么就可以删了）
            # 那么每次循环都要备份
            """

            # 多次多倍化受到同源的影响 比如4套  那么1套可能找到3套都比较好
            ###### 推断第一次加倍的染色体
            wgd = wgds.pop(-1)
            wgd_chr, wgd_color_order = [], []
            num = 0
            # print('-------------------------------------------------')
            # print('start')
            midd_com_types = copy.deepcopy(com_types)
            while num <10:
                combin_chrs = list(product(midd_com_types.keys(), repeat=2))
                ratio = []
                for chr_combin in combin_chrs:
                    ref_chr, com_chr = chr_combin
                    # # 自己和自己的不确定
                    if ref_chr == com_chr:
                        continue
                    ref = midd_com_types[ref_chr]
                    com = midd_com_types[com_chr]
                    # if (not ref) or (not com):
                        # continue
                    ref_mds = set([x[3] for x in ref])
                    com_mds = set([x[3] for x in com])
                    comm_ty = list(ref_mds & com_mds)
                    remain_ref = ref_mds-set(comm_ty)
                    remain_com = com_mds-set(comm_ty)
                    cover_ref = min(sum([x[4] for x in ref if x[3] in comm_ty]), 1)
                    cover_com = min(sum([x[4] for x in com if x[3] in comm_ty]), 1)
                    ref_comm_ratio = len(comm_ty)/len(ref_mds)
                    com_comm_ratio = len(comm_ty)/len(com_mds)

                    ratio.append([ref_chr, com_chr, ref_comm_ratio, com_comm_ratio, cover_ref, cover_com, comm_ty, list(remain_ref), list(remain_com)])
                new_ratio = sorted(ratio, key=lambda x: int(x[0]))
                group_ratio = groupby(new_ratio, lambda x:int(x[0]))

                midd = {}
                for name, group in group_ratio:
                    # new_group = sorted(list(group), key=lambda x:(float(x[2]), float(x[3])), reverse=True)
                    new_group = sorted(list(group), key=lambda x:(float(x[2]), float(x[3])), reverse=True)
                    new_group = [x for x in new_group if x[2] > 0 and x[3] > 0]
                    # print(new_group)
                    pair_chr = new_group[:sub_num[wgd]]
                    # result_gr = [x for x in pair_chr if x[2] > 0 and x[3] > 0]
                    # print([x[4] for x in pair_chr if x[2] > 0 and x[3] > 0])
                    # print(result_gr)

                    midd[name] = pair_chr
                # print('-------------------------------')

                # 先把最大交集的取完,然后把里面最大交集的删除 再进行下一轮交集  直到这次加倍所有的块取完
                # 取完之后 或者有一些易位之类的（设置条件忽略  大部分都删除了，只有这个没有 那么就可以删了）
                # 那么每次循环都要备份

                midd_wgd_color_order = []
                for chr_key in midd:
                    lines = midd[chr_key]

                    for lis in lines:
                        tyss = lines[0][6]

                        new_order = []
                        ty_order = [x[3] for x in old_types[str(chr_key)]]
                        new_order = [[ty_order.index(x), x] for x in tyss]
                        new_order = sorted(new_order, key= lambda x: x[0])
                        new_order = [x[1] for x in new_order]

                        # 需要记录已经删除的 并且这些没有再记录里 才能记录
                        midd_wgd_color_order.append([str(chr_key), lis[1], lis[2], lis[3], new_order])


                # 对列表中的列表进行去重
                flag, new_midd_wgd_color_order = [], []
                for lis in midd_wgd_color_order:
                    lkk = f'{lis[0]}_{lis[1]}'
                    lkk1 = f'{lis[1]}_{lis[0]}'
                    if lkk not in flag and lkk1 not in flag:
                        new_midd_wgd_color_order.append(lis)
                        flag.append(lkk)
                        flag.append(lkk1)

                sort_midd_wgd_color_order = sorted(new_midd_wgd_color_order, key=lambda x: (len(x[4]), float(x[2]), float(x[3])), reverse=True)


                # 把涉及这一连接模式的都删去
                new_midd_wgd_color_order = []
                for lk in sort_midd_wgd_color_order:
                    if lk[0] not in midd_com_types or lk[1] not in midd_com_types:
                        continue
                    onetys = [t[3] for t in midd_com_types[lk[0]]]
                    twotys = [t[3] for t in midd_com_types[lk[1]]]
                    fg1 = [1 if t in onetys else 0 for t in lk[4]]
                    fg2 = [1 if t in twotys else 0 for t in lk[4]]

                    if all((fg1)) and all((fg2)):
                        wgd_color_order.append(lk)
                        del_chr = []
                        # 将涉及这两条染色体的部分删除
                        for k in lk[:2]:
                            remove_tt, new_tt = [], []
                            
                            # 记录坐标 只删除第一个
                            old_lks = midd_com_types[k]
                            old_lk = [x[3] for x in old_lks]
                            for t in lk[4]:
                                if t in old_lk:
                                    pos = old_lk.index(t)
                                    remove_tt.append(pos)

                            for i, v in enumerate(old_lks):
                                if i in remove_tt:
                                    continue
                                new_tt.append(v)
                            if new_tt:
                                midd_com_types[k] = new_tt
                            else:
                                del_chr.append(k)
                        for c in del_chr:
                            midd_com_types.pop(c)

                num+=1


            num_chrs = [c for x in wgd_color_order for c in x[:2]]
            dit_num = Counter(num_chrs)

            result_wgd_color_order = []
            midd_left, midd_right, eq_midd = [], [], []
            for i in range(len(wgd_color_order)):
                row = wgd_color_order[i]
                if dit_num[row[0]] < dit_num[row[1]]:
                    midd_left.append(row[1])
                    midd_right.append(row[0])
                    new_row = [row[1], row[0], row[3], row[2], row[4]]
                    result_wgd_color_order.append(new_row)
                elif dit_num[row[0]] > dit_num[row[1]]:
                    midd_left.append(row[0])
                    midd_right.append(row[1])
                    result_wgd_color_order.append(row)
                else:
                    eq_midd.append(row)
            for row in eq_midd:
                if row[0] in midd_left:
                    result_wgd_color_order.append(row)
                    midd_left.append(row[0])
                    midd_right.append(row[1])
                elif row[1] in midd_left:
                    new_row = [row[1], row[0], row[3], row[2], row[4]]
                    result_wgd_color_order.append(new_row)
                    midd_left.append(row[1])
                    midd_right.append(row[0])
                elif row[0] in midd_right:
                    midd_left.append(row[1])
                    midd_right.append(row[0])
                    new_row = [row[1], row[0], row[3], row[2], row[4]]
                    result_wgd_color_order.append(new_row)
                elif row[1] in midd_right:
                    midd_left.append(row[0])
                    midd_right.append(row[1])
                    result_wgd_color_order.append(row)
                else:
                    result_wgd_color_order.append(row)

            all_wgds_chromosome.append(result_wgd_color_order)

            # chrnum1 = 1
            # for lis in result_wgd_color_order:
            #     link_model = '_'.join(lis[4])
            #     print(f'chromosome:{chrnum1}\t{lis[0]}\t{lis[1]}\t{link_model}')
            #     chrnum1 += 1
            # if flagpr:
            #     print(f'WGD {wgd} Before duplication: {len(wgd_color_order)} chromosome')

            # print('-------------------------------------------------')
            deal_chr = [x[0] for x in result_wgd_color_order]
            for dc in deal_chr:
            # for dc in ['2', '10', '19', '5', '26', '3']:
                if dc in com_types:
                    com_types.pop(dc)


            # 之后进行下一次循环
        return all_wgds_chromosome

    def define_ancestral_color(self, all_wgds_chromosome, color_pos, chr_spec_types):
        """
        加倍前后生成祖先染色体颜色
        这时候是不是每条有相同的CSR名称，这时候考虑二维坐标系

        如果一个染色体内存在相同名称CSR
        """
        define_ty_color, colors = {}, list(self.colors.values())
        ancestral = all_wgds_chromosome.pop()
        for chromosome in ancestral:
            color = colors[color_pos]
            for csr in chromosome[4]:
                # key = f'{csr}_{chromosome[0]}_{chromosome[1]}'
                # define_ty_color[key] = color
                define_ty_color[csr] = color
            color_pos += 1

        all_CSR = list(define_ty_color.keys())
        for key_chr in chr_spec_types.keys():
            middle_color = []
            chr_all_CSR = [ x[3] for x in chr_spec_types[key_chr]]
            # chr_all_CSR = [ f'{x[3]}_{x[4]}_{key_chr}' for x in chr_spec_types[key_chr]]

            nomatch_CSR = list(set(chr_all_CSR)-set(all_CSR))
            match_CSR = list(set(chr_all_CSR)-set(nomatch_CSR))
            if nomatch_CSR:
                middle_color = [define_ty_color[x] for x in match_CSR]
                max_cor = Counter(middle_color).most_common(1)[0]
                for Csr in nomatch_CSR:
                    define_ty_color[Csr] = max_cor[0]

        return define_ty_color

    def count_csr_onchromosome(self, modules_infos, lens_dit):
        """
        首先判断比较物种中有一个物种的染色体满足对应的类型(颜色)为整条染色体，那么这个为祖先染色体，
        (也需要结合加倍，如果有一次二倍化)  还需要进一步讨论（满足于有一条完整还是每个物种都需要有）
        
        可以计算一下module 占compair物种整条染色体的比例

        然后统计每一条染色体内的所有module(可以统计类型，或者module的范围), 之后取哪两条染色体交集最大
        物种间的交集是推断祖先(会有合并的情况，也就是一个物种的一条染色体，对应另外物种的多条染色体)
        """

        # 将没有在两个物种中共同存在的module 去掉
        no_exist_CSR = []
        all_key = list(modules_infos.keys())
        if len(all_key) > 1:  # 如果只有一个物种不用去除
            ll_CSR = [x[14] for x in modules_infos[all_key[0]]]
            yy_CSR = [x[14] for x in modules_infos[all_key[1]]]
            no_exist_CSR.extend(list(set(ll_CSR)-set(yy_CSR)))
            no_exist_CSR.extend(list(set(yy_CSR)-set(ll_CSR)))

        # 统计每一条染色体内的所有module 类型以及所占染色体的覆盖度
        temp_list,chr_all_types = [], {}
        for key in modules_infos.keys():
            spec_md_infos = modules_infos[key]
            spec_md_infos = sorted(spec_md_infos, key=lambda x: x[1])
            temp_list = groupby(spec_md_infos, lambda x: x[1])
            midd_chr_all_type, midd_chr_all_ls_type = {}, {}
            for name, group in temp_list:
                group_list = list(group)
                chr_all_type = []
                for lis in group_list:
                    if lis[-1] in no_exist_CSR:
                        continue
                    specs = key.split('_')
                    chr_len = lens_dit[specs[1]][int(lis[1])]
                    # coverage = abs(int(lis[5])-int(lis[4]))/chr_len
                    coverage = abs(int(lis[7]))/chr_len

                    # module com_start com_end CSRnum ref_chr
                    chr_all_type.append([lis[12], lis[4], lis[5], lis[-1], coverage])
                # 将每个染色体的module 按顺序排列
                chr_all_type = sorted(chr_all_type, key=lambda x: int(x[1]))
                midd_chr_all_type[name] = chr_all_type
            chr_all_types[key] = midd_chr_all_type

        return chr_all_types

    def read_blockinfo(self, block_info_file, bk_num, hocv, hocv_depth):
        new_block, hocv_pos = [], 10
        hocv_pos += hocv_depth
        block_list = [x.strip().split(',') for x in open(block_info_file)]

        for bk in block_list[1:]:
            if not int(bk[8]) >= bk_num:
                continue
            if float(bk[hocv_pos]) >= hocv:
                new_block.append(bk)

        return new_block

    def read_spec_all_module(self, block_list):
        block_info = []

        ##########
        # 根据module 的序号进行分组
        temp_list = []
        block_list = sorted(block_list, key=lambda x: x[25])
        temp_list = groupby(block_list, lambda x: x[25])

        # 生成module时记录下具体的block位置，因为block之间有空缺，在算覆盖度时不准确
        for name, group in temp_list:
            group_list = list(group)
            chr1 = list(set([x[1] for x in group_list]))[0]
            chr2 = list(set([x[2] for x in group_list]))[0]
            top_lens = sum([abs(int(x[5]) - int(x[6]))+1 for x in group_list])
            left_lens = sum([abs(int(x[3]) - int(x[4]))+1 for x in group_list])
            start1, end1 = min([x[3] for x in group_list]+[x[4] for x in group_list], key=lambda x: int(x)), max([x[3] for x in group_list]+[x[4] for x in group_list], key=lambda x: int(x))
            start2, end2 = min([x[5] for x in group_list]+[x[6] for x in group_list], key=lambda x: int(x)), max([x[5] for x in group_list]+[x[6] for x in group_list], key=lambda x: int(x))
            all_ks = [float(v) for x in group_list for v in x[21].split('_')]
            new_all_ks = list(filter(lambda x: x >= 0, all_ks))
            ks_median, ks_average = np.median(new_all_ks), np.mean(new_all_ks)
            length_block = sum([int(x[8]) for x in group_list])

            CSR, module = list(set([x[24] for x in group_list]))[0], list(set([x[25] for x in group_list]))[0]
            Ts = list(set([x[26] for x in group_list]))[0]
            if Ts == '1':
                continue
            block_info.append(
                [chr1, chr2, start1, end1, start2, end2, left_lens, top_lens, length_block, ks_median, ks_average,
                 CSR, module])  # 11 12

        return block_info

    def get_new_color_list(self, color_list, len_dit):
        # 根据block 对色块进行补全
        new_color_list = []

        # 去重
        lk_color = ['_'.join([str(x) for x in lis]) for lis in color_list]
        lk_color = list(set(lk_color))
        color_list = []
        for x in lk_color:
            a = x.split('_')
            color_list.append([a[0], int(a[1]), int(a[2]), a[3]])

        
        color_order_list = sorted(color_list, key=lambda x: x[0])
        temp_list = groupby(color_order_list, lambda x: x[0])
        for name, group in temp_list:

            group_list = list(group)
            group_order_list = sorted(group_list, key=lambda x: (int(x[1]), x[3]))

            chr_len = len_dit[int(name)]
            color_old, start_old, end_old = '', 0, 0

            midd_color_list = []
            for index, row in enumerate(group_order_list):
                if index == 0:
                    start_old = min(1, row[1])
                    end_old = row[2]
                    color_old = row[3]
                    continue
                if row[3] == color_old:
                    end_old = max(row[2], end_old)
                    continue
                else:
                    last_row = group_order_list[index-1]
                    # 可能这个末尾 也小于end_old
                    if end_old >= row[2]:
                        midd_color_list.append([name, start_old, row[1]-1, color_old])
                        midd_color_list.append([name, row[1], row[2], row[3]])
                        start_old = row[2]+1
                        continue
                    overlap = row[1]-end_old
                    if overlap < 0:
                        t = row[1]
                        row[1] = end_old
                        end_old = t
                    overlap = abs(overlap)
                    if overlap == 1:
                        midd_color_list.append([name, start_old, end_old, color_old])
                        start_old = row[1]
                        end_old = row[2]
                        color_old = row[3]
                    elif overlap >= 2:
                        last_len = abs(last_row[2]-last_row[1])+1
                        now_len = abs(row[2]-row[1])+1
                        last_add = round(last_len/(last_len+now_len)*(overlap-1))
                        now_dele = round(now_len/(last_len+now_len)*(overlap-1))
                        if last_add == 0 and now_dele == 0:
                            last_add += 1
                            # now_dele-=1
                        midd_color_list.append([name, start_old, end_old+last_add, color_old])
                        start_old = row[1]-now_dele
                        end_old = row[2]
                        color_old = row[3]
                    elif overlap == 0:
                        midd_color_list.append([name, start_old, end_old, color_old])
                        start_old = row[1]+1
                        end_old = row[2]
                        color_old = row[3]
            midd_color_list.append([name, start_old, end_old, color_old])
            last_row = midd_color_list[-1]
            max_end = max(chr_len, last_row[2])
            midd_color_list[-1][2] = max_end

            new_color_list.extend(midd_color_list)

        return new_color_list

    def get_new_color_list1(self, color_list, len_dit):
        # 根据block 对色块进行补全
        new_color_list = []

        # 去重
        lk_color = ['_'.join([str(x) for x in lis]) for lis in color_list]
        lk_color = list(set(lk_color))
        color_list = []
        for x in lk_color:
            a = x.split('_')
            color_list.append([a[0], int(a[1]), int(a[2]), a[3]])

        
        color_order_list = sorted(color_list, key=lambda x: x[0])
        temp_list = groupby(color_order_list, lambda x: x[0])
        for name, group in temp_list:

            group_list = list(group)
            group_order_list = sorted(group_list, key=lambda x: int(x[1]))

            chr_len = len_dit[int(name)]
            color_old, start_old, end_old = '', 0, 0

            midd_color_list = []
            for index, row in enumerate(group_order_list):
                if index == 0:
                    start_old = min(1, row[1])
                    end_old = row[2]
                    color_old = row[3]
                    continue
                if row[3] == color_old:
                    end_old = max(row[2], end_old)
                    continue
                else:
                    last_row = group_order_list[index-1]
                    # 可能这个末尾 也小于end_old
                    if end_old >= row[2]:
                        midd_color_list.append([name, start_old, row[1]-1, color_old])
                        midd_color_list.append([name, row[1], row[2], row[3]])
                        start_old = row[2]+1
                        continue
                    overlap = row[1]-end_old
                    if overlap < 0:
                        t = row[1]
                        row[1] = end_old
                        end_old = t
                    overlap = abs(overlap)
                    if overlap == 1:
                        midd_color_list.append([name, start_old, end_old, color_old])
                        start_old = row[1]
                        end_old = row[2]
                        color_old = row[3]
                    elif overlap >= 2:
                        last_len = abs(last_row[2]-last_row[1])+1
                        now_len = abs(row[2]-row[1])+1
                        last_add = round(last_len/(last_len+now_len)*(overlap-1))
                        now_dele = round(now_len/(last_len+now_len)*(overlap-1))
                        if last_add == 0 and now_dele == 0:
                            last_add += 1
                            # now_dele-=1
                        midd_color_list.append([name, start_old, end_old+last_add, color_old])
                        start_old = row[1]-now_dele
                        end_old = row[2]
                        color_old = row[3]
                    elif overlap == 0:
                        midd_color_list.append([name, start_old, end_old, color_old])
                        start_old = row[1]+1
                        end_old = row[2]
                        color_old = row[3]
            midd_color_list.append([name, start_old, end_old, color_old])
            last_row = midd_color_list[-1]
            max_end = max(chr_len, last_row[2])
            midd_color_list[-1][2] = max_end

            new_color_list.extend(midd_color_list)

        return new_color_list

    def get_module_block_color(self, block_list, ancestor_color, lens_dit, name):
        spec1, spec2 = name.split('_')
        top_color_list, left_color_list = [], []
        # Extraction of blocks
        for lis in block_list:
            chr1, chr2 = lis[1], lis[2]
            start1, end1 = min([int(lis[3]), int(lis[4])]), max([int(lis[3]), int(lis[4])])
            start2, end2 = min([int(lis[5]), int(lis[6])]), max([int(lis[5]), int(lis[6])])
            ty_name = 'CSR'+lis[24]
            color = ancestor_color[ty_name]
            module_order = lis[25]
            top_color_list.append([chr2, int(start2), int(end2), color])
            left_color_list.append([chr1, int(start1), int(end1), color])
        # Expansion of block color range
        new_top_color_list = self.get_new_color_list(top_color_list, lens_dit[spec2])
        new_left_color_list = self.get_new_color_list(left_color_list, lens_dit[spec1])
        new_top_color_list = sorted(new_top_color_list, key=lambda x: int(x[0]))
        new_left_color_list = sorted(new_left_color_list, key=lambda x: int(x[0]))

        return new_top_color_list, new_left_color_list

    def get_define_ancestral_color(self, ancestral_color_file):
        define_ty_color = {}
        for li in open(ancestral_color_file):
            lis = li.strip().split()
            define_ty_color[lis[0]] = lis[1]
        return define_ty_color

    def run(self):
        specs = self.species.split('_')
        mdinfos, block_list_dit = {}, {}
        for key, bkfile in self.bk_files.items():
            block_list = self.read_blockinfo(bkfile, self.block_num, self.hocv, self.hocv_depths[key])
            md_info = self.read_spec_all_module(block_list)
            mdinfos[key] = md_info
            block_list_dit[key] = block_list

        spec_lens_dit, lens_list = self.read_lens(self.len_files)

        # Add CSR names to each module
        new_mdinfos, modules_colorinfo = self.map_spec_modeule(self.species, mdinfos)

        # Counting CSRs on each chromosome
        chr_all_types = self.count_csr_onchromosome(new_mdinfos, spec_lens_dit)

        # The total length of identical CSRs in all chromosomes of each species was counted.
        # For inferring the proportion of ancestral chromosomes

        # The proportion of identical CSRs in each chromosome in each species was counted and
        # used to infer the proportion of chromosomes before and after WGD
        CSR_top_all_lens, CSR_top_each_len = self.count_CSR_lens(new_mdinfos)

        # commond polyploidization to infer ancestral karyotypes
        print()
        print()
        # print(f'Get the most likely ancestral karyotype \033[1;31m......\033[0m')
        print(f'Get the most likely ancestral karyotype......')
        print()
        print()
        if self.common_wgd.upper() == 'TRUE':
            kkey = f'{specs[0]}_{self.conserved_spec}'
            all_wgds_chromosome, CSR_order = self.infer_ancestral_karyotype_common(self.species, self.conserved_spec, self.wgds, copy.deepcopy(chr_all_types), self.recentwgdchr[self.conserved_spec])
            # all_wgds_chromosome = self.get_adjacent_within_module(self.species, self.conserved_spec, self.wgds, copy.deepcopy(chr_all_types), self.recentwgdchr[self.conserved_spec], False)  # '2_3' 一次2倍一次三倍
            define_color_module = self.define_ancestral_color(all_wgds_chromosome, self.color_pos, chr_all_types[kkey])
            color_order = {define_color_module[c]:CSR_order[c] for c in CSR_order}
        # Inferring ancestral karyotype
        else:
            define_color_module, color_order = self.infer_ancestral_karyotype(self.species, self.intergenomicratio, copy.deepcopy(chr_all_types), self.color_pos)
            # define_color_module, ancestral_chr = self.get_adjacent_between_module(self.species, self.intergenomicratio, copy.deepcopy(chr_all_types), self.color_pos)


        expansion_color, spec_color = {}, {}
        if len(specs) == 3:
            names = [f'{specs[0]}_{specs[1]}', f'{specs[0]}_{specs[2]}']
        else:
            names = [f'{specs[0]}_{specs[1]}']

        for name in names:
            spec1, spec2 = name.split('_')
            top_color, left_color = self.get_module_block_color(block_list_dit[name], define_color_module, spec_lens_dit, name)
            expansion_color[name] = [top_color, left_color]
            spec_color[spec2] = top_color
            if spec1 in spec_color:
                if len(left_color) < len(spec_color[spec1]):
                    spec_color[spec1] = left_color
            else:
                spec_color[spec1] = left_color

        # ancestral_chromosome_conf = self.product_ancestral_chromosome(self.species, define_color_module, CSR_top_all_lens, 1)
        ancestral_chromosome_conf = self.product_ancestral_chromosome(spec_color[self.select_ref_ancestor_spec], color_order)
        
        import pickle
        with open(f'{self.save_path}/ancestor_color_order', 'wb') as save:
            pickle.dump(color_order, save)


        sf_ags = open(f'{self.save_path}/A.AKRUP-ags-select_{self.select_ref_ancestor_spec}.pep.Construct_ancestral_genomes.conf.txt', 'w')
        new_ags_color = ['\t'.join([str(v) for v in x]) for x in spec_color[self.select_ref_ancestor_spec]]
        sf_ags.write('\n'.join(new_ags_color))

        # save ancestral file
        sf_anc_chromosome = open(f'{self.save_path}/A-select_{self.select_ref_ancestor_spec}-ancestral_chromosome_conf.txt', 'w')
        ancestral_chromosome_conf = sorted(ancestral_chromosome_conf, key=lambda x: int(x[0]))
        new_anc_color = ['\t'.join([str(v) for v in x]) for x in ancestral_chromosome_conf]
        sf_anc_chromosome.write('\n'.join(new_anc_color))

        # Save the ancestral color of the CSR
        sf_anc_CSR_color = open(f'{self.save_path}/A-{self.species}-ancestral_CSR-color_conf.txt', 'w')
        for CSR, cor in define_color_module.items():
            sf_anc_CSR_color.write(f'{CSR}\t{cor}\n')

        
        # import threading
        event = threading.Event()
        print()
        content = 'Drawings are in progress'
        th = threading.Thread(target=await_run, args=(content,event))
        th.start()

        try:
            # for name in self.infer_name:
            # for name in [f'{specs[0]}_{specs[1]}', f'{specs[0]}_{specs[2]}']:
            for name in names:
                spec1, spec2 = name.split('_')

                # Karyotype color block for expansion
                top_color, left_color = expansion_color[name]
                # top_color, left_color = self.get_module_block_color(block_list_dit[name], define_color_module, spec_lens_dit, name)
                
                # save Karyotype color
                sf_top = open(f'{self.save_path}/A.{name}.top.color.pos.txt', 'w')
                sf_left = open(f'{self.save_path}/A.{name}.left.color.pos.txt', 'w')
                new_top_color = ['\t'.join([str(v) for v in x]) for x in top_color]
                new_left_color = ['\t'.join([str(v) for v in x]) for x in left_color]
                sf_top.write('\n'.join(new_top_color))
                sf_left.write('\n'.join(new_left_color))

                if self.infer_wgd_flag.upper() == 'TRUE' and name in self.infer_name:
                    # Inferring karyotypes before and after polyploidization
                    target_spec = spec2
                    all_wgds_chromosome = self.get_adjacent_within_module(self.species, target_spec, self.wgds, copy.deepcopy(chr_all_types), self.recentwgdchr[target_spec])  # '2_3' 一次2倍一次三倍
                    
                    self.print_ancestral_before_after(name, self.wgds, all_wgds_chromosome)
                    wgds_names, wgds_chr_colors = self.product_wgd_chromosome(self.species, target_spec, define_color_module, CSR_top_each_len, all_wgds_chromosome, self.wgds, 1)

                    for wgdna, wgds_anc in zip(wgds_names, wgds_chr_colors):
                        sf_color_chromosome = open(f'{self.save_path}/A-{spec2}-chromosome_ancestral_color_pos.{wgdna}-WGD-before.txt', 'w')
                        for lis in wgds_anc:
                            newlis = [str(x) for x in lis]
                            sf_color_chromosome.write('\t'.join(newlis)+'\n')
                        sf_color_chromosome.close()

                    # Karyotype mapping
                    p = DotplotBlock(name, lens_list, block_list_dit[name], top_color, left_color, self.latin_name, define_color_module, self.save_path)
                    p.main(ancestral_chromosome_conf, wgds_chr_colors)
                elif self.infer_wgd_flag.upper() == 'FALSE':
                    p = DotplotBlock(name, lens_list, block_list_dit[name], top_color, left_color, self.latin_name, define_color_module, self.save_path)
                    p.main(ancestral_chromosome_conf, [])
            event.set()
            import time
            time.sleep(1)
        except Exception as e:
            print('\n')
            print('\033[1;31mPlease check the input file and parameters!!!\033[0m')
            event.set()
        else:
            print('\n')
            print('\033[0;32mCompleted.........\033[0m')  


if __name__ == '__main__':

    # conf = 'conf/ancestral_karyotype_iteration.conf'
    conf = 'ancestral_karyotype_iteration.conf'
    
    options = load_conf(conf, 'ancestral')
    p = ConservedAncestralRegions(options)
    p.run()


"""
module手动定义好之后  需要对这些module定义相同的颜色 初步可以绘图
之后统计每个染色体上的module 并且根据染色体的顺序对module进行排序 初步定义染色体的连接模式
连接模式也有助于判断 祖先状态

之后在根据最大交集的思路推断


祖先染色体核型重构主要是利用代表性物种之间的比较，其核心思想是：对任意两个基因组及外类群（或者参考基因组),
鉴定两基因组最近祖先节点下每个基因组中最大独立的共线性区域（ISR: independent syntenic region）。
假如ISR相对外类群（参考基因组）基因组具有一致的block链接式样，则ISR为极大似然的祖先染色体；
假如ISR相对外类群（参考基因组）具有不一致的链接式样，那么与外类群链接式样一致的ISR为祖先染色体。
基于这一思想，不仅可以获得物种分歧之前的染色体核型结构；并且可以获得一个物种多倍化事件前后的染色体。


对任意两个基因组及外类群（或者参考基因组），鉴定两基因组最近祖先节点下每个基因组中最大连续共线性区域
（CSRs: contiguous syntenic regions），然后比较CSR在被考虑基因组中的保守祖先区域（CARs：conservative ancestral regions）。
在两个被考虑的基因组中，那些具有一致基因组共线性区域的最大CSRs，被推定为CARs；
那些具有不一致共线性区域的最大CSRs中，与外类群（参考基因组）一致的CSRs，被推定为CARs。
基于这一思路，最终推定两基因组最近祖先的CARs，即为极大似然的祖先染色体结构；
并且对于一个物种特定多倍化事件前祖先染色体核型也可以基于此原则推定。


如果有多个物种，那么我们是根据那个树的关系总会确定出一个符合: 两个基因组及外类群  Vvi_Aya_Tsi 
                +--Aya
            +--|
            |   +--Tsi
        +---|
        |   +------Vvi
    ----|   
        +----------Bsi

这种树的关系可以根据输入的物种定义 Vvi_Aya_Tsi  树形由内像外定义一个顺序
如果保守区域 映射到参考/比较物种染色体上，如果都在在一条染色体上那么是合并的、否则是分开的
在参考中是分开的， 在比较物种中只要有一个相同的结构  那就就是那种结构
"""

