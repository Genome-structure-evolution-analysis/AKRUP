import copy
import random
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
        self.infer_name = 'refspec_spec1,refspec_spec2'
        self.intergenomicratio = 1 # orthologous synteny ratio (default: 1)
        self.block_num = 5
        self.hocv = -1 # (-1< hocv <1)

        self.common_wgd = False
        self.infer_wgd_flag = True
        self.Conserved_spec = sepc
        
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
        self.infer_name = self.infer_name.split(',')
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

    def product_ancestral_chromosome(self, species, define_color, CSR_lens, alpla):
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
            ancestral_color_lens[ancestral_order] = [1, ancestral_color_lens[ancestral_order], color, alpla]
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
                    conf_line = [c, order, order+row[0], row[1], alpla]
                    midd_chromosome_conf.append(conf_line)
                    order += row[0]


            wgds_chrs_color.append(midd_chromosome_conf)


        return wgds_names, wgds_chrs_color
        # print(spec_each_lens)

    def map_spec_modeule(self, spec_names, mdinfos):

        modules_colorinfo, modules_infos = {}, {}
        specs = spec_names.split('_')
        ref_spec = specs[0]
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

    def estimation_confidence(self, species, intergenomicratio, chr_types, color_pos):  # estimation Confidence
        specs = species.split('_')
        repeat_results, CSR_link = [], []
        for i in range(1, 1001):
            print('Inferring ancestral karyotype %sst time'%(i))
            color_module, ancestral_chr = self.get_adjacent_between_module(species, intergenomicratio, copy.deepcopy(chr_types), color_pos)
            ancestral_chr = [[str(x[0]), x[1][0], x[2]] for x in ancestral_chr]
            repeat_results.append([color_module, ancestral_chr])
            CSR_link.append(ancestral_chr)


        # 统计每种连接方式出现的次数
        t = lambda x: tuple(map(str.strip, x))
        c = Counter(map(t, chain.from_iterable(CSR_link)))

        ancs = []
        flag, ancnum = repeat_results[0], 1

        for anc in repeat_results[1:]:
            if flag[0] == anc[0] and flag[1] == anc[1]:
                ancnum += 1
            else:
                ancs.append([ancnum, flag[0], flag[1]])
                flag = anc
                ancnum = 1
        ancs.append([ancnum, flag[0], flag[1]])

        new_ancs = sorted(ancs, key=lambda x: int(x[0]), reverse=True)

        print()
        print()
        print()
        print('Get the most likely ancestral karyotype......')
        best_anc = new_ancs[0]

        define_color_module = best_anc[1]
        total_cof = best_anc[0]/1000*100
        print(f'Ancestral karyotype confidence \033[1;31m{total_cof}%\033[0m')

        print()
        print(f'chromosome order\t{specs[1]}_chr\t{specs[2]}_chr\tConnection pattern\tConfidence')
        chrnum = 1
        for lis in best_anc[2]:
            times = c[tuple(lis)]
            cof = float(times)/1000*100
            print(f'chromosome:{chrnum}\t{lis[0]}\t{lis[1]}\t{lis[2]}\t{cof}%')
            chrnum += 1

        return define_color_module

    @staticmethod
    def p_random(arr1, arr2):
        assert len(arr1) == len(arr2), "Length does not match."
        assert sum(arr2) == 1 , "Total rate is not 1."

        sup_list = [len(str(i).split(".")[-1]) for i in arr2]
        top = 10 ** max(sup_list)
        new_rate = [int(i*top) for i in arr2]
        
        rate_arr = []
        for i in range(1,len(new_rate)+1):
            rate_arr.append(sum(new_rate[:i]))
        rand = random.randint(1,top)
        
        data = None
        for i in range(len(rate_arr)):
            if rand <= rate_arr[i]:
                data = arr1[i]
                break
        return data

    def get_adjacent_between_module(self, species, sub_num, chr_all_types, color_pos):
        
        define_type_module, colors = {}, list(self.colors.values())

        specs = species.split('_')
        ref_spec = specs[0]
        init_key = ref_spec+'_'+specs[1]
        init_types = chr_all_types[init_key] 

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
                    select_CSR = self.p_random(a, b)

                    # max_ty = Counter(CSRs).most_common(1)[0]
                else:
                    continue
                if select_CSR:
                

                    ancestral_chr.append([chr_key, link_dit[select_CSR], select_CSR])

                    tyss = select_CSR.split('_')
                    color, flag = colors[color_pos], 1
                    if any(([(t in define_type_module) for t in tyss])):
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

        return define_type_module, ancestral_chr

    def sort_method(self, args):
        a = sorted([args[2], args[3]])
        b = sorted([args[4], args[5]])
        return a+b

    def compair_type_distance(self, a, b):
        if a == 0 or b == 0:
            min_dis = 99999
        else:
            distance = []
            Balanced = [0.5, 0.5]
            
            x = np.array(Balanced)
            y = np.array([a, b])
            d = np.sqrt( np.sum ( np.square ( x - y ) ) )
            distance.append(d)
            min_dis = min(distance)

        return min_dis

    def get_adjacent_within_module(self, sepcies, target_spec, wgd_num, chr_all_types, define_recent_wgd_file, flagpr=True):
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
        print(f'chromosome order\tref_chr\tcom_chr\tConnection pattern')
        chrnum = 1
        for lis in wgd_recent_order:
            link_model = '_'.join(lis[4])
            print(f'chromosome:{chrnum}\t{lis[0]}\t{lis[1]}\t{link_model}')
            chrnum += 1
        if flagpr:
            print(f'WGD {na} Before duplication: {len(wgd_recent_order)} chromosome')
        deal_chr = sorted(list(set([x[0] for x in wgd_recent_order])), key=lambda x: int(x))

        for dc in deal_chr:
            if dc in com_types:
                com_types.pop(dc)
                
        while wgds:

            wgd = wgds.pop(-1)
            wgd_chr, wgd_color_order = [], []
            num = 0
            print('-------------------------------------------------')
            print('start')
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

                        midd_wgd_color_order.append([str(chr_key), lis[1], lis[2], lis[3], new_order])

                flag, new_midd_wgd_color_order = [], []
                for lis in midd_wgd_color_order:
                    lkk = f'{lis[0]}_{lis[1]}'
                    lkk1 = f'{lis[1]}_{lis[0]}'
                    if lkk not in flag and lkk1 not in flag:
                        new_midd_wgd_color_order.append(lis)
                        flag.append(lkk)
                        flag.append(lkk1)

                sort_midd_wgd_color_order = sorted(new_midd_wgd_color_order, key=lambda x: (len(x[4]), float(x[2]), float(x[3])), reverse=True)


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
                        for k in lk[:2]:
                            remove_tt, new_tt = [], []
                            
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

            chrnum1 = 1
            for lis in result_wgd_color_order:
                link_model = '_'.join(lis[4])
                print(f'chromosome:{chrnum1}\t{lis[0]}\t{lis[1]}\t{link_model}')
                chrnum1 += 1
            if flagpr:
                print(f'WGD {wgd} Before duplication: {len(wgd_color_order)} chromosome')

            print('-------------------------------------------------')
            deal_chr = [x[0] for x in result_wgd_color_order]
            for dc in deal_chr:
            # for dc in ['2', '10', '19', '5', '26', '3']:
                if dc in com_types:
                    com_types.pop(dc)


        return all_wgds_chromosome

    def define_ancestral_color(self, all_wgds_chromosome, color_pos, chr_spec_types):
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

    def infer_ancestral_karyotype(self, modules_infos, lens_dit):

        no_exist_CSR = []
        all_key = list(modules_infos.keys())
        if len(all_key) > 1:  
            ll_CSR = [x[14] for x in modules_infos[all_key[0]]]
            yy_CSR = [x[14] for x in modules_infos[all_key[1]]]
            no_exist_CSR.extend(list(set(ll_CSR)-set(yy_CSR)))
            no_exist_CSR.extend(list(set(yy_CSR)-set(ll_CSR)))

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

        temp_list = []
        block_list = sorted(block_list, key=lambda x: x[25])
        temp_list = groupby(block_list, lambda x: x[25])

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
        new_color_list = []

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

    def main(self):
        ref_spec = self.species.split('_')[0]
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
        chr_all_types = self.infer_ancestral_karyotype(new_mdinfos, spec_lens_dit)


        # The total length of identical CSRs in all chromosomes of each species was counted.
        # For inferring the proportion of ancestral chromosomes

        # The proportion of identical CSRs in each chromosome in each species was counted and
        # used to infer the proportion of chromosomes before and after WGD
        CSR_top_all_lens, CSR_top_each_len = self.count_CSR_lens(new_mdinfos)


        # commond polyploidization to infer ancestral karyotypes
        if self.common_wgd.upper() == 'TRUE':
            kkey = f'{ref_spec}_{self.Conserved_spec}'
            all_wgds_chromosome = self.get_adjacent_within_module(self.species, self.Conserved_spec, self.wgds, copy.deepcopy(chr_all_types), self.recentwgdchr[self.Conserved_spec], False)  # '2_3' 一次2倍一次三倍
            define_color_module = self.define_ancestral_color(all_wgds_chromosome, self.color_pos, chr_all_types[kkey])
        # Inferring ancestral karyotype
        else:
            define_color_module = self.estimation_confidence(self.species, self.intergenomicratio, copy.deepcopy(chr_all_types), self.color_pos)
            # define_color_module, ancestral_chr = self.get_adjacent_between_module(self.species, self.intergenomicratio, copy.deepcopy(chr_all_types), self.color_pos)

        ancestral_chromosome_conf = self.product_ancestral_chromosome(self.species, define_color_module, CSR_top_all_lens, 1)


        # save ancestral file
        sf_anc_chromosome = open(f'A-{self.species}-ancestral_chromosome_conf.txt', 'w')
        ancestral_chromosome_conf = sorted(ancestral_chromosome_conf, key=lambda x: int(x[0]))
        new_anc_color = ['\t'.join([str(v) for v in x]) for x in ancestral_chromosome_conf]
        sf_anc_chromosome.write('\n'.join(new_anc_color))

        # Save the ancestral color of the CSR
        sf_anc_CSR_color = open(f'A-{self.species}-ancestral_CSR-color_conf.txt', 'w')
        for CSR, cor in define_color_module.items():
            sf_anc_CSR_color.write(f'{CSR}\t{cor}\n')

        
        # import threading
        event = threading.Event()
        print()
        content = 'Drawings are in progress'
        th = threading.Thread(target=await_run, args=(content,event))
        th.start()

        for name in self.infer_name:
            spec1, spec2 = name.split('_')

            # Karyotype color block for expansion
            top_color, left_color = self.get_module_block_color(block_list_dit[name], define_color_module, spec_lens_dit, name)
            
            # save Karyotype color
            sf_top = open(f'A.{name}.top.color.pos.txt', 'w')
            sf_left = open(f'A.{name}.left.color.pos.txt', 'w')
            new_top_color = ['\t'.join([str(v) for v in x]) for x in top_color]
            new_left_color = ['\t'.join([str(v) for v in x]) for x in left_color]
            sf_top.write('\n'.join(new_top_color))
            sf_left.write('\n'.join(new_left_color))

            if self.infer_wgd_flag.upper() == 'TRUE':
                # Inferring karyotypes before and after polyploidization
                target_spec = spec2
                all_wgds_chromosome = self.get_adjacent_within_module(self.species, target_spec, self.wgds, copy.deepcopy(chr_all_types), self.recentwgdchr[target_spec])  # '2_3' 一次2倍一次三倍
                wgds_names, wgds_chr_colors = self.product_wgd_chromosome(self.species, target_spec, define_color_module, CSR_top_each_len, all_wgds_chromosome, self.wgds, 1)

                # Karyotype mapping
                p = DotplotBlock(name, lens_list, block_list_dit[name], top_color, left_color, self.latin_name, define_color_module)
                p.main(ancestral_chromosome_conf, wgds_chr_colors)
            elif self.infer_wgd_flag.upper() == 'FALSE':
                p = DotplotBlock(name, lens_list, block_list_dit[name], top_color, left_color, self.latin_name, define_color_module)
                p.main(ancestral_chromosome_conf, [])
        event.set()

