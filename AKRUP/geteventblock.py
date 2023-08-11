from scipy.optimize import curve_fit
from scipy.stats import gaussian_kde, linregress

from AKRUP.funcbase import *
from AKRUP.dotplot_CSR_event import DotplotCsrKs


class KidFit:
    def __init__(self, name, all_ks_list, corret=1, peaks=1, area=[0, 3], bins_number=200):
        self.name = name
        self.bins_number = bins_number
        self.area = area
        self.peaks = peaks
        self.corret = corret
        self.all_ks_list = all_ks_list
        sf_distribute_file = f'{name}.ks_distribute.txt'
        self.sf_ks_handle = open(sf_distribute_file, 'w')
        self.sf_ks_handle.write(',a1,mu1,sigma1\n')

    def gaussian_fuc(self, x, *params):
        y = np.zeros_like(x)
        for i in range(0, len(params), 3):
            amp = float(params[i])
            ctr = float(params[i+1])
            wid = float(params[i+2])
            y = y + amp * np.exp(-((x - ctr)/wid)**2)
        return y

    def kde_fit(self, data, x, name):
        kde = gaussian_kde(data)
        kde.set_bandwidth(bw_method=kde.factor/3.)
        p = kde(x)
        guess = [1, 1, 1]*self.peaks
        popt, pcov = curve_fit(self.gaussian_fuc, x, p, guess, maxfev=80000)
        popt = [abs(k) for k in popt]
        data = []
        y = self.gaussian_fuc(x, *popt)
        for i in range(0, len(popt), 3):
            array = [popt[i], popt[i+1], popt[i+2]]
            data.append(self.gaussian_fuc(x, *array))
        slope, intercept, r_value, p_value, std_err = linregress(p, y)
        print("\nR-square: "+str(r_value**2))
        print("The gaussian fitting curve parameters are :")
        for i in range(0, len(popt), 3):
            popt[i+2] = popt[i+2]/sqrt(2)
            popt[i] = popt[i]*sqrt(2*pi)*popt[i+2]
            popt[i+2] = round(popt[i+2], 4)
            popt[i] = round(popt[i], 4)
            popt[i+1] = round(popt[i+1], 4)
            print('\t'.join([str(popt[i]), str(popt[i+1]), str(popt[i+2])]))
            self.sf_ks_handle.write(','.join([name, str(popt[i]), str(popt[i+1]), str(popt[i+2])]))

        return popt

    def main(self):
        data = [ks*float(self.corret) for ks in self.all_ks_list]
        data = [k for k in data if self.area[0] <= k <= self.area[1]]
        x = np.linspace(self.area[0], self.area[1], self.bins_number)
        peaks = self.kde_fit(data, x, self.name)
        peaks.insert(0, self.name)

        return peaks

class EventBlock:
    def __init__(self, options):

        self.name = 'leftname_topname'
        self.block_info = 'blockinfo_file'
        self.save_file = 'save_file (*.EventRelate_block.information.csv)'
        self.lens_file1 = 'lens1 file'
        self.lens_file2 = 'lens2 file'
        self.range_k = 0.15
        self.pkcolor = 'orange'
        self.pk_hocv = 0.8
        self.pk_block_num = 30
        self.block_num = 5
        self.hocv_depth = 1
        self.hocv = -1
        self.dpi = 300

        for k, v in options:
            if str(k) in ['range_k', 'pk_hocv', 'pk_block_num', 'block_num', 'hocv_depth', 'hocv', 'dpi']:
                setattr(self, str(k), float(v))
            else:
                setattr(self, str(k), v)
        self.hocv_depth = int(self.hocv_depth)

    def select_blockinfo(self, block_list, bk_num, hocv, hocv_depth, peaks, range_k):
        new_block, hocv_pos = [], 10
        hocv_pos += hocv_depth
        for bk in block_list[1:]:
            if not int(bk[8]) >= bk_num:
                continue
            if peaks-range_k <= float(bk[9]) <= peaks+range_k and float(bk[hocv_pos]) >= hocv:
                new_block.append(bk)

        return new_block

    def classific_event_bestpart(self, blockinfos, block_num, hocv, hocv_depth=1):
        ks_s, ks_e = 0, 3
        new_blockinfos, hocv_pos = [], 10
        hocv_pos += hocv_depth
        num = 1
        for b_lines in blockinfos[1:]:
            if int(b_lines[8]) >= block_num and float(b_lines[hocv_pos]) >= hocv:
                b_lines[0] = str(num)
                new_blockinfos.append(b_lines)
                num += 1

        return new_blockinfos

    @staticmethod
    def write_blockinfo(save_handle, block_infos):
        save_handle.write(','.join(['num', 'chr1', 'chr2', 'start1', 'end1', 'start2', 'end2', 'pvalue', 'length',
                             'ks_median', 'ks_average', 'hocv1', 'hocv2', 'hocv3',
                             'hocv4', 'hocv5', 'hocv6', 'hocv7', 'hocv8', 'block1', 'block2', 'ks', 'density1', 'density2', 'CSR', 'module', 'Ts\n']))

        for rows in block_infos:
            save_handle.write(','.join(rows) + '\n')

    def run(self):
        block_list = [x.strip().split(',') for x in open(self.block_info)]
        peaks_blocks = self.classific_event_bestpart(block_list, self.pk_block_num, self.pk_hocv, self.hocv_depth)
        all_ks_list = [float(x[9]) for x in peaks_blocks if float(x[9]) > 0]
        ks_p = KidFit(self.name, all_ks_list)
        peaks = ks_p.main()

        new_blockinfos = self.select_blockinfo(block_list, self.block_num, self.hocv, self.hocv_depth, peaks[2], self.range_k)
        self.save_file_handle = open(self.save_file, 'w')
        self.write_blockinfo(self.save_file_handle, new_blockinfos)

        left_name, top_name = self.name.split('_')
        savefile = f'{self.name}.Ks-event.middle.dotplot.png'
        conf = dict(block_info=self.save_file, hocv_depth=self.hocv_depth, hocv=self.hocv, peaks=','.join([str(x) for x in peaks[1:]]), 
            range_k=self.range_k, genome1_name=left_name, genome2_name=top_name, pkcolor=self.pkcolor, lens_file1=self.lens_file1, 
            lens_file2=self.lens_file2, block_num=self.block_num, savefile=savefile, dpi=self.dpi)
        p = DotplotCsrKs(conf.items(), 'eventdotplot')
        p.run()

