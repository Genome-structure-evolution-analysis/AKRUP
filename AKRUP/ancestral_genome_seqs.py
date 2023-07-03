from Bio import SeqIO
from AKRUP.funcbase import *


class ancestral_karyotype:
    def __init__(self, options):
        self.mark = 'seq name'
        self.gff = 'gff file'
        self.ancestor_conf = 'ancestor conf'
        self.pep_file = 'pep file'
        self.ancestor_pep =  'ancestor pep file'
        self.ancestor_gff =  'ancestor gff file'
        self.ancestor_lenstxt = 'ancestor lenstxt'
        self.ancestor_lens =  'ancestor lens'
        self.ancestor_file = 'ancestor file'

        for k, v in options:
            setattr(self, str(k), v)

    def newgff(self, file):
        gff = pd.read_csv(file, sep="\t", header=None, index_col=5)
        gff.rename(columns={0: 'chr', 1: 'start',
                            2: 'end', 3: 'stand', 6: 'order'}, inplace=True)
        gff['chr'] = gff['chr'].astype(str)
        gff['start'] = gff['start'].astype(np.int64)
        gff['end'] = gff['end'].astype(np.int64)
        gff['stand'] = gff['stand'].astype(str)
        gff['order'] = gff['order'].astype(int)
        return gff

    def read_calassfication(self, file):
        classification = pd.read_csv(file, sep="\t", header=None)
        classification[0] = classification[0].astype(str)
        classification[1] = classification[1].astype(int)
        classification[2] = classification[2].astype(int)
        classification[3] = classification[3].astype(str)
        return classification

    def run(self):
        gff = self.newgff(self.gff)
        ancestor = self.read_calassfication(self.ancestor_conf)
        gff = gff[gff['chr'].isin(ancestor[0].values.tolist())]
        newgff = gff.copy()
        data,num = [],1
        chr_arr = ancestor[3].drop_duplicates().to_list()
        chr_dict = dict(zip(chr_arr, range(1, len(chr_arr)+1)))
        ancestor['order'] = ancestor[3].map(chr_dict)
        dict1, dict2 = {}, {}
        for order, group in ancestor.groupby(['order'], sort=[False, False]):
            for index, row in group.iterrows():
                index1 = gff[(gff['chr'] == row[0]) & (
                    gff['order'] >= row[1]) & (gff['order'] <= row[2])].index
                newgff.loc[index1, 'chr'] = str(num)
                for k in index1:
                    data.append(newgff.loc[k, :].values.tolist()+[k])
            # dict1[str(num)] = cla
            dict2[str(num)] = group[3].values[0]
            num+=1
        df = pd.DataFrame(data)
        pep = SeqIO.to_dict(SeqIO.parse(self.pep_file, "fasta"))
        df = df[df[6].isin(pep.keys())]
        for name, group in df.groupby([0]):
            df.loc[group.index, 'order'] = list(range(1, len(group)+1))
            df.loc[group.index, 'newname'] = list(
                [str(self.mark)+str(name)+'g'+str(i).zfill(5) for i in range(1, len(group)+1)])
        df['order'] = df['order'].astype('int')
        df = df[[0, 1, 2, 3, 6, 'newname', 'order']]
        df = df.sort_values(by=[0, 'order'])
        df.to_csv(self.ancestor_gff, sep="\t", index=False, header=None)
        lens = df.groupby(0).max()[[2, 'order']]
        lens1 = df.groupby(0).max()[['order']]
        lens1.to_csv(self.ancestor_lens, sep="\t", header=None)
        lens.to_csv(self.ancestor_lenstxt, sep="\t", header=None)
        lens[1] = 1
        lens['color'] = lens.index.map(dict2)
        # lens['class'] = lens.index.map(dict1)
        lens[[1, 'order', 'color']].to_csv(
            self.ancestor_file, sep="\t", header=None)
        id_dict = df.set_index(6).to_dict()['newname']
        seqs = []
        for seq_record in SeqIO.parse(self.pep_file, "fasta"):
            if seq_record.id in id_dict:
                seq_record.id = id_dict[seq_record.id]
            else:
                continue
            seqs.append(seq_record)
        SeqIO.write(seqs, self.ancestor_pep, "fasta")
