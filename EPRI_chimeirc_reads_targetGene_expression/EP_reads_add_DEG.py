import os,sys,math
def load_data1(file_path):
    with open(file_path) as f:
        gene_data = {}
        datalines = f.readlines()
        for line in datalines:
            line = line.strip()
            items = line.split('\t')
            EP_Fold = math.log2((int(items[1])+1)/((int(items[0])+1)*1.34))
            key = items[6]
            gene_data[key] = EP_Fold
    return gene_data


def load_data2(file_path):
    with open(file_path) as f:
        DEG_data = {}
        datalines = f.readlines()
        for line in datalines[1:]:
            line = line.strip()
            items = line.split('\t')
            gene_Fold = items[2]
            padj = items[-1]
            pvalue = items[-2]
            key = items[0]
            DEG_data[key] = gene_Fold+'@'+pvalue+'@'+padj
    return DEG_data

def EP_reads_add_DEG(file1_path,file2_path):
    gene_EP_info = load_data1(file1_path)
    DEG_info = load_data2(file2_path)
    basename = os.path.basename(file1_path).split('.')[0]
    output = open(basename+'_add_DEG.txt','w')
    output.write('Gene_ID\tlog2EP_chimeric_reads_fold\tlog2DEG_fold\tSignificant\n')
    for key in gene_EP_info:
        if key in DEG_info:
            gene_Fold = DEG_info[key].split('@')[0]
            padj = DEG_info[key].split('@')[2]
            pvalue = DEG_info[key].split('@')[1]
            try:
                if float(gene_Fold) >= 1 and float(pvalue) <0.05:
                    output.write(key+'\t'+str(gene_EP_info[key])+'\t'+str(gene_Fold)+'\t'+'up'+'\n')
                elif float(gene_Fold) <= -1 and float(pvalue) <0.05:
                    output.write(key+'\t'+str(gene_EP_info[key])+'\t'+str(gene_Fold)+'\t'+'down'+'\n')
            #if abs(float(gene_Fold)) <= 1 and float(padj)>=0.05:
                else:
                    output.write(key+'\t'+str(gene_EP_info[key])+'\t'+str(gene_Fold)+'\t'+'none'+'\n')
            except ValueError:
                # 处理无法转换为浮点数的情况
                output.write(key + '\t' + str(gene_EP_info[key]) + '\t' + gene_Fold + '\t' + 'NA'+'\n')

if __name__ == "__main__":

    EP_reads_add_DEG(sys.argv[1],sys.argv[2])

