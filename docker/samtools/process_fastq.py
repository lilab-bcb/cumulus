import sys

def convert_to_10x_fastq(in_file, out_name):
    with open(in_file, 'r') as fin, open(f"{out_name}_R1.fastq", 'w') as fout1, open(f"{out_name}_R2.fastq", 'w') as fout2:
        while True:
            line1 = fin.readline()
            if not line1:
                break

            # Read 1 - Line 1
            l1_list = line1.strip().split('\t')

            r1_read_list = ':'.join(l1_list[1].split(':')[2:]).split('@')
            r1_l1_str = '@' + r1_read_list[0] + ' 1' + r1_read_list[1][1:]
            fout1.write(r1_l1_str + '\n')

            # Read 1 - Line 2
            cr = l1_list[-2]
            assert cr[0:2] == 'CR'
            qx = l1_list[2]
            assert qx[0:2] == 'QX'
            fout1.write(''.join([cr.split(':')[2], qx.split(':')[2]]) + '\n')

            # Read 1 - Line 3
            fout1.write('+\n')

            # Read 1 - Line 4
            cy = l1_list[-1]
            assert cy[0:2] == 'CY'
            oq = l1_list[-3]
            assert oq[0:2] == 'OQ'
            fout1.write(''.join([cy.split(':')[2], oq.split(':')[2]]) + '\n')

            # Read 2 - Line 1
            r2_read_list = ':'.join(l1_list[1].split(':')[2:]).split('@')
            r2_l1_str = '@' + r2_read_list[0] + ' ' + r2_read_list[1]
            fout2.write(r2_l1_str + '\n')

            line2 = fin.readline()
            fout2.write(line2)

            line3 = fin.readline()
            fout2.write(line3)

            line4 = fin.readline()
            fout2.write(line4)
            


if __name__ == '__main__':
    in_file = sys.argv[1]
    out_name = sys.argv[2]

    convert_to_10x_fastq(in_file, out_name)
