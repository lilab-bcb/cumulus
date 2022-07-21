import pegasus as pg
from matplotlib.backends.backend_pdf import PdfPages
import sys

def chunk_list(var_names, chunk_size):
    for i in range(0, len(var_names), chunk_size):
        yield var_names[i:i + chunk_size]

def load_generate_plot(count_matrix_file,count_matrix,output_plot_pdf):
    data = pg.read_input(count_matrix_file)
    pg.arcsinh(data)
    with PdfPages(output_plot_pdf) as ridgep_figs:
        firstPage = plt.figure(figsize=(10, 10))
        firstPage.clf()
        txt = count_matrix
        firstPage.text(0.5,0.5,txt, transform=firstPage.transFigure, size=24, ha="center")
        ridgep_figs.savefig()
        plt.close()
        for feat_list in chunk_list(data.var_names,8):
            fig = pg.ridgeplot(data,feat_list[0:len(feat_list)],return_fig=True)
            fig.set_figheight(10)
            fig.set_figwidth(10)
            ridgep_figs.savefig(fig, pad_inches=1, bbox_inches="tight")

if __name__ == "__main__":
    count_matrix_file=sys.argv[1]
    count_matrix=sys.argv[2]
    output_plot_pdf=sys.argv[3]
    load_generate_plot(count_matrix_file,count_matrix,output_plot_pdf)