from fpdf import FPDF
import sys

def convert_txt_pdf(report_file,output_file):
    pdf=FPDF()
    pdf.add_page()
    with open(report_file,"r") as rh:
        for x in rh:
            if "Section" in x or "Total number of reads" in x:
                pdf.set_font('Arial', 'B', 15)
                pdf.cell(2, 10, txt = x, ln = 1, align = 'L')
            else:
                pdf.set_font('Arial', '', 10)
                pdf.cell(2, 5, txt = x, ln = 1, align = 'L')
        pdf.output(output_file)

if __name__ == "__main__":
    report_file=sys.argv[1]
    output_file=sys.argv[2]
    convert_txt_pdf(report_file,output_file)