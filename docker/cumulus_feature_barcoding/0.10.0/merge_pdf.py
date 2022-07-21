from PyPDF2 import PdfFileMerger
import sys

def merge_pdf(files_to_merge,output_file):
    merger = PdfFileMerger()
    for pdf in files_to_merge.split(","):
        merger.append(pdf)
    with open(output_file, "wb") as new_file:
        merger.write(new_file) 

if __name__ == "__main__":
    files_to_merge=sys.argv[1]
    output_file=sys.argv[2]
    merge_pdf(files_to_merge,output_file)