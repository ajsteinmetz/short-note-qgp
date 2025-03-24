import PyPDF2

def merge_pdfs(pdf_files, output_path):
    # Create a PdfMerger object
    merger = PyPDF2.PdfMerger()
    
    # Append each PDF file in the list to the merger
    for pdf in pdf_files:
        merger.append(pdf)
        print(f"Added {pdf}")
    
    # Write the merged PDF to the output file
    with open(output_path, 'wb') as output_file:
        merger.write(output_file)
    print(f"Merged PDF saved as {output_path}")
    
    # Close the merger
    merger.close()

if __name__ == "__main__":
    # List of PDF files to merge; change these to your file names
    pdf_files = ["revisions_letter\\revisions_letter.pdf", "revisions_letter\\revised_short_note_qgp.pdf"]
    
    # Output file name for the merged PDF
    output_path = "revisions_letter_package.pdf"
    
    merge_pdfs(pdf_files, output_path)
