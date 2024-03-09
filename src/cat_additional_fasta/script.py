#!/usr/bin/env python3
"""
Read a custom fasta file and create a custom GTF containing each entry
"""
from itertools import groupby
import logging
import os
import sys

# Create a logger
logging.basicConfig(format="%(name)s - %(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__file__)
logger.setLevel(logging.INFO)


def fasta_iter(fasta_name):
    """
    modified from Brent Pedersen
    Correct Way To Parse A Fasta File In Python
    given a fasta file. yield tuples of header, sequence

    Fasta iterator from https://www.biostars.org/p/710/#120760
    """
    with open(fasta_name) as fh:
        # ditch the boolean (x[0]) and just keep the header or sequence since
        # we know they alternate.
        faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
        for header in faiter:
            # drop the ">"
            headerStr = header.__next__()[1:].strip()

            # join all sequence lines to one.
            seq = "".join(s.strip() for s in faiter.__next__())

            yield (headerStr, seq)


def fasta2gtf(fasta, output, biotype):
    fiter = fasta_iter(fasta)
    # GTF output lines
    lines = []
    attributes = 'exon_id "{name}.1"; exon_number "1";{biotype} gene_id "{name}_gene"; gene_name "{name}_gene"; gene_source "custom"; transcript_id "{name}_gene"; transcript_name "{name}_gene";\n'
    line_template = "{name}\ttransgene\texon\t1\t{length}\t.\t+\t.\t" + attributes

    for ff in fiter:
        name, seq = ff
        # Use first ID as separated by spaces as the "sequence name"
        # (equivalent to "chromosome" in other cases)
        seqname = name.split()[0]
        # Remove all spaces
        name = seqname.replace(" ", "_")
        length = len(seq)
        biotype_attr = ""
        if biotype:
            biotype_attr = f' {biotype} "transgene";'
        line = line_template.format(name=name, length=length, biotype=biotype_attr)
        lines.append(line)

    with open(output, "w") as f:
        f.write("".join(lines))


if __name__ == "__main__":
    add_name = os.path.basename(par['additional_fasta'])
    output = os.path.splitext(add_name)[0] + ".gtf"
    fasta2gtf(par['additional_fasta'], output, par['biotype'])

    with open(par['fasta'], 'r') as f1:
        content1 = f1.read()
    with open(par['additional_fasta'], 'r') as f2:
        content2 = f2.read()
    with open(par['fasta_output'], 'w') as f_out:
        f_out.write(content1 + content2)
    with open(par['gtf'], 'r') as g1:
        g_content1 = g1.read()
    with open(output, 'r') as g2:
        g_content2 = g2.read()
    with open(par['gtf_output'], 'w') as g_out:
        g_out.write(g_content1 + g_content2)

    text = f"{meta['functionality_name']}:\n  python: {sys.version.split()[0]}"

    if par['versions'] and os.path.exists(par['versions']):
        with open(par['versions'], 'a') as f:
            f.write(text + '\n')
        os.rename(par['versions'], par['updated_versions'])
    else:
        with open(par['updated_versions'], 'w') as f:
            f.write(text + '\n')
