#!/usr/bin/env Rscript

require(docopt)
require(PepBedR)

'Usage:
build_pepbed_report.R [-i <input> -ref <genome reference> -o <output>]

Options:
-i Input peptide bed file
-ref Reference genome ensembl version [default: hg38]
-o Output PDF file [default: report.pdf]

]' -> doc

opts <- docopt(doc)

PepBedR::generate_pepbed_report(bedfile = opts$i, ref_genome = opts$ref, output_file = opts$o)
