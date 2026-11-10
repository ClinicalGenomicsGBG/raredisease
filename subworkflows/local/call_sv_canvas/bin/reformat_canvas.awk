BEGIN {
  FS = "\t"
  OFS = "\t"
  if (expected_x_cn == 1) {
    gonosome = "chr[XY]"
    autosome = "chr[0-9]+"
  } else {
    # Treat X as an autosome if we expect CN 2 on it
    gonosome = "chrY"
    autosome = "chr[0-9]+|chrX"
  }

}

# Insert additional ALT lines into the header of a Canvas VCF file.
# Should appear below DUP allele because we can reuse that
/^##ALT=<ID=DUP/ {
  print
  print "##ALT=<ID=DEL, Description=\"Deletion\">"
  print "##ALT=<ID=DEL:HEMI, Description=\"Hemizygous deletion\">"
  print "##ALT=<ID=DEL:HOM, Description=\"Homozygous deletion\">"
  print "##ALT=<ID=LOH, Description=\"Loss of heterozygosity\">"
  next
}

# Skip all other symbolic ALT alleles
/^##ALT/ {
  next
}

# Print all header lines verbatim
/^#/ {
  print
  next
}

# Now we've only got records left
# Set up variables for the relevant fields for each line
{
  chr = $1
  format =  $9
  sample_data = $10
  split(format, format_fields, ":")
  for (i in format_fields) {
    if (format_fields[i] == "CN") {
      cn_index = i
    } else if (format_fields[i] == "MCC") {
      mcc_index = i
    }
  }
  split(sample_data, sample_fields, ":")

  cn = sample_fields[cn_index]
  mcc = sample_fields[mcc_index]
  # ALT will be <REF> by default, but we will change it if we meet certain conditions
  new_alt = "<REF>"
}

# Assign new ALT values based on CN and MCC values, and whether the chromosome is a gonosome or autosome

cn == 0 && chr ~ autosome {
  new_alt = "<DEL:HOM>"
}

cn == 0 && chr ~ gonosome {
  new_alt = "<DEL:HEMI>"
}

cn == 1 && chr ~ autosome {
  new_alt = "<DEL>"
}

cn == 2 && mcc == 2 {
  new_alt = "<LOH>"
}

cn > 2 || (cn == 2 && chr ~ gonosome ) {
  new_alt = "<DUP>"
}

# Print the modified line with the new ALT fields
{
  $5 = new_alt
  print
}
