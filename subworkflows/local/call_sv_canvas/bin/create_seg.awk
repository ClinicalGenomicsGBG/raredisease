# Set up track for IGV
BEGIN {
  print "#track graphType=points maxHeightPixels=300:300:300 color=0,0,220 altColor=220,0,0"
}

# Only output values for locations with coverage
$4 > 0 && $7 > 0 {
  cnv = $4
  ncov = $7
  # log is natural log but we want log2.
  cnvlog = log(cnv) / log(2)
  covlog = log(ncov) / log(2)

  # 0 should indicate no CNV call in our track.
  # On the gonosomes, we expect CN 1, whose log is already 0
  # but on the autosomes, CN 2 is expected, whose log is 1, so we adjust accordingly
  if ($1 !~ /X|Y/) {
    cnvlog = cnvlog - 1
    covlog = covlog - 1
  }

  print "Observed_CNVs\t" $1 "\t" $2 "\t" $3 "\t" covlog
  print "Called CNVs\t" $1 "\t" $2 "\t" $3 "\t" cnvlog
}
