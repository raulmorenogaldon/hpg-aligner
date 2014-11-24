#!/bin/sh

$(dirname $0)/extract_recal_graphs.awk $1 > chart.tmp
$(dirname $0)/extract_recal_graphs.awk $2 > chart2.tmp

export chart_filename=$3 
$(dirname $0)/plot_recal_graphs.gp 

#rm chart.tmp
#rm chart2.tmp
