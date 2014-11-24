#!/bin/awk -f
BEGIN {
#	sumr = 0;
#	sume = 0;
#	cont = 0;
#	sumr2 = 0;
#	sume2 = 0;
#	summul = 0;
	found = 0; 
	field = "REPORTED QUALITY VS EMPIRICAL:"; 
}
{
	if(found == 1){
		# End of data
		if($0 ~ /=/){
			exit 0;
		}

		# Data found
		if($2 != 0){
			print $0;
#			sumr += $1;
#			sume += $2;
#			cont++;
#			sumr2 += $1 * $1;
#			sume2 += $2 * $2;
#			summul += $1 * $2;
		}
	}
	else{
		# Searching for data
		if($0 == field){
			found = 1;
		}
	}
}
END {
#	avgr = sumr / cont;
#	avge = sume / cont;
#	cov = (summul / cont) - (avgr * avge); 
#	devr = sqrt( (sumr2 / cont) - (avgr * avgr) );
#	deve = sqrt( (sume2 / cont) - (avge * avge) );
#	r = cov / (devr * deve);
#	print r " 0";
}
