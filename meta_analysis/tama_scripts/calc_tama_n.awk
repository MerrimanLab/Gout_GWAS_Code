# Short script to calculate total N from the pre-TAMA summary data:

NR == 1 {
	print $1, $2, $3, "N_AFR", "N_EAS", "N_EUR", "N_LAT", "N_Total";
}

NR > 1 {
	print $1, $2, $3, $7, $12, $17, $22, $7 + $12 + $17 + $22;
}
