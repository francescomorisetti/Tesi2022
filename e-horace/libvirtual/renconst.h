* renconst.h
* the declarations for renconst.F
* generated by WriteRenConst 7 Mar 2005 17:55

	double complex dMWsq1, dMZsq1, dSW1, dZAA1, dZe1, dZW1, dZZA1
	double complex dZfL1(4,1,1)
	common /renconst/ dMWsq1, dMZsq1, dSW1, dZAA1, dZe1, dZW1
	common /renconst/ dZZA1, dZfL1

	integer sizeof_rc
	parameter (sizeof_rc = 11)
	double complex rc(sizeof_rc)
	equivalence (dMWsq1, rc)

