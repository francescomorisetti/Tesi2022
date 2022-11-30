#include "model.h"
#include "looptools.h"
#include "renconst.h"

	double precision S, T, U
	common /kinvars/ S, T, U

	integer Hel1, Hel2, Hel3, Hel4
	common /kinvars/ Hel1, Hel2, Hel3, Hel4

	double complex Pair3, Pair4, Pair9, Pair5, Pair8, Pair2, Pair7
	double complex Pair1, Pair6, Eps4, Eps14, Eps8, Eps10, Eps15
	double complex Eps13, Eps9, Eps6, Eps7, Eps5, Eps12, Eps11
	double complex Eps2, Eps3, Eps1, Abb21, Abb33, Abb36, Abb12
	double complex Abb22, Abb41, Abb49, Abb42, Abb23, Abb24, Abb1
	double complex Abb13, Abb47, Abb37, Abb38, Abb8, Abb2, Abb25
	double complex Abb14, Abb26, Abb50, Abb27, Abb53, Abb9, Abb3
	double complex Abb43, Abb48, Abb4, Abb39, Abb34, Abb40, Abb5
	double complex Abb10, Abb44, Abb11, Abb45, Abb20, Abb16, Abb35
	double complex Abb51, Abb52, Abb6, Abb46, Abb28, Abb15, Abb29
	double complex Abb17, Abb7, Abb30, Abb31, Abb32, Abb18, Abb19
	double complex AbbSum83, AbbSum13, AbbSum18, AbbSum58
	double complex AbbSum46, AbbSum76, AbbSum2, AbbSum20, AbbSum29
	double complex AbbSum72, AbbSum19, AbbSum73, AbbSum84
	double complex AbbSum95, AbbSum33, AbbSum56, AbbSum32
	double complex AbbSum36, AbbSum89, AbbSum60, AbbSum39
	double complex AbbSum93, AbbSum64, AbbSum53, AbbSum66
	double complex AbbSum92, AbbSum59, AbbSum75, AbbSum37
	double complex AbbSum38, AbbSum94, AbbSum70, AbbSum87
	double complex AbbSum91, AbbSum47, AbbSum28, AbbSum25, AbbSum7
	double complex AbbSum80, AbbSum10, AbbSum31, AbbSum69
	double complex AbbSum68, AbbSum44, AbbSum50, AbbSum26
	double complex AbbSum35, AbbSum34, AbbSum51, AbbSum54
	double complex AbbSum52, AbbSum74, AbbSum57, AbbSum71
	double complex AbbSum55, AbbSum41, AbbSum62, AbbSum48, AbbSum5
	double complex AbbSum23, AbbSum27, AbbSum24, AbbSum43, AbbSum9
	double complex AbbSum8, AbbSum82, AbbSum49, AbbSum11, AbbSum12
	double complex AbbSum79, AbbSum77, AbbSum81, AbbSum30
	double complex AbbSum88, AbbSum78, AbbSum14, AbbSum67
	double complex AbbSum65, AbbSum85, AbbSum90, AbbSum15
	double complex AbbSum45, AbbSum17, AbbSum16, AbbSum86
	double complex AbbSum42, AbbSum63, AbbSum1, AbbSum6, AbbSum21
	double complex AbbSum3, AbbSum4, AbbSum61, AbbSum40, AbbSum22
	common /abbrev/ Pair3, Pair4, Pair9, Pair5, Pair8, Pair2
	common /abbrev/ Pair7, Pair1, Pair6, Eps4, Eps14, Eps8, Eps10
	common /abbrev/ Eps15, Eps13, Eps9, Eps6, Eps7, Eps5, Eps12
	common /abbrev/ Eps11, Eps2, Eps3, Eps1, Abb21, Abb33, Abb36
	common /abbrev/ Abb12, Abb22, Abb41, Abb49, Abb42, Abb23
	common /abbrev/ Abb24, Abb1, Abb13, Abb47, Abb37, Abb38, Abb8
	common /abbrev/ Abb2, Abb25, Abb14, Abb26, Abb50, Abb27, Abb53
	common /abbrev/ Abb9, Abb3, Abb43, Abb48, Abb4, Abb39, Abb34
	common /abbrev/ Abb40, Abb5, Abb10, Abb44, Abb11, Abb45, Abb20
	common /abbrev/ Abb16, Abb35, Abb51, Abb52, Abb6, Abb46, Abb28
	common /abbrev/ Abb15, Abb29, Abb17, Abb7, Abb30, Abb31, Abb32
	common /abbrev/ Abb18, Abb19, AbbSum83, AbbSum13, AbbSum18
	common /abbrev/ AbbSum58, AbbSum46, AbbSum76, AbbSum2
	common /abbrev/ AbbSum20, AbbSum29, AbbSum72, AbbSum19
	common /abbrev/ AbbSum73, AbbSum84, AbbSum95, AbbSum33
	common /abbrev/ AbbSum56, AbbSum32, AbbSum36, AbbSum89
	common /abbrev/ AbbSum60, AbbSum39, AbbSum93, AbbSum64
	common /abbrev/ AbbSum53, AbbSum66, AbbSum92, AbbSum59
	common /abbrev/ AbbSum75, AbbSum37, AbbSum38, AbbSum94
	common /abbrev/ AbbSum70, AbbSum87, AbbSum91, AbbSum47
	common /abbrev/ AbbSum28, AbbSum25, AbbSum7, AbbSum80
	common /abbrev/ AbbSum10, AbbSum31, AbbSum69, AbbSum68
	common /abbrev/ AbbSum44, AbbSum50, AbbSum26, AbbSum35
	common /abbrev/ AbbSum34, AbbSum51, AbbSum54, AbbSum52
	common /abbrev/ AbbSum74, AbbSum57, AbbSum71, AbbSum55
	common /abbrev/ AbbSum41, AbbSum62, AbbSum48, AbbSum5
	common /abbrev/ AbbSum23, AbbSum27, AbbSum24, AbbSum43
	common /abbrev/ AbbSum9, AbbSum8, AbbSum82, AbbSum49, AbbSum11
	common /abbrev/ AbbSum12, AbbSum79, AbbSum77, AbbSum81
	common /abbrev/ AbbSum30, AbbSum88, AbbSum78, AbbSum14
	common /abbrev/ AbbSum67, AbbSum65, AbbSum85, AbbSum90
	common /abbrev/ AbbSum15, AbbSum45, AbbSum17, AbbSum16
	common /abbrev/ AbbSum86, AbbSum42, AbbSum63, AbbSum1, AbbSum6
	common /abbrev/ AbbSum21, AbbSum3, AbbSum4, AbbSum61, AbbSum40
	common /abbrev/ AbbSum22

	double complex cint1, cint2, cint3, cint4, cint5, cint6, cint7
	double complex cint8, cint9, cint10, cint11, cint12, cint13
	double complex cint14, cint15, cint16, cint17, cint18, cint19
	double complex cint20, cint21, cint22, cint23, cint24, cint25
	double complex cint26, cint27, cint28, cint29, cint30, cint31
	double complex cint32, cint33, cint34, cint35, cint36, cint37
	double complex cint38, cint39, cint40, cint41, cint42, cint43
	common /loopint/ cint1, cint2, cint3, cint4, cint5, cint6
	common /loopint/ cint7, cint8, cint9, cint10, cint11, cint12
	common /loopint/ cint13, cint14, cint15, cint16, cint17
	common /loopint/ cint18, cint19, cint20, cint21, cint22
	common /loopint/ cint23, cint24, cint25, cint26, cint27
	common /loopint/ cint28, cint29, cint30, cint31, cint32
	common /loopint/ cint33, cint34, cint35, cint36, cint37
	common /loopint/ cint38, cint39, cint40, cint41, cint42
	common /loopint/ cint43

	integer*8 iint1, iint2, iint3, iint4, iint5, iint6, iint7, iint8
	integer*8 iint9, iint10, iint11, iint12, iint13, iint14, iint15
	integer*8 iint16, iint17, iint18, iint19
	common /loopint/ iint1, iint2, iint3, iint4, iint5, iint6
	common /loopint/ iint7, iint8, iint9, iint10, iint11, iint12
	common /loopint/ iint13, iint14, iint15, iint16, iint17
	common /loopint/ iint18, iint19

	double complex Cloop(18), Ctree(1), MatF(18,1)
	common /coeff/ Cloop, Ctree, MatF
