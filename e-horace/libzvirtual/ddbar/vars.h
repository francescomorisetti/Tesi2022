#include "model.h"
#include "looptools.h"
#include "renconst.h"

	double precision S, T, U
	common /ddkinvars/ S, T, U

	integer Hel1, Hel2, Hel3, Hel4
	common /ddkinvars/ Hel1, Hel2, Hel3, Hel4

	double complex Pair5, Pair6, Pair14, Pair13, Pair7, Pair11
	double complex Pair3, Pair4, Pair12, Pair9, Pair2, Pair10
	double complex Pair1, Pair8, Eps25, Eps29, Eps30, Eps31, Eps32
	double complex Eps35, Eps21, Eps33, Eps18, Eps19, Eps20, Eps27
	double complex Eps16, Eps14, Eps28, Eps34, Eps15, Eps12, Eps26
	double complex Eps13, Eps22, Eps23, Eps7, Eps24, Eps8, Eps9
	double complex Eps10, Eps5, Eps6, Eps11, Eps2, Eps3, Eps4
	double complex Eps17, Eps1, Abb273, Abb286, Abb16, Abb216
	double complex Abb230, Abb259, Abb159, Abb203, Abb221, Abb239
	double complex Abb211, Abb192, Abb217, Abb200, Abb160, Abb233
	double complex Abb206, Abb245, Abb268, Abb246, Abb231, Abb290
	double complex Abb288, Abb262, Abb56, Abb223, Abb234, Abb256
	double complex Abb92, Abb267, Abb107, Abb269, Abb235, Abb247
	double complex Abb25, Abb18, Abb342, Abb294, Abb20, Abb57
	double complex Abb58, Abb26, Abb97, Abb114, Abb41, Abb59
	double complex Abb108, Abb115, Abb60, Abb36, Abb295, Abb104
	double complex Abb343, Abb61, Abb54, Abb212, Abb122, Abb224
	double complex Abb218, Abb219, Abb95, Abb335, Abb55, Abb157
	double complex Abb283, Abb207, Abb163, Abb201, Abb254, Abb123
	double complex Abb257, Abb270, Abb164, Abb195, Abb196, Abb250
	double complex Abb236, Abb260, Abb1, Abb334, Abb347, Abb318
	double complex Abb319, Abb336, Abb85, Abb384, Abb385, Abb27
	double complex Abb21, Abb320, Abb105, Abb337, Abb62, Abb274
	double complex Abb193, Abb278, Abb263, Abb63, Abb225, Abb284
	double complex Abb248, Abb93, Abb213, Abb109, Abb167, Abb161
	double complex Abb124, Abb28, Abb2, Abb307, Abb376, Abb375
	double complex Abb316, Abb296, Abb24, Abb321, Abb53, Abb297
	double complex Abb64, Abb86, Abb133, Abb360, Abb304, Abb134
	double complex Abb29, Abb365, Abb162, Abb125, Abb126, Abb197
	double complex Abb127, Abb3, Abb322, Abb204, Abb128, Abb240
	double complex Abb165, Abb271, Abb191, Abb198, Abb251, Abb208
	double complex Abb275, Abb291, Abb241, Abb226, Abb129, Abb255
	double complex Abb121, Abb348, Abb356, Abb358, Abb135, Abb338
	double complex Abb37, Abb323, Abb90, Abb136, Abb388, Abb98
	double complex Abb116, Abb50, Abb137, Abb4, Abb138, Abb209
	double complex Abb303, Abb139, Abb130, Abb7, Abb317, Abb361
	double complex Abb305, Abb140, Abb339, Abb99, Abb366, Abb65
	double complex Abb252, Abb141, Abb10, Abb142, Abb66, Abb205
	double complex Abb67, Abb289, Abb42, Abb272, Abb94, Abb199
	double complex Abb264, Abb43, Abb315, Abb100, Abb214, Abb227
	double complex Abb68, Abb69, Abb87, Abb349, Abb357, Abb344
	double complex Abb143, Abb351, Abb40, Abb324, Abb91, Abb292
	double complex Abb144, Abb345, Abb215, Abb70, Abb51, Abb145
	double complex Abb5, Abb350, Abb372, Abb368, Abb146, Abb44
	double complex Abb110, Abb174, Abb71, Abb8, Abb293, Abb362
	double complex Abb306, Abb147, Abb352, Abb279, Abb367, Abb72
	double complex Abb281, Abb148, Abb11, Abb149, Abb73, Abb166
	double complex Abb194, Abb202, Abb131, Abb168, Abb52, Abb389
	double complex Abb6, Abb386, Abb117, Abb9, Abb308, Abb74
	double complex Abb12, Abb45, Abb46, Abb158, Abb132, Abb232
	double complex Abb261, Abb302, Abb265, Abb210, Abb242, Abb101
	double complex Abb243, Abb118, Abb253, Abb249, Abb170, Abb346
	double complex Abb175, Abb298, Abb38, Abb176, Abb75, Abb96
	double complex Abb390, Abb13, Abb355, Abb377, Abb325, Abb326
	double complex Abb30, Abb76, Abb84, Abb177, Abb299, Abb77
	double complex Abb285, Abb22, Abb380, Abb34, Abb378, Abb382
	double complex Abb150, Abb102, Abb276, Abb151, Abb244, Abb363
	double complex Abb78, Abb35, Abb379, Abb381, Abb383, Abb391
	double complex Abb152, Abb103, Abb277, Abb153, Abb280, Abb364
	double complex Abb79, Abb17, Abb31, Abb266, Abb222, Abb171
	double complex Abb119, Abb172, Abb258, Abb220, Abb287, Abb111
	double complex Abb228, Abb112, Abb237, Abb32, Abb14, Abb329
	double complex Abb314, Abb169, Abb340, Abb369, Abb178, Abb353
	double complex Abb330, Abb310, Abb23, Abb179, Abb106, Abb180
	double complex Abb80, Abb19, Abb341, Abb181, Abb300, Abb313
	double complex Abb331, Abb309, Abb88, Abb182, Abb387, Abb370
	double complex Abb183, Abb47, Abb312, Abb184, Abb120, Abb113
	double complex Abb154, Abb89, Abb327, Abb185, Abb354, Abb186
	double complex Abb332, Abb311, Abb229, Abb373, Abb371, Abb81
	double complex Abb282, Abb374, Abb359, Abb39, Abb301, Abb173
	double complex Abb82, Abb15, Abb328, Abb333, Abb187, Abb188
	double complex Abb238, Abb83, Abb33, Abb48, Abb189, Abb155
	double complex Abb49, Abb392, Abb190, Abb156, AbbSum60
	double complex AbbSum120, AbbSum892, AbbSum891, AbbSum1812
	double complex AbbSum1020, AbbSum579, AbbSum4, AbbSum1683
	double complex AbbSum1206, AbbSum1167, AbbSum640, AbbSum180
	double complex AbbSum1363, AbbSum1457, AbbSum504, AbbSum812
	double complex AbbSum61, AbbSum455, AbbSum1841, AbbSum1218
	double complex AbbSum815, AbbSum769, AbbSum1826, AbbSum538
	double complex AbbSum1963, AbbSum766, AbbSum74, AbbSum171
	double complex AbbSum77, AbbSum1965, AbbSum457, AbbSum770
	double complex AbbSum956, AbbSum440, AbbSum12, AbbSum1770
	double complex AbbSum1100, AbbSum1657, AbbSum95, AbbSum1033
	double complex AbbSum1569, AbbSum147, AbbSum107, AbbSum1129
	double complex AbbSum845, AbbSum138, AbbSum1343, AbbSum1500
	double complex AbbSum122, AbbSum818, AbbSum2038, AbbSum2043
	double complex AbbSum2214, AbbSum1750, AbbSum1232, AbbSum202
	double complex AbbSum762, AbbSum1386, AbbSum251, AbbSum341
	double complex AbbSum1744, AbbSum1632, AbbSum355, AbbSum1945
	double complex AbbSum1864, AbbSum420, AbbSum236, AbbSum789
	double complex AbbSum1376, AbbSum829, AbbSum878, AbbSum1039
	double complex AbbSum595, AbbSum246, AbbSum285, AbbSum1352
	double complex AbbSum800, AbbSum83, AbbSum1048, AbbSum278
	double complex AbbSum399, AbbSum997, AbbSum1108, AbbSum1860
	double complex AbbSum2102, AbbSum303, AbbSum893, AbbSum1300
	double complex AbbSum2188, AbbSum310, AbbSum377, AbbSum480
	double complex AbbSum1833, AbbSum316, AbbSum856, AbbSum1309
	double complex AbbSum1621, AbbSum1843, AbbSum332, AbbSum1141
	double complex AbbSum1666, AbbSum1429, AbbSum1204, AbbSum945
	double complex AbbSum1610, AbbSum1123, AbbSum1494, AbbSum1980
	double complex AbbSum615, AbbSum1981, AbbSum723, AbbSum500
	double complex AbbSum1130, AbbSum1692, AbbSum1689, AbbSum1501
	double complex AbbSum1075, AbbSum1532, AbbSum1780, AbbSum1922
	double complex AbbSum652, AbbSum581, AbbSum948, AbbSum1613
	double complex AbbSum1699, AbbSum1207, AbbSum1431, AbbSum1684
	double complex AbbSum1866, AbbSum1846, AbbSum1932, AbbSum1759
	double complex AbbSum582, AbbSum653, AbbSum654, AbbSum583
	double complex AbbSum986, AbbSum1412, AbbSum526, AbbSum701
	double complex AbbSum528, AbbSum655, AbbSum584, AbbSum527
	double complex AbbSum470, AbbSum548, AbbSum699, AbbSum525
	double complex AbbSum1456, AbbSum1166, AbbSum565, AbbSum722
	double complex AbbSum1660, AbbSum673, AbbSum499, AbbSum1715
	double complex AbbSum1688, AbbSum1691, AbbSum2008, AbbSum1979
	double complex AbbSum1928, AbbSum721, AbbSum1765, AbbSum498
	double complex AbbSum501, AbbSum874, AbbSum873, AbbSum512
	double complex AbbSum250, AbbSum271, AbbSum248, AbbSum269
	double complex AbbSum18, AbbSum1776, AbbSum2081, AbbSum1045
	double complex AbbSum17, AbbSum1046, AbbSum2083, AbbSum823
	double complex AbbSum1042, AbbSum1110, AbbSum2082, AbbSum1043
	double complex AbbSum1047, AbbSum510, AbbSum729, AbbSum514
	double complex AbbSum692, AbbSum2213, AbbSum2042, AbbSum1321
	double complex AbbSum867, AbbSum524, AbbSum698, AbbSum529
	double complex AbbSum535, AbbSum705, AbbSum1603, AbbSum938
	double complex AbbSum1227, AbbSum1399, AbbSum1335, AbbSum832
	double complex AbbSum1226, AbbSum1336, AbbSum833, AbbSum834
	double complex AbbSum1337, AbbSum2185, AbbSum1290, AbbSum881
	double complex AbbSum1114, AbbSum1486, AbbSum1551, AbbSum1006
	double complex AbbSum1565, AbbSum1819, AbbSum1027, AbbSum1891
	double complex AbbSum1455, AbbSum1165, AbbSum1025, AbbSum909
	double complex AbbSum2118, AbbSum1277, AbbSum2293, AbbSum2116
	double complex AbbSum1821, AbbSum1563, AbbSum1893, AbbSum1023
	double complex AbbSum1026, AbbSum2292, AbbSum1564, AbbSum1820
	double complex AbbSum2117, AbbSum1024, AbbSum1892, AbbSum580
	double complex AbbSum651, AbbSum659, AbbSum2123, AbbSum590
	double complex AbbSum2296, AbbSum616, AbbSum1870, AbbSum1419
	double complex AbbSum1829, AbbSum1194, AbbSum622, AbbSum644
	double complex AbbSum786, AbbSum1144, AbbSum1438, AbbSum1919
	double complex AbbSum1777, AbbSum1547, AbbSum1000, AbbSum397
	double complex AbbSum272, AbbSum1224, AbbSum1397, AbbSum1924
	double complex AbbSum708, AbbSum1761, AbbSum475, AbbSum883
	double complex AbbSum1292, AbbSum405, AbbSum252, AbbSum1334
	double complex AbbSum831, AbbSum880, AbbSum557, AbbSum665
	double complex AbbSum1289, AbbSum1333, AbbSum830, AbbSum963
	double complex AbbSum882, AbbSum2084, AbbSum1291, AbbSum2300
	double complex AbbSum1581, AbbSum1883, AbbSum386, AbbSum290
	double complex AbbSum1806, AbbSum84, AbbSum156, AbbSum682
	double complex AbbSum544, AbbSum1398, AbbSum1225, AbbSum1049
	double complex AbbSum1514, AbbSum1658, AbbSum373, AbbSum1484
	double complex AbbSum1714, AbbSum305, AbbSum1112, AbbSum2301
	double complex AbbSum187, AbbSum884, AbbSum20, AbbSum2085
	double complex AbbSum1293, AbbSum1053, AbbSum1518, AbbSum1052
	double complex AbbSum1517, AbbSum1051, AbbSum1516, AbbSum1050
	double complex AbbSum1515, AbbSum2023, AbbSum2232, AbbSum795
	double complex AbbSum1380, AbbSum936, AbbSum1267, AbbSum1349
	double complex AbbSum851, AbbSum885, AbbSum2014, AbbSum1988
	double complex AbbSum1136, AbbSum1505, AbbSum899, AbbSum1306
	double complex AbbSum927, AbbSum981, AbbSum2176, AbbSum1601
	double complex AbbSum2248, AbbSum992, AbbSum1028, AbbSum1054
	double complex AbbSum1796, AbbSum1073, AbbSum1088, AbbSum1539
	double complex AbbSum1090, AbbSum1092, AbbSum1540, AbbSum1089
	double complex AbbSum1713, AbbSum1674, AbbSum2276, AbbSum2139
	double complex AbbSum1541, AbbSum1091, AbbSum1093, AbbSum2212
	double complex AbbSum2052, AbbSum1771, AbbSum1470, AbbSum1180
	double complex AbbSum1472, AbbSum1182, AbbSum2159, AbbSum2158
	double complex AbbSum1184, AbbSum1183, AbbSum1474, AbbSum1473
	double complex AbbSum1471, AbbSum1181, AbbSum1185, AbbSum2073
	double complex AbbSum1247, AbbSum2198, AbbSum1415, AbbSum1648
	double complex AbbSum568, AbbSum1656, AbbSum1731, AbbSum1982
	double complex AbbSum1999, AbbSum2316, AbbSum2182, AbbSum1994
	double complex AbbSum2224, AbbSum2037, AbbSum2119, AbbSum206
	double complex AbbSum836, AbbSum2264, AbbSum2266, AbbSum2155
	double complex AbbSum2157, AbbSum2265, AbbSum2156, AbbSum2160
	double complex AbbSum88, AbbSum915, AbbSum2061, AbbSum483
	double complex AbbSum1149, AbbSum973, AbbSum1675, AbbSum1701
	double complex AbbSum926, AbbSum1259, AbbSum1238, AbbSum2170
	double complex AbbSum1407, AbbSum2242, AbbSum1594, AbbSum980
	double complex AbbSum465, AbbSum2094, AbbSum737, AbbSum2186
	double complex AbbSum1724, AbbSum1647, AbbSum760, AbbSum1384
	double complex AbbSum1274, AbbSum906, AbbSum1058, AbbSum1522
	double complex AbbSum1897, AbbSum1782, AbbSum923, AbbSum924
	double complex AbbSum1373, AbbSum785, AbbSum976, AbbSum1592
	double complex AbbSum1721, AbbSum1643, AbbSum1258, AbbSum925
	double complex AbbSum1608, AbbSum943, AbbSum1738, AbbSum1625
	double complex AbbSum1966, AbbSum2018, AbbSum447, AbbSum748
	double complex AbbSum336, AbbSum351, AbbSum2029, AbbSum2216
	double complex AbbSum563, AbbSum671, AbbSum91, AbbSum160
	double complex AbbSum1295, AbbSum837, AbbSum1697, AbbSum601
	double complex AbbSum1959, AbbSum1201, AbbSum1152, AbbSum1836
	double complex AbbSum320, AbbSum484, AbbSum1153, AbbSum2193
	double complex AbbSum191, AbbSum1681, AbbSum628, AbbSum1426
	double complex AbbSum1444, AbbSum1850, AbbSum363, AbbSum711
	double complex AbbSum1876, AbbSum1445, AbbSum2066, AbbSum28
	double complex AbbSum783, AbbSum782, AbbSum1372, AbbSum1257
	double complex AbbSum922, AbbSum1371, AbbSum781, AbbSum784
	double complex AbbSum1154, AbbSum1593, AbbSum977, AbbSum1758
	double complex AbbSum647, AbbSum1931, AbbSum570, AbbSum1939
	double complex AbbSum1753, AbbSum208, AbbSum428, AbbSum1644
	double complex AbbSum230, AbbSum416, AbbSum979, AbbSum978
	double complex AbbSum1646, AbbSum1723, AbbSum1443, AbbSum1151
	double complex AbbSum1722, AbbSum1645, AbbSum2, AbbSum571
	double complex AbbSum1661, AbbSum521, AbbSum572, AbbSum97
	double complex AbbSum1811, AbbSum574, AbbSum308, AbbSum573
	double complex AbbSum575, AbbSum516, AbbSum2115, AbbSum1813
	double complex AbbSum1814, AbbSum1164, AbbSum235, AbbSum1018
	double complex AbbSum275, AbbSum503, AbbSum57, AbbSum517
	double complex AbbSum276, AbbSum520, AbbSum101, AbbSum65
	double complex AbbSum31, AbbSum2068, AbbSum2131, AbbSum2132
	double complex AbbSum2133, AbbSum279, AbbSum100, AbbSum776
	double complex AbbSum778, AbbSum335, AbbSum229, AbbSum1752
	double complex AbbSum1844, AbbSum779, AbbSum916, AbbSum780
	double complex AbbSum914, AbbSum1635, AbbSum26, AbbSum2075
	double complex AbbSum1640, AbbSum975, AbbSum913, AbbSum27
	double complex AbbSum1237, AbbSum482, AbbSum1641, AbbSum2092
	double complex AbbSum1148, AbbSum1639, AbbSum1134, AbbSum1504
	double complex AbbSum1703, AbbSum1677, AbbSum1135, AbbSum1690
	double complex AbbSum1693, AbbSum73, AbbSum170, AbbSum1993
	double complex AbbSum2181, AbbSum2315, AbbSum1998, AbbSum1824
	double complex AbbSum1896, AbbSum1930, AbbSum1769, AbbSum1662
	double complex AbbSum1716, AbbSum533, AbbSum549, AbbSum689
	double complex AbbSum534, AbbSum1673, AbbSum586, AbbSum311
	double complex AbbSum378, AbbSum588, AbbSum1729, AbbSum1654
	double complex AbbSum739, AbbSum471, AbbSum1600, AbbSum989
	double complex AbbSum1781, AbbSum1923, AbbSum280, AbbSum400
	double complex AbbSum1627, AbbSum3, AbbSum1118, AbbSum794
	double complex AbbSum1278, AbbSum910, AbbSum1133, AbbSum764
	double complex AbbSum935, AbbSum1131, AbbSum1502, AbbSum2071
	double complex AbbSum1265, AbbSum933, AbbSum1801, AbbSum1913
	double complex AbbSum2021, AbbSum2230, AbbSum792, AbbSum1378
	double complex AbbSum1895, AbbSum1823, AbbSum2135, AbbSum2274
	double complex AbbSum2137, AbbSum2136, AbbSum1086, AbbSum1084
	double complex AbbSum1087, AbbSum1083, AbbSum1081, AbbSum1082
	double complex AbbSum144, AbbSum104, AbbSum1132, AbbSum1503
	double complex AbbSum793, AbbSum1379, AbbSum2022, AbbSum2231
	double complex AbbSum1538, AbbSum1085, AbbSum1266, AbbSum934
	double complex AbbSum658, AbbSum589, AbbSum704, AbbSum532
	double complex AbbSum1672, AbbSum1712, AbbSum2048, AbbSum2208
	double complex AbbSum1537, AbbSum1080, AbbSum1760, AbbSum1933
	double complex AbbSum991, AbbSum990, AbbSum1730, AbbSum1655
	double complex AbbSum146, AbbSum106, AbbSum587, AbbSum657
	double complex AbbSum1599, AbbSum988, AbbSum145, AbbSum105
	double complex AbbSum1243, AbbSum1413, AbbSum2245, AbbSum2173
	double complex AbbSum987, AbbSum1598, AbbSum508, AbbSum727
	double complex AbbSum1986, AbbSum2012, AbbSum506, AbbSum725
	double complex AbbSum2010, AbbSum1984, AbbSum2049, AbbSum2209
	double complex AbbSum216, AbbSum433, AbbSum1757, AbbSum1943
	double complex AbbSum898, AbbSum1305, AbbSum643, AbbSum621
	double complex AbbSum2261, AbbSum1177, AbbSum2152, AbbSum1467
	double complex AbbSum2262, AbbSum2153, AbbSum1178, AbbSum1468
	double complex AbbSum1847, AbbSum1867, AbbSum1245, AbbSum1246
	double complex AbbSum2175, AbbSum2247, AbbSum1414, AbbSum2174
	double complex AbbSum1244, AbbSum2246, AbbSum1948, AbbSum1862
	double complex AbbSum1839, AbbSum1879, AbbSum367, AbbSum325
	double complex AbbSum2210, AbbSum2050, AbbSum1348, AbbSum850
	double complex AbbSum2036, AbbSum2223, AbbSum2011, AbbSum726
	double complex AbbSum1985, AbbSum507, AbbSum1911, AbbSum1799
	double complex AbbSum1651, AbbSum1727, AbbSum1076, AbbSum1533
	double complex AbbSum2134, AbbSum2273, AbbSum1534, AbbSum1077
	double complex AbbSum2020, AbbSum2229, AbbSum790, AbbSum1377
	double complex AbbSum2069, AbbSum2195, AbbSum931, AbbSum1263
	double complex AbbSum1800, AbbSum1912, AbbSum1535, AbbSum1078
	double complex AbbSum1894, AbbSum1822, AbbSum142, AbbSum102
	double complex AbbSum2122, AbbSum1568, AbbSum1032, AbbSum1802
	double complex AbbSum1914, AbbSum181, AbbSum62, AbbSum868
	double complex AbbSum123, AbbSum139, AbbSum1322, AbbSum1364
	double complex AbbSum813, AbbSum390, AbbSum296, AbbSum1854
	double complex AbbSum1954, AbbSum409, AbbSum259, AbbSum1853
	double complex AbbSum1956, AbbSum1687, AbbSum1694, AbbSum1567
	double complex AbbSum1030, AbbSum2295, AbbSum2121, AbbSum2275
	double complex AbbSum2138, AbbSum1079, AbbSum1536, AbbSum703
	double complex AbbSum531, AbbSum702, AbbSum530, AbbSum1710
	double complex AbbSum1670, AbbSum1711, AbbSum1671, AbbSum2009
	double complex AbbSum724, AbbSum505, AbbSum1983, AbbSum585
	double complex AbbSum656, AbbSum2096, AbbSum2310, AbbSum1987
	double complex AbbSum509, AbbSum728, AbbSum2013, AbbSum143
	double complex AbbSum103, AbbSum1652, AbbSum1168, AbbSum2143
	double complex AbbSum2252, AbbSum1458, AbbSum1728, AbbSum1653
	double complex AbbSum2072, AbbSum2197, AbbSum1031, AbbSum949
	double complex AbbSum1614, AbbSum1700, AbbSum1686, AbbSum1960
	double complex AbbSum1851, AbbSum1432, AbbSum1208, AbbSum1967
	double complex AbbSum2019, AbbSum751, AbbSum451, AbbSum675
	double complex AbbSum566, AbbSum94, AbbSum161, AbbSum32
	double complex AbbSum192, AbbSum641, AbbSum617, AbbSum791
	double complex AbbSum238, AbbSum421, AbbSum894, AbbSum1301
	double complex AbbSum2294, AbbSum2120, AbbSum1029, AbbSum1566
	double complex AbbSum2051, AbbSum2211, AbbSum2263, AbbSum1179
	double complex AbbSum1469, AbbSum2154, AbbSum2070, AbbSum2196
	double complex AbbSum932, AbbSum1264, AbbSum1745, AbbSum1634
	double complex AbbSum1946, AbbSum1865, AbbSum343, AbbSum356
	double complex AbbSum1344, AbbSum846, AbbSum2219, AbbSum2032
	double complex AbbSum2144, AbbSum1169, AbbSum2253, AbbSum1459
	double complex AbbSum68, AbbSum166, AbbSum1989, AbbSum2177
	double complex AbbSum2311, AbbSum1995, AbbSum70, AbbSum1991
	double complex AbbSum2179, AbbSum2313, AbbSum1996, AbbSum1990
	double complex AbbSum69, AbbSum2178, AbbSum2312, AbbSum167
	double complex AbbSum1170, AbbSum2145, AbbSum2254, AbbSum1460
	double complex AbbSum1174, AbbSum2149, AbbSum2258, AbbSum1464
	double complex AbbSum642, AbbSum618, AbbSum1172, AbbSum2147
	double complex AbbSum2256, AbbSum1462, AbbSum619, AbbSum620
	double complex AbbSum1171, AbbSum2146, AbbSum2255, AbbSum1461
	double complex AbbSum2207, AbbSum2047, AbbSum896, AbbSum1303
	double complex AbbSum895, AbbSum1302, AbbSum1175, AbbSum1465
	double complex AbbSum2259, AbbSum2150, AbbSum897, AbbSum1304
	double complex AbbSum2257, AbbSum1173, AbbSum2148, AbbSum1463
	double complex AbbSum1346, AbbSum848, AbbSum2221, AbbSum2034
	double complex AbbSum1347, AbbSum849, AbbSum2222, AbbSum2035
	double complex AbbSum847, AbbSum1345, AbbSum2033, AbbSum2220
	double complex AbbSum71, AbbSum168, AbbSum1992, AbbSum2180
	double complex AbbSum2314, AbbSum1997, AbbSum2260, AbbSum1176
	double complex AbbSum1466, AbbSum2151, AbbSum592, AbbSum593
	double complex AbbSum1138, AbbSum282, AbbSum245, AbbSum442
	double complex AbbSum2101, AbbSum301, AbbSum298, AbbSum994
	double complex AbbSum476, AbbSum912, AbbSum1230, AbbSum477
	double complex AbbSum197, AbbSum313, AbbSum200, AbbSum479
	double complex AbbSum1003, AbbSum903, AbbSum2074, AbbSum1842
	double complex AbbSum328, AbbSum330, AbbSum610, AbbSum2077
	double complex AbbSum213, AbbSum1976, AbbSum2041, AbbSum1845
	double complex AbbSum45, AbbSum2039, AbbSum2076, AbbSum1659
	double complex AbbSum96, AbbSum108, AbbSum46, AbbSum2040
	double complex AbbSum49, AbbSum1011, AbbSum112, AbbSum53
	double complex AbbSum1016, AbbSum841, AbbSum608, AbbSum606
	double complex AbbSum607, AbbSum1059, AbbSum2045, AbbSum1783
	double complex AbbSum1061, AbbSum1062, AbbSum1012, AbbSum1063
	double complex AbbSum116, AbbSum494, AbbSum2129, AbbSum1789
	double complex AbbSum1974, AbbSum1975, AbbSum983, AbbSum255
	double complex AbbSum488, AbbSum52, AbbSum256, AbbSum293
	double complex AbbSum292, AbbSum1767, AbbSum1906, AbbSum1793
	double complex AbbSum1262, AbbSum1529, AbbSum1069, AbbSum930
	double complex AbbSum2194, AbbSum2067, AbbSum117, AbbSum865
	double complex AbbSum1319, AbbSum134, AbbSum1070, AbbSum1794
	double complex AbbSum1907, AbbSum1122, AbbSum1493, AbbSum808
	double complex AbbSum1359, AbbSum389, AbbSum295, AbbSum258
	double complex AbbSum407, AbbSum366, AbbSum324, AbbSum984
	double complex AbbSum1596, AbbSum2007, AbbSum1942, AbbSum2006
	double complex AbbSum1978, AbbSum1756, AbbSum1977, AbbSum1709
	double complex AbbSum1448, AbbSum1669, AbbSum1157, AbbSum1298
	double complex AbbSum889, AbbSum1340, AbbSum2218, AbbSum842
	double complex AbbSum2031, AbbSum55, AbbSum2027, AbbSum1947
	double complex AbbSum1878, AbbSum496, AbbSum1017, AbbSum1852
	double complex AbbSum497, AbbSum2206, AbbSum176, AbbSum2228
	double complex AbbSum1861, AbbSum1838, AbbSum719, AbbSum720
	double complex AbbSum1955, AbbSum1559, AbbSum2046, AbbSum1905
	double complex AbbSum1792, AbbSum1528, AbbSum1068, AbbSum1702
	double complex AbbSum1676, AbbSum750, AbbSum449, AbbSum614
	double complex AbbSum214, AbbSum432, AbbSum1241, AbbSum1410
	double complex AbbSum1072, AbbSum2172, AbbSum2244, AbbSum1530
	double complex AbbSum1071, AbbSum638, AbbSum613, AbbSum1795
	double complex AbbSum611, AbbSum1908, AbbSum636, AbbSum2272
	double complex AbbSum2130, AbbSum612, AbbSum1791, AbbSum637
	double complex AbbSum1904, AbbSum396, AbbSum681, AbbSum1602
	double complex AbbSum81, AbbSum82, AbbSum1193, AbbSum1580
	double complex AbbSum597, AbbSum1579, AbbSum372, AbbSum1918
	double complex AbbSum80, AbbSum1577, AbbSum1576, AbbSum153
	double complex AbbSum959, AbbSum958, AbbSum304, AbbSum1775
	double complex AbbSum543, AbbSum937, AbbSum154, AbbSum155
	double complex AbbSum1418, AbbSum962, AbbSum445, AbbSum961
	double complex AbbSum270, AbbSum37, AbbSum1221, AbbSum1307
	double complex AbbSum1475, AbbSum1324, AbbSum819, AbbSum124
	double complex AbbSum1394, AbbSum854, AbbSum1101, AbbSum221
	double complex AbbSum1420, AbbSum1584, AbbSum1582, AbbSum1585
	double complex AbbSum967, AbbSum964, AbbSum966, AbbSum326
	double complex AbbSum1195, AbbSum1520, AbbSum683, AbbSum757
	double complex AbbSum21, AbbSum902, AbbSum2163, AbbSum1937
	double complex AbbSum558, AbbSum1251, AbbSum911, AbbSum1056
	double complex AbbSum461, AbbSum1381, AbbSum85, AbbSum1271
	double complex AbbSum2235, AbbSum1749, AbbSum666, AbbSum22
	double complex AbbSum1732, AbbSum758, AbbSum1382, AbbSum86
	double complex AbbSum1617, AbbSum1485, AbbSum189, AbbSum2282
	double complex AbbSum464, AbbSum560, AbbSum1005, AbbSum600
	double complex AbbSum1405, AbbSum188, AbbSum190, AbbSum2192
	double complex AbbSum599, AbbSum1294, AbbSum2065, AbbSum626
	double complex AbbSum835, AbbSum1113, AbbSum24, AbbSum2105
	double complex AbbSum736, AbbSum668, AbbSum1550, AbbSum627
	double complex AbbSum1235, AbbSum23, AbbSum25, AbbSum411
	double complex AbbSum141, AbbSum1421, AbbSum1004, AbbSum273
	double complex AbbSum970, AbbSum223, AbbSum99, AbbSum1196
	double complex AbbSum1588, AbbSum1549, AbbSum374, AbbSum344
	double complex AbbSum816, AbbSum76, AbbSum539, AbbSum1189
	double complex AbbSum1416, AbbSum218, AbbSum1280, AbbSum152
	double complex AbbSum676, AbbSum1527, AbbSum928, AbbSum133
	double complex AbbSum1428, AbbSum864, AbbSum944, AbbSum1628
	double complex AbbSum1297, AbbSum1060, AbbSum1296, AbbSum131
	double complex AbbSum863, AbbSum2288, AbbSum1595, AbbSum132
	double complex AbbSum862, AbbSum2287, AbbSum114, AbbSum1316
	double complex AbbSum2110, AbbSum982, AbbSum1260, AbbSum115
	double complex AbbSum1203, AbbSum1318, AbbSum1609, AbbSum1740
	double complex AbbSum888, AbbSum1523, AbbSum887, AbbSum113
	double complex AbbSum1317, AbbSum2111, AbbSum1067, AbbSum1558
	double complex AbbSum406, AbbSum294, AbbSum839, AbbSum840
	double complex AbbSum175, AbbSum173, AbbSum1557, AbbSum2109
	double complex AbbSum1556, AbbSum174, AbbSum1358, AbbSum2268
	double complex AbbSum1013, AbbSum51, AbbSum807, AbbSum2125
	double complex AbbSum2286, AbbSum1014, AbbSum1338, AbbSum1339
	double complex AbbSum54, AbbSum50, AbbSum1015, AbbSum257
	double complex AbbSum388, AbbSum714, AbbSum1971, AbbSum1972
	double complex AbbSum1926, AbbSum1927, AbbSum2030, AbbSum2205
	double complex AbbSum1925, AbbSum1492, AbbSum1385, AbbSum605
	double complex AbbSum1409, AbbSum1969, AbbSum1968, AbbSum1970
	double complex AbbSum467, AbbSum604, AbbSum431, AbbSum430
	double complex AbbSum2203, AbbSum1275, AbbSum929, AbbSum1766
	double complex AbbSum564, AbbSum29, AbbSum718, AbbSum749
	double complex AbbSum716, AbbSum717, AbbSum1755, AbbSum1155
	double complex AbbSum365, AbbSum1707, AbbSum1156, AbbSum1708
	double complex AbbSum713, AbbSum2187, AbbSum1787, AbbSum603
	double complex AbbSum1526, AbbSum1525, AbbSum2269, AbbSum712
	double complex AbbSum1786, AbbSum602, AbbSum1408, AbbSum429
	double complex AbbSum1120, AbbSum2270, AbbSum1754, AbbSum1877
	double complex AbbSum364, AbbSum2171, AbbSum1784, AbbSum1785
	double complex AbbSum1788, AbbSum1790, AbbSum1973, AbbSum609
	double complex AbbSum715, AbbSum1524, AbbSum2271, AbbSum2128
	double complex AbbSum1064, AbbSum2005, AbbSum635, AbbSum491
	double complex AbbSum1239, AbbSum210, AbbSum1903, AbbSum1902
	double complex AbbSum1491, AbbSum2127, AbbSum1940, AbbSum1837
	double complex AbbSum322, AbbSum2243, AbbSum1898, AbbSum1899
	double complex AbbSum631, AbbSum1066, AbbSum487, AbbSum1900
	double complex AbbSum1065, AbbSum2126, AbbSum632, AbbSum1901
	double complex AbbSum1121, AbbSum761, AbbSum634, AbbSum1240
	double complex AbbSum2001, AbbSum2000, AbbSum2002, AbbSum687
	double complex AbbSum633, AbbSum212, AbbSum211, AbbSum2059
	double complex AbbSum907, AbbSum1261, AbbSum672, AbbSum92
	double complex AbbSum495, AbbSum448, AbbSum492, AbbSum493
	double complex AbbSum1941, AbbSum1446, AbbSum323, AbbSum1667
	double complex AbbSum1447, AbbSum1668, AbbSum489, AbbSum1764
	double complex AbbSum2003, AbbSum2004, AbbSum2217, AbbSum2044
	double complex AbbSum1762, AbbSum1763, AbbSum490, AbbSum518
	double complex AbbSum890, AbbSum700, AbbSum1021, AbbSum375
	double complex AbbSum376, AbbSum30, AbbSum674, AbbSum1560
	double complex AbbSum1778, AbbSum1612, AbbSum1779, AbbSum519
	double complex AbbSum1562, AbbSum1611, AbbSum1863, AbbSum946
	double complex AbbSum1944, AbbSum1022, AbbSum469, AbbSum1561
	double complex AbbSum307, AbbSum309, AbbSum93, AbbSum450
	double complex AbbSum1019, AbbSum1920, AbbSum947, AbbSum1921
	double complex AbbSum695, AbbSum694, AbbSum1299, AbbSum1282
	double complex AbbSum286, AbbSum1103, AbbSum299, AbbSum1859
	double complex AbbSum870, AbbSum247, AbbSum1478, AbbSum369
	double complex AbbSum1953, AbbSum393, AbbSum63, AbbSum1773
	double complex AbbSum125, AbbSum1858, AbbSum1281, AbbSum1951
	double complex AbbSum1950, AbbSum511, AbbSum1949, AbbSum732
	double complex AbbSum1350, AbbSum798, AbbSum1855, AbbSum459
	double complex AbbSum266, AbbSum162, AbbSum1916, AbbSum690
	double complex AbbSum1856, AbbSum1857, AbbSum38, AbbSum1952
	double complex AbbSum869, AbbSum1519, AbbSum127, AbbSum1586
	double complex AbbSum1002, AbbSum478, AbbSum64, AbbSum2164
	double complex AbbSum387, AbbSum1873, AbbSum361, AbbSum201
	double complex AbbSum2304, AbbSum1936, AbbSum965, AbbSum2162
	double complex AbbSum1934, AbbSum199, AbbSum1402, AbbSum1228
	double complex AbbSum198, AbbSum1935, AbbSum1831, AbbSum312
	double complex AbbSum1830, AbbSum314, AbbSum1229, AbbSum2302
	double complex AbbSum2303, AbbSum2161, AbbSum2233, AbbSum2086
	double complex AbbSum2087, AbbSum1747, AbbSum1401, AbbSum1872
	double complex AbbSum359, AbbSum1871, AbbSum360, AbbSum1583
	double complex AbbSum1746, AbbSum423, AbbSum1231, AbbSum1400
	double complex AbbSum422, AbbSum2234, AbbSum1748, AbbSum1055
	double complex AbbSum41, AbbSum968, AbbSum1548, AbbSum163
	double complex AbbSum709, AbbSum2088, AbbSum2236, AbbSum253
	double complex AbbSum424, AbbSum1832, AbbSum315, AbbSum1618
	double complex AbbSum203, AbbSum225, AbbSum1885, AbbSum2089
	double complex AbbSum1619, AbbSum645, AbbSum128, AbbSum684
	double complex AbbSum222, AbbSum1807, AbbSum969, AbbSum347
	double complex AbbSum462, AbbSum1423, AbbSum224, AbbSum1422
	double complex AbbSum1272, AbbSum1365, AbbSum327, AbbSum1605
	double complex AbbSum331, AbbSum226, AbbSum1620, AbbSum1604
	double complex AbbSum1735, AbbSum939, AbbSum1197, AbbSum904
	double complex AbbSum771, AbbSum414, AbbSum346, AbbSum940
	double complex AbbSum348, AbbSum1198, AbbSum412, AbbSum317
	double complex AbbSum413, AbbSum1808, AbbSum2305, AbbSum1733
	double complex AbbSum1809, AbbSum42, AbbSum545, AbbSum410
	double complex AbbSum1884, AbbSum734, AbbSum329, AbbSum1587
	double complex AbbSum567, AbbSum1734, AbbSum920, AbbSum773
	double complex AbbSum772, AbbSum775, AbbSum129, AbbSum921
	double complex AbbSum801, AbbSum569, AbbSum1150, AbbSum1834
	double complex AbbSum1835, AbbSum1200, AbbSum1607, AbbSum334
	double complex AbbSum774, AbbSum227, AbbSum349, AbbSum1848
	double complex AbbSum941, AbbSum362, AbbSum1886, AbbSum2090
	double complex AbbSum2091, AbbSum710, AbbSum1146, AbbSum1145
	double complex AbbSum1636, AbbSum1370, AbbSum1199, AbbSum228
	double complex AbbSum1849, AbbSum2093, AbbSum2308, AbbSum1424
	double complex AbbSum415, AbbSum1958, AbbSum1717, AbbSum777
	double complex AbbSum481, AbbSum1440, AbbSum1439, AbbSum1369
	double complex AbbSum43, AbbSum1256, AbbSum1310, AbbSum646
	double complex AbbSum1442, AbbSum1874, AbbSum1875, AbbSum1425
	double complex AbbSum942, AbbSum350, AbbSum1368, AbbSum318
	double complex AbbSum333, AbbSum1957, AbbSum1606, AbbSum319
	double complex AbbSum2306, AbbSum2307, AbbSum1810, AbbSum1255
	double complex AbbSum1367, AbbSum1366, AbbSum1314, AbbSum838
	double complex AbbSum804, AbbSum1312, AbbSum806, AbbSum861
	double complex AbbSum803, AbbSum1315, AbbSum1354, AbbSum1357
	double complex AbbSum886, AbbSum1355, AbbSum858, AbbSum860
	double complex AbbSum149, AbbSum1571, AbbSum2053, AbbSum1268
	double complex AbbSum150, AbbSum540, AbbSum662, AbbSum1868
	double complex AbbSum357, AbbSum151, AbbSum1574, AbbSum148
	double complex AbbSum7, AbbSum955, AbbSum2199, AbbSum752
	double complex AbbSum8, AbbSum677, AbbSum554, AbbSum1827
	double complex AbbSum194, AbbSum75, AbbSum9, AbbSum953
	double complex AbbSum1487, AbbSum209, AbbSum630, AbbSum337
	double complex AbbSum485, AbbSum1739, AbbSum1202, AbbSum1554
	double complex AbbSum1116, AbbSum1488, AbbSum1009, AbbSum321
	double complex AbbSum486, AbbSum231, AbbSum629, AbbSum1626
	double complex AbbSum1427, AbbSum1115, AbbSum291, AbbSum1490
	double complex AbbSum1311, AbbSum1356, AbbSum274, AbbSum2202
	double complex AbbSum802, AbbSum1313, AbbSum686, AbbSum2285
	double complex AbbSum1489, AbbSum2024, AbbSum130, AbbSum44
	double complex AbbSum1008, AbbSum2025, AbbSum47, AbbSum1010
	double complex AbbSum1007, AbbSum2284, AbbSum48, AbbSum2026
	double complex AbbSum2283, AbbSum2106, AbbSum111, AbbSum2227
	double complex AbbSum2107, AbbSum1552, AbbSum2226, AbbSum110
	double complex AbbSum1555, AbbSum1553, AbbSum254, AbbSum1119
	double complex AbbSum857, AbbSum805, AbbSum306, AbbSum2058
	double complex AbbSum1353, AbbSum859, AbbSum1117, AbbSum2225
	double complex AbbSum109, AbbSum172, AbbSum2108, AbbSum466
	double complex AbbSum1818, AbbSum1375, AbbSum179, AbbSum1362
	double complex AbbSum1450, AbbSum1452, AbbSum178, AbbSum1361
	double complex AbbSum1454, AbbSum639, AbbSum1242, AbbSum1449
	double complex AbbSum1451, AbbSum1816, AbbSum985, AbbSum1374
	double complex AbbSum2095, AbbSum1682, AbbSum1430, AbbSum1817
	double complex AbbSum1815, AbbSum576, AbbSum578, AbbSum1360
	double complex AbbSum177, AbbSum577, AbbSum649, AbbSum648
	double complex AbbSum650, AbbSum809, AbbSum56, AbbSum1887
	double complex AbbSum1889, AbbSum1698, AbbSum1205, AbbSum788
	double complex AbbSum59, AbbSum811, AbbSum1159, AbbSum1161
	double complex AbbSum58, AbbSum810, AbbSum1163, AbbSum502
	double complex AbbSum1411, AbbSum1158, AbbSum1160, AbbSum1597
	double complex AbbSum787, AbbSum2309, AbbSum1888, AbbSum1890
	double complex AbbSum2249, AbbSum1633, AbbSum215, AbbSum408
	double complex AbbSum2250, AbbSum72, AbbSum1768, AbbSum1650
	double complex AbbSum1798, AbbSum164, AbbSum67, AbbSum2251
	double complex AbbSum1797, AbbSum1649, AbbSum1909, AbbSum1725
	double complex AbbSum169, AbbSum1929, AbbSum1726, AbbSum1910
	double complex AbbSum66, AbbSum165, AbbSum2140, AbbSum2141
	double complex AbbSum1685, AbbSum342, AbbSum237, AbbSum2142
	double complex AbbSum561, AbbSum1938, AbbSum2267, AbbSum427
	double complex AbbSum562, AbbSum905, AbbSum1622, AbbSum1591
	double complex AbbSum204, AbbSum1680, AbbSum205, AbbSum1678
	double complex AbbSum2201, AbbSum1737, AbbSum971, AbbSum972
	double complex AbbSum1720, AbbSum463, AbbSum685, AbbSum1057
	double complex AbbSum559, AbbSum446, AbbSum87, AbbSum90
	double complex AbbSum2241, AbbSum89, AbbSum1719, AbbSum2237
	double complex AbbSum1383, AbbSum1254, AbbSum1147, AbbSum2239
	double complex AbbSum1718, AbbSum1253, AbbSum2189, AbbSum2191
	double complex AbbSum1679, AbbSum2215, AbbSum1236, AbbSum1234
	double complex AbbSum2240, AbbSum2238, AbbSum1233, AbbSum1252
	double complex AbbSum2190, AbbSum2063, AbbSum1403, AbbSum917
	double complex AbbSum1404, AbbSum1406, AbbSum2168, AbbSum2166
	double complex AbbSum2064, AbbSum1696, AbbSum2028, AbbSum759
	double complex AbbSum919, AbbSum1637, AbbSum918, AbbSum2062
	double complex AbbSum2167, AbbSum1441, AbbSum2165, AbbSum207
	double complex AbbSum670, AbbSum1273, AbbSum1736, AbbSum974
	double complex AbbSum425, AbbSum1624, AbbSum426, AbbSum1695
	double complex AbbSum2057, AbbSum1623, AbbSum1589, AbbSum1590
	double complex AbbSum1642, AbbSum735, AbbSum546, AbbSum1521
	double complex AbbSum1638, AbbSum667, AbbSum747, AbbSum157
	double complex AbbSum159, AbbSum2169, AbbSum158, AbbSum2124
	double complex AbbSum1751, AbbSum669, AbbSum1125, AbbSum398
	double complex AbbSum137, AbbSum844, AbbSum136, AbbSum1531
	double complex AbbSum1126, AbbSum843, AbbSum1127, AbbSum763
	double complex AbbSum468, AbbSum2289, AbbSum418, AbbSum2290
	double complex AbbSum688, AbbSum908, AbbSum419, AbbSum353
	double complex AbbSum2060, AbbSum1630, AbbSum2114, AbbSum354
	double complex AbbSum1631, AbbSum1124, AbbSum866, AbbSum135
	double complex AbbSum696, AbbSum352, AbbSum1629, AbbSum417
	double complex AbbSum1162, AbbSum1499, AbbSum697, AbbSum693
	double complex AbbSum515, AbbSum1128, AbbSum523, AbbSum1320
	double complex AbbSum118, AbbSum522, AbbSum1453, AbbSum232
	double complex AbbSum338, AbbSum1741, AbbSum1495, AbbSum121
	double complex AbbSum1342, AbbSum119, AbbSum1074, AbbSum1497
	double complex AbbSum1341, AbbSum1498, AbbSum1387, AbbSum340
	double complex AbbSum1743, AbbSum2112, AbbSum233, AbbSum2113
	double complex AbbSum547, AbbSum1276, AbbSum234, AbbSum339
	double complex AbbSum2204, AbbSum1742, AbbSum2291, AbbSum738
	double complex AbbSum277, AbbSum1496, AbbSum552, AbbSum731
	double complex AbbSum368, AbbSum182, AbbSum553, AbbSum591
	double complex AbbSum742, AbbSum1573, AbbSum379, AbbSum743
	double complex AbbSum11, AbbSum730, AbbSum952, AbbSum2017
	double complex AbbSum1250, AbbSum740, AbbSum663, AbbSum741
	double complex AbbSum1390, AbbSum2016, AbbSum1249, AbbSum1389
	double complex AbbSum1279, AbbSum2015, AbbSum1248, AbbSum1570
	double complex AbbSum1391, AbbSum1388, AbbSum1214, AbbSum1217
	double complex AbbSum951, AbbSum814, AbbSum1961, AbbSum765
	double complex AbbSum1215, AbbSum1216, AbbSum555, AbbSum437
	double complex AbbSum1962, AbbSum767, AbbSum436, AbbSum456
	double complex AbbSum261, AbbSum10, AbbSum661, AbbSum472
	double complex AbbSum438, AbbSum240, AbbSum439, AbbSum183
	double complex AbbSum454, AbbSum1572, AbbSum1964, AbbSum768
	double complex AbbSum954, AbbSum660, AbbSum297, AbbSum817
	double complex AbbSum1506, AbbSum402, AbbSum140, AbbSum1476
	double complex AbbSum594, AbbSum284, AbbSum382, AbbSum263
	double complex AbbSum1323, AbbSum1034, AbbSum243, AbbSum623
	double complex AbbSum98, AbbSum1102, AbbSum1393, AbbSum458
	double complex AbbSum1507, AbbSum381, AbbSum358, AbbSum403
	double complex AbbSum1508, AbbSum13, AbbSum1037, AbbSum744
	double complex AbbSum380, AbbSum1664, AbbSum1140, AbbSum401
	double complex AbbSum1510, AbbSum1704, AbbSum1139, AbbSum1392
	double complex AbbSum1706, AbbSum1137, AbbSum1219, AbbSum1433
	double complex AbbSum1665, AbbSum242, AbbSum1434, AbbSum1663
	double complex AbbSum1038, AbbSum1220, AbbSum678, AbbSum1035
	double complex AbbSum283, AbbSum195, AbbSum78, AbbSum1509
	double complex AbbSum441, AbbSum1435, AbbSum281, AbbSum1705
	double complex AbbSum1036, AbbSum244, AbbSum541, AbbSum1774
	double complex AbbSum268, AbbSum14, AbbSum1105, AbbSum265
	double complex AbbSum1104, AbbSum1269, AbbSum1480, AbbSum1192
	double complex AbbSum219, AbbSum267, AbbSum300, AbbSum302
	double complex AbbSum1543, AbbSum1575, AbbSum1481, AbbSum2200
	double complex AbbSum2279, AbbSum2277, AbbSum2278, AbbSum1544
	double complex AbbSum1772, AbbSum264, AbbSum1542, AbbSum391
	double complex AbbSum993, AbbSum1915, AbbSum957, AbbSum2054
	double complex AbbSum1107, AbbSum2100, AbbSum2098, AbbSum2099
	double complex AbbSum996, AbbSum79, AbbSum1479, AbbSum392
	double complex AbbSum753, AbbSum1106, AbbSum1477, AbbSum995
	double complex AbbSum1417, AbbSum345, AbbSum394, AbbSum370
	double complex AbbSum371, AbbSum395, AbbSum1917, AbbSum679
	double complex AbbSum1222, AbbSum1882, AbbSum385, AbbSum754
	double complex AbbSum1109, AbbSum999, AbbSum1483, AbbSum460
	double complex AbbSum680, AbbSum2079, AbbSum2078, AbbSum998
	double complex AbbSum1545, AbbSum2298, AbbSum2297, AbbSum1395
	double complex AbbSum1805, AbbSum289, AbbSum1270, AbbSum1482
	double complex AbbSum1546, AbbSum1111, AbbSum733, AbbSum542
	double complex AbbSum1283, AbbSum249, AbbSum1284, AbbSum876
	double complex AbbSum39, AbbSum799, AbbSum828, AbbSum691
	double complex AbbSum1326, AbbSum1142, AbbSum1325, AbbSum1143
	double complex AbbSum1223, AbbSum383, AbbSum2103, AbbSum746
	double complex AbbSum875, AbbSum877, AbbSum1828, AbbSum186
	double complex AbbSum185, AbbSum706, AbbSum707, AbbSum624
	double complex AbbSum1803, AbbSum288, AbbSum598, AbbSum825
	double complex AbbSum1327, AbbSum1308, AbbSum827, AbbSum826
	double complex AbbSum2080, AbbSum960, AbbSum824, AbbSum745
	double complex AbbSum2104, AbbSum1288, AbbSum664, AbbSum184
	double complex AbbSum1804, AbbSum1512, AbbSum1513, AbbSum1511
	double complex AbbSum1040, AbbSum1041, AbbSum1044, AbbSum1328
	double complex AbbSum1881, AbbSum879, AbbSum556, AbbSum15
	double complex AbbSum2281, AbbSum443, AbbSum822, AbbSum855
	double complex AbbSum1331, AbbSum1330, AbbSum1578, AbbSum2299
	double complex AbbSum1329, AbbSum126, AbbSum1351, AbbSum1332
	double complex AbbSum513, AbbSum821, AbbSum1436, AbbSum820
	double complex AbbSum1437, AbbSum1396, AbbSum1880, AbbSum384
	double complex AbbSum625, AbbSum2280, AbbSum444, AbbSum596
	double complex AbbSum2183, AbbSum1285, AbbSum1287, AbbSum1869
	double complex AbbSum19, AbbSum16, AbbSum473, AbbSum474
	double complex AbbSum287, AbbSum1286, AbbSum404, AbbSum872
	double complex AbbSum2184, AbbSum871, AbbSum1094, AbbSum35
	double complex AbbSum1095, AbbSum2097, AbbSum36, AbbSum796
	double complex AbbSum797, AbbSum262, AbbSum1096, AbbSum1097
	double complex AbbSum1098, AbbSum1190, AbbSum1099, AbbSum1191
	double complex AbbSum241, AbbSum852, AbbSum853, AbbSum1615
	double complex AbbSum40, AbbSum1001, AbbSum1616, AbbSum1
	double complex AbbSum220, AbbSum900, AbbSum901, AbbSum2055
	double complex AbbSum2056, AbbSum755, AbbSum756, AbbSum196
	double complex AbbSum1825, AbbSum1209, AbbSum5, AbbSum1210
	double complex AbbSum434, AbbSum435, AbbSum239, AbbSum1186
	double complex AbbSum1211, AbbSum33, AbbSum1187, AbbSum1212
	double complex AbbSum1213, AbbSum1188, AbbSum193, AbbSum550
	double complex AbbSum551, AbbSum950, AbbSum6, AbbSum1840
	double complex AbbSum260, AbbSum536, AbbSum537, AbbSum34
	double complex AbbSum452, AbbSum453, AbbSum217
	common /ddabbrev/ Pair5, Pair6, Pair14, Pair13, Pair7, Pair11
	common /ddabbrev/ Pair3, Pair4, Pair12, Pair9, Pair2, Pair10
	common /ddabbrev/ Pair1, Pair8, Eps25, Eps29, Eps30, Eps31
	common /ddabbrev/ Eps32, Eps35, Eps21, Eps33, Eps18, Eps19
	common /ddabbrev/ Eps20, Eps27, Eps16, Eps14, Eps28, Eps34
	common /ddabbrev/ Eps15, Eps12, Eps26, Eps13, Eps22, Eps23, Eps7
	common /ddabbrev/ Eps24, Eps8, Eps9, Eps10, Eps5, Eps6, Eps11
	common /ddabbrev/ Eps2, Eps3, Eps4, Eps17, Eps1, Abb273, Abb286
	common /ddabbrev/ Abb16, Abb216, Abb230, Abb259, Abb159, Abb203
	common /ddabbrev/ Abb221, Abb239, Abb211, Abb192, Abb217, Abb200
	common /ddabbrev/ Abb160, Abb233, Abb206, Abb245, Abb268, Abb246
	common /ddabbrev/ Abb231, Abb290, Abb288, Abb262, Abb56, Abb223
	common /ddabbrev/ Abb234, Abb256, Abb92, Abb267, Abb107, Abb269
	common /ddabbrev/ Abb235, Abb247, Abb25, Abb18, Abb342, Abb294
	common /ddabbrev/ Abb20, Abb57, Abb58, Abb26, Abb97, Abb114
	common /ddabbrev/ Abb41, Abb59, Abb108, Abb115, Abb60, Abb36
	common /ddabbrev/ Abb295, Abb104, Abb343, Abb61, Abb54, Abb212
	common /ddabbrev/ Abb122, Abb224, Abb218, Abb219, Abb95, Abb335
	common /ddabbrev/ Abb55, Abb157, Abb283, Abb207, Abb163, Abb201
	common /ddabbrev/ Abb254, Abb123, Abb257, Abb270, Abb164, Abb195
	common /ddabbrev/ Abb196, Abb250, Abb236, Abb260, Abb1, Abb334
	common /ddabbrev/ Abb347, Abb318, Abb319, Abb336, Abb85, Abb384
	common /ddabbrev/ Abb385, Abb27, Abb21, Abb320, Abb105, Abb337
	common /ddabbrev/ Abb62, Abb274, Abb193, Abb278, Abb263, Abb63
	common /ddabbrev/ Abb225, Abb284, Abb248, Abb93, Abb213, Abb109
	common /ddabbrev/ Abb167, Abb161, Abb124, Abb28, Abb2, Abb307
	common /ddabbrev/ Abb376, Abb375, Abb316, Abb296, Abb24, Abb321
	common /ddabbrev/ Abb53, Abb297, Abb64, Abb86, Abb133, Abb360
	common /ddabbrev/ Abb304, Abb134, Abb29, Abb365, Abb162, Abb125
	common /ddabbrev/ Abb126, Abb197, Abb127, Abb3, Abb322, Abb204
	common /ddabbrev/ Abb128, Abb240, Abb165, Abb271, Abb191, Abb198
	common /ddabbrev/ Abb251, Abb208, Abb275, Abb291, Abb241, Abb226
	common /ddabbrev/ Abb129, Abb255, Abb121, Abb348, Abb356, Abb358
	common /ddabbrev/ Abb135, Abb338, Abb37, Abb323, Abb90, Abb136
	common /ddabbrev/ Abb388, Abb98, Abb116, Abb50, Abb137, Abb4
	common /ddabbrev/ Abb138, Abb209, Abb303, Abb139, Abb130, Abb7
	common /ddabbrev/ Abb317, Abb361, Abb305, Abb140, Abb339, Abb99
	common /ddabbrev/ Abb366, Abb65, Abb252, Abb141, Abb10, Abb142
	common /ddabbrev/ Abb66, Abb205, Abb67, Abb289, Abb42, Abb272
	common /ddabbrev/ Abb94, Abb199, Abb264, Abb43, Abb315, Abb100
	common /ddabbrev/ Abb214, Abb227, Abb68, Abb69, Abb87, Abb349
	common /ddabbrev/ Abb357, Abb344, Abb143, Abb351, Abb40, Abb324
	common /ddabbrev/ Abb91, Abb292, Abb144, Abb345, Abb215, Abb70
	common /ddabbrev/ Abb51, Abb145, Abb5, Abb350, Abb372, Abb368
	common /ddabbrev/ Abb146, Abb44, Abb110, Abb174, Abb71, Abb8
	common /ddabbrev/ Abb293, Abb362, Abb306, Abb147, Abb352, Abb279
	common /ddabbrev/ Abb367, Abb72, Abb281, Abb148, Abb11, Abb149
	common /ddabbrev/ Abb73, Abb166, Abb194, Abb202, Abb131, Abb168
	common /ddabbrev/ Abb52, Abb389, Abb6, Abb386, Abb117, Abb9
	common /ddabbrev/ Abb308, Abb74, Abb12, Abb45, Abb46, Abb158
	common /ddabbrev/ Abb132, Abb232, Abb261, Abb302, Abb265, Abb210
	common /ddabbrev/ Abb242, Abb101, Abb243, Abb118, Abb253, Abb249
	common /ddabbrev/ Abb170, Abb346, Abb175, Abb298, Abb38, Abb176
	common /ddabbrev/ Abb75, Abb96, Abb390, Abb13, Abb355, Abb377
	common /ddabbrev/ Abb325, Abb326, Abb30, Abb76, Abb84, Abb177
	common /ddabbrev/ Abb299, Abb77, Abb285, Abb22, Abb380, Abb34
	common /ddabbrev/ Abb378, Abb382, Abb150, Abb102, Abb276, Abb151
	common /ddabbrev/ Abb244, Abb363, Abb78, Abb35, Abb379, Abb381
	common /ddabbrev/ Abb383, Abb391, Abb152, Abb103, Abb277, Abb153
	common /ddabbrev/ Abb280, Abb364, Abb79, Abb17, Abb31, Abb266
	common /ddabbrev/ Abb222, Abb171, Abb119, Abb172, Abb258, Abb220
	common /ddabbrev/ Abb287, Abb111, Abb228, Abb112, Abb237, Abb32
	common /ddabbrev/ Abb14, Abb329, Abb314, Abb169, Abb340, Abb369
	common /ddabbrev/ Abb178, Abb353, Abb330, Abb310, Abb23, Abb179
	common /ddabbrev/ Abb106, Abb180, Abb80, Abb19, Abb341, Abb181
	common /ddabbrev/ Abb300, Abb313, Abb331, Abb309, Abb88, Abb182
	common /ddabbrev/ Abb387, Abb370, Abb183, Abb47, Abb312, Abb184
	common /ddabbrev/ Abb120, Abb113, Abb154, Abb89, Abb327, Abb185
	common /ddabbrev/ Abb354, Abb186, Abb332, Abb311, Abb229, Abb373
	common /ddabbrev/ Abb371, Abb81, Abb282, Abb374, Abb359, Abb39
	common /ddabbrev/ Abb301, Abb173, Abb82, Abb15, Abb328, Abb333
	common /ddabbrev/ Abb187, Abb188, Abb238, Abb83, Abb33, Abb48
	common /ddabbrev/ Abb189, Abb155, Abb49, Abb392, Abb190, Abb156
	common /ddabbrev/ AbbSum60, AbbSum120, AbbSum892, AbbSum891
	common /ddabbrev/ AbbSum1812, AbbSum1020, AbbSum579, AbbSum4
	common /ddabbrev/ AbbSum1683, AbbSum1206, AbbSum1167, AbbSum640
	common /ddabbrev/ AbbSum180, AbbSum1363, AbbSum1457, AbbSum504
	common /ddabbrev/ AbbSum812, AbbSum61, AbbSum455, AbbSum1841
	common /ddabbrev/ AbbSum1218, AbbSum815, AbbSum769, AbbSum1826
	common /ddabbrev/ AbbSum538, AbbSum1963, AbbSum766, AbbSum74
	common /ddabbrev/ AbbSum171, AbbSum77, AbbSum1965, AbbSum457
	common /ddabbrev/ AbbSum770, AbbSum956, AbbSum440, AbbSum12
	common /ddabbrev/ AbbSum1770, AbbSum1100, AbbSum1657, AbbSum95
	common /ddabbrev/ AbbSum1033, AbbSum1569, AbbSum147, AbbSum107
	common /ddabbrev/ AbbSum1129, AbbSum845, AbbSum138, AbbSum1343
	common /ddabbrev/ AbbSum1500, AbbSum122, AbbSum818, AbbSum2038
	common /ddabbrev/ AbbSum2043, AbbSum2214, AbbSum1750, AbbSum1232
	common /ddabbrev/ AbbSum202, AbbSum762, AbbSum1386, AbbSum251
	common /ddabbrev/ AbbSum341, AbbSum1744, AbbSum1632, AbbSum355
	common /ddabbrev/ AbbSum1945, AbbSum1864, AbbSum420, AbbSum236
	common /ddabbrev/ AbbSum789, AbbSum1376, AbbSum829, AbbSum878
	common /ddabbrev/ AbbSum1039, AbbSum595, AbbSum246, AbbSum285
	common /ddabbrev/ AbbSum1352, AbbSum800, AbbSum83, AbbSum1048
	common /ddabbrev/ AbbSum278, AbbSum399, AbbSum997, AbbSum1108
	common /ddabbrev/ AbbSum1860, AbbSum2102, AbbSum303, AbbSum893
	common /ddabbrev/ AbbSum1300, AbbSum2188, AbbSum310, AbbSum377
	common /ddabbrev/ AbbSum480, AbbSum1833, AbbSum316, AbbSum856
	common /ddabbrev/ AbbSum1309, AbbSum1621, AbbSum1843, AbbSum332
	common /ddabbrev/ AbbSum1141, AbbSum1666, AbbSum1429, AbbSum1204
	common /ddabbrev/ AbbSum945, AbbSum1610, AbbSum1123, AbbSum1494
	common /ddabbrev/ AbbSum1980, AbbSum615, AbbSum1981, AbbSum723
	common /ddabbrev/ AbbSum500, AbbSum1130, AbbSum1692, AbbSum1689
	common /ddabbrev/ AbbSum1501, AbbSum1075, AbbSum1532, AbbSum1780
	common /ddabbrev/ AbbSum1922, AbbSum652, AbbSum581, AbbSum948
	common /ddabbrev/ AbbSum1613, AbbSum1699, AbbSum1207, AbbSum1431
	common /ddabbrev/ AbbSum1684, AbbSum1866, AbbSum1846, AbbSum1932
	common /ddabbrev/ AbbSum1759, AbbSum582, AbbSum653, AbbSum654
	common /ddabbrev/ AbbSum583, AbbSum986, AbbSum1412, AbbSum526
	common /ddabbrev/ AbbSum701, AbbSum528, AbbSum655, AbbSum584
	common /ddabbrev/ AbbSum527, AbbSum470, AbbSum548, AbbSum699
	common /ddabbrev/ AbbSum525, AbbSum1456, AbbSum1166, AbbSum565
	common /ddabbrev/ AbbSum722, AbbSum1660, AbbSum673, AbbSum499
	common /ddabbrev/ AbbSum1715, AbbSum1688, AbbSum1691, AbbSum2008
	common /ddabbrev/ AbbSum1979, AbbSum1928, AbbSum721, AbbSum1765
	common /ddabbrev/ AbbSum498, AbbSum501, AbbSum874, AbbSum873
	common /ddabbrev/ AbbSum512, AbbSum250, AbbSum271, AbbSum248
	common /ddabbrev/ AbbSum269, AbbSum18, AbbSum1776, AbbSum2081
	common /ddabbrev/ AbbSum1045, AbbSum17, AbbSum1046, AbbSum2083
	common /ddabbrev/ AbbSum823, AbbSum1042, AbbSum1110, AbbSum2082
	common /ddabbrev/ AbbSum1043, AbbSum1047, AbbSum510, AbbSum729
	common /ddabbrev/ AbbSum514, AbbSum692, AbbSum2213, AbbSum2042
	common /ddabbrev/ AbbSum1321, AbbSum867, AbbSum524, AbbSum698
	common /ddabbrev/ AbbSum529, AbbSum535, AbbSum705, AbbSum1603
	common /ddabbrev/ AbbSum938, AbbSum1227, AbbSum1399, AbbSum1335
	common /ddabbrev/ AbbSum832, AbbSum1226, AbbSum1336, AbbSum833
	common /ddabbrev/ AbbSum834, AbbSum1337, AbbSum2185, AbbSum1290
	common /ddabbrev/ AbbSum881, AbbSum1114, AbbSum1486, AbbSum1551
	common /ddabbrev/ AbbSum1006, AbbSum1565, AbbSum1819, AbbSum1027
	common /ddabbrev/ AbbSum1891, AbbSum1455, AbbSum1165, AbbSum1025
	common /ddabbrev/ AbbSum909, AbbSum2118, AbbSum1277, AbbSum2293
	common /ddabbrev/ AbbSum2116, AbbSum1821, AbbSum1563, AbbSum1893
	common /ddabbrev/ AbbSum1023, AbbSum1026, AbbSum2292, AbbSum1564
	common /ddabbrev/ AbbSum1820, AbbSum2117, AbbSum1024, AbbSum1892
	common /ddabbrev/ AbbSum580, AbbSum651, AbbSum659, AbbSum2123
	common /ddabbrev/ AbbSum590, AbbSum2296, AbbSum616, AbbSum1870
	common /ddabbrev/ AbbSum1419, AbbSum1829, AbbSum1194, AbbSum622
	common /ddabbrev/ AbbSum644, AbbSum786, AbbSum1144, AbbSum1438
	common /ddabbrev/ AbbSum1919, AbbSum1777, AbbSum1547, AbbSum1000
	common /ddabbrev/ AbbSum397, AbbSum272, AbbSum1224, AbbSum1397
	common /ddabbrev/ AbbSum1924, AbbSum708, AbbSum1761, AbbSum475
	common /ddabbrev/ AbbSum883, AbbSum1292, AbbSum405, AbbSum252
	common /ddabbrev/ AbbSum1334, AbbSum831, AbbSum880, AbbSum557
	common /ddabbrev/ AbbSum665, AbbSum1289, AbbSum1333, AbbSum830
	common /ddabbrev/ AbbSum963, AbbSum882, AbbSum2084, AbbSum1291
	common /ddabbrev/ AbbSum2300, AbbSum1581, AbbSum1883, AbbSum386
	common /ddabbrev/ AbbSum290, AbbSum1806, AbbSum84, AbbSum156
	common /ddabbrev/ AbbSum682, AbbSum544, AbbSum1398, AbbSum1225
	common /ddabbrev/ AbbSum1049, AbbSum1514, AbbSum1658, AbbSum373
	common /ddabbrev/ AbbSum1484, AbbSum1714, AbbSum305, AbbSum1112
	common /ddabbrev/ AbbSum2301, AbbSum187, AbbSum884, AbbSum20
	common /ddabbrev/ AbbSum2085, AbbSum1293, AbbSum1053, AbbSum1518
	common /ddabbrev/ AbbSum1052, AbbSum1517, AbbSum1051, AbbSum1516
	common /ddabbrev/ AbbSum1050, AbbSum1515, AbbSum2023, AbbSum2232
	common /ddabbrev/ AbbSum795, AbbSum1380, AbbSum936, AbbSum1267
	common /ddabbrev/ AbbSum1349, AbbSum851, AbbSum885, AbbSum2014
	common /ddabbrev/ AbbSum1988, AbbSum1136, AbbSum1505, AbbSum899
	common /ddabbrev/ AbbSum1306, AbbSum927, AbbSum981, AbbSum2176
	common /ddabbrev/ AbbSum1601, AbbSum2248, AbbSum992, AbbSum1028
	common /ddabbrev/ AbbSum1054, AbbSum1796, AbbSum1073, AbbSum1088
	common /ddabbrev/ AbbSum1539, AbbSum1090, AbbSum1092, AbbSum1540
	common /ddabbrev/ AbbSum1089, AbbSum1713, AbbSum1674, AbbSum2276
	common /ddabbrev/ AbbSum2139, AbbSum1541, AbbSum1091, AbbSum1093
	common /ddabbrev/ AbbSum2212, AbbSum2052, AbbSum1771, AbbSum1470
	common /ddabbrev/ AbbSum1180, AbbSum1472, AbbSum1182, AbbSum2159
	common /ddabbrev/ AbbSum2158, AbbSum1184, AbbSum1183, AbbSum1474
	common /ddabbrev/ AbbSum1473, AbbSum1471, AbbSum1181, AbbSum1185
	common /ddabbrev/ AbbSum2073, AbbSum1247, AbbSum2198, AbbSum1415
	common /ddabbrev/ AbbSum1648, AbbSum568, AbbSum1656, AbbSum1731
	common /ddabbrev/ AbbSum1982, AbbSum1999, AbbSum2316, AbbSum2182
	common /ddabbrev/ AbbSum1994, AbbSum2224, AbbSum2037, AbbSum2119
	common /ddabbrev/ AbbSum206, AbbSum836, AbbSum2264, AbbSum2266
	common /ddabbrev/ AbbSum2155, AbbSum2157, AbbSum2265, AbbSum2156
	common /ddabbrev/ AbbSum2160, AbbSum88, AbbSum915, AbbSum2061
	common /ddabbrev/ AbbSum483, AbbSum1149, AbbSum973, AbbSum1675
	common /ddabbrev/ AbbSum1701, AbbSum926, AbbSum1259, AbbSum1238
	common /ddabbrev/ AbbSum2170, AbbSum1407, AbbSum2242, AbbSum1594
	common /ddabbrev/ AbbSum980, AbbSum465, AbbSum2094, AbbSum737
	common /ddabbrev/ AbbSum2186, AbbSum1724, AbbSum1647, AbbSum760
	common /ddabbrev/ AbbSum1384, AbbSum1274, AbbSum906, AbbSum1058
	common /ddabbrev/ AbbSum1522, AbbSum1897, AbbSum1782, AbbSum923
	common /ddabbrev/ AbbSum924, AbbSum1373, AbbSum785, AbbSum976
	common /ddabbrev/ AbbSum1592, AbbSum1721, AbbSum1643, AbbSum1258
	common /ddabbrev/ AbbSum925, AbbSum1608, AbbSum943, AbbSum1738
	common /ddabbrev/ AbbSum1625, AbbSum1966, AbbSum2018, AbbSum447
	common /ddabbrev/ AbbSum748, AbbSum336, AbbSum351, AbbSum2029
	common /ddabbrev/ AbbSum2216, AbbSum563, AbbSum671, AbbSum91
	common /ddabbrev/ AbbSum160, AbbSum1295, AbbSum837, AbbSum1697
	common /ddabbrev/ AbbSum601, AbbSum1959, AbbSum1201, AbbSum1152
	common /ddabbrev/ AbbSum1836, AbbSum320, AbbSum484, AbbSum1153
	common /ddabbrev/ AbbSum2193, AbbSum191, AbbSum1681, AbbSum628
	common /ddabbrev/ AbbSum1426, AbbSum1444, AbbSum1850, AbbSum363
	common /ddabbrev/ AbbSum711, AbbSum1876, AbbSum1445, AbbSum2066
	common /ddabbrev/ AbbSum28, AbbSum783, AbbSum782, AbbSum1372
	common /ddabbrev/ AbbSum1257, AbbSum922, AbbSum1371, AbbSum781
	common /ddabbrev/ AbbSum784, AbbSum1154, AbbSum1593, AbbSum977
	common /ddabbrev/ AbbSum1758, AbbSum647, AbbSum1931, AbbSum570
	common /ddabbrev/ AbbSum1939, AbbSum1753, AbbSum208, AbbSum428
	common /ddabbrev/ AbbSum1644, AbbSum230, AbbSum416, AbbSum979
	common /ddabbrev/ AbbSum978, AbbSum1646, AbbSum1723, AbbSum1443
	common /ddabbrev/ AbbSum1151, AbbSum1722, AbbSum1645, AbbSum2
	common /ddabbrev/ AbbSum571, AbbSum1661, AbbSum521, AbbSum572
	common /ddabbrev/ AbbSum97, AbbSum1811, AbbSum574, AbbSum308
	common /ddabbrev/ AbbSum573, AbbSum575, AbbSum516, AbbSum2115
	common /ddabbrev/ AbbSum1813, AbbSum1814, AbbSum1164, AbbSum235
	common /ddabbrev/ AbbSum1018, AbbSum275, AbbSum503, AbbSum57
	common /ddabbrev/ AbbSum517, AbbSum276, AbbSum520, AbbSum101
	common /ddabbrev/ AbbSum65, AbbSum31, AbbSum2068, AbbSum2131
	common /ddabbrev/ AbbSum2132, AbbSum2133, AbbSum279, AbbSum100
	common /ddabbrev/ AbbSum776, AbbSum778, AbbSum335, AbbSum229
	common /ddabbrev/ AbbSum1752, AbbSum1844, AbbSum779, AbbSum916
	common /ddabbrev/ AbbSum780, AbbSum914, AbbSum1635, AbbSum26
	common /ddabbrev/ AbbSum2075, AbbSum1640, AbbSum975, AbbSum913
	common /ddabbrev/ AbbSum27, AbbSum1237, AbbSum482, AbbSum1641
	common /ddabbrev/ AbbSum2092, AbbSum1148, AbbSum1639, AbbSum1134
	common /ddabbrev/ AbbSum1504, AbbSum1703, AbbSum1677, AbbSum1135
	common /ddabbrev/ AbbSum1690, AbbSum1693, AbbSum73, AbbSum170
	common /ddabbrev/ AbbSum1993, AbbSum2181, AbbSum2315, AbbSum1998
	common /ddabbrev/ AbbSum1824, AbbSum1896, AbbSum1930, AbbSum1769
	common /ddabbrev/ AbbSum1662, AbbSum1716, AbbSum533, AbbSum549
	common /ddabbrev/ AbbSum689, AbbSum534, AbbSum1673, AbbSum586
	common /ddabbrev/ AbbSum311, AbbSum378, AbbSum588, AbbSum1729
	common /ddabbrev/ AbbSum1654, AbbSum739, AbbSum471, AbbSum1600
	common /ddabbrev/ AbbSum989, AbbSum1781, AbbSum1923, AbbSum280
	common /ddabbrev/ AbbSum400, AbbSum1627, AbbSum3, AbbSum1118
	common /ddabbrev/ AbbSum794, AbbSum1278, AbbSum910, AbbSum1133
	common /ddabbrev/ AbbSum764, AbbSum935, AbbSum1131, AbbSum1502
	common /ddabbrev/ AbbSum2071, AbbSum1265, AbbSum933, AbbSum1801
	common /ddabbrev/ AbbSum1913, AbbSum2021, AbbSum2230, AbbSum792
	common /ddabbrev/ AbbSum1378, AbbSum1895, AbbSum1823, AbbSum2135
	common /ddabbrev/ AbbSum2274, AbbSum2137, AbbSum2136, AbbSum1086
	common /ddabbrev/ AbbSum1084, AbbSum1087, AbbSum1083, AbbSum1081
	common /ddabbrev/ AbbSum1082, AbbSum144, AbbSum104, AbbSum1132
	common /ddabbrev/ AbbSum1503, AbbSum793, AbbSum1379, AbbSum2022
	common /ddabbrev/ AbbSum2231, AbbSum1538, AbbSum1085, AbbSum1266
	common /ddabbrev/ AbbSum934, AbbSum658, AbbSum589, AbbSum704
	common /ddabbrev/ AbbSum532, AbbSum1672, AbbSum1712, AbbSum2048
	common /ddabbrev/ AbbSum2208, AbbSum1537, AbbSum1080, AbbSum1760
	common /ddabbrev/ AbbSum1933, AbbSum991, AbbSum990, AbbSum1730
	common /ddabbrev/ AbbSum1655, AbbSum146, AbbSum106, AbbSum587
	common /ddabbrev/ AbbSum657, AbbSum1599, AbbSum988, AbbSum145
	common /ddabbrev/ AbbSum105, AbbSum1243, AbbSum1413, AbbSum2245
	common /ddabbrev/ AbbSum2173, AbbSum987, AbbSum1598, AbbSum508
	common /ddabbrev/ AbbSum727, AbbSum1986, AbbSum2012, AbbSum506
	common /ddabbrev/ AbbSum725, AbbSum2010, AbbSum1984, AbbSum2049
	common /ddabbrev/ AbbSum2209, AbbSum216, AbbSum433, AbbSum1757
	common /ddabbrev/ AbbSum1943, AbbSum898, AbbSum1305, AbbSum643
	common /ddabbrev/ AbbSum621, AbbSum2261, AbbSum1177, AbbSum2152
	common /ddabbrev/ AbbSum1467, AbbSum2262, AbbSum2153, AbbSum1178
	common /ddabbrev/ AbbSum1468, AbbSum1847, AbbSum1867, AbbSum1245
	common /ddabbrev/ AbbSum1246, AbbSum2175, AbbSum2247, AbbSum1414
	common /ddabbrev/ AbbSum2174, AbbSum1244, AbbSum2246, AbbSum1948
	common /ddabbrev/ AbbSum1862, AbbSum1839, AbbSum1879, AbbSum367
	common /ddabbrev/ AbbSum325, AbbSum2210, AbbSum2050, AbbSum1348
	common /ddabbrev/ AbbSum850, AbbSum2036, AbbSum2223, AbbSum2011
	common /ddabbrev/ AbbSum726, AbbSum1985, AbbSum507, AbbSum1911
	common /ddabbrev/ AbbSum1799, AbbSum1651, AbbSum1727, AbbSum1076
	common /ddabbrev/ AbbSum1533, AbbSum2134, AbbSum2273, AbbSum1534
	common /ddabbrev/ AbbSum1077, AbbSum2020, AbbSum2229, AbbSum790
	common /ddabbrev/ AbbSum1377, AbbSum2069, AbbSum2195, AbbSum931
	common /ddabbrev/ AbbSum1263, AbbSum1800, AbbSum1912, AbbSum1535
	common /ddabbrev/ AbbSum1078, AbbSum1894, AbbSum1822, AbbSum142
	common /ddabbrev/ AbbSum102, AbbSum2122, AbbSum1568, AbbSum1032
	common /ddabbrev/ AbbSum1802, AbbSum1914, AbbSum181, AbbSum62
	common /ddabbrev/ AbbSum868, AbbSum123, AbbSum139, AbbSum1322
	common /ddabbrev/ AbbSum1364, AbbSum813, AbbSum390, AbbSum296
	common /ddabbrev/ AbbSum1854, AbbSum1954, AbbSum409, AbbSum259
	common /ddabbrev/ AbbSum1853, AbbSum1956, AbbSum1687, AbbSum1694
	common /ddabbrev/ AbbSum1567, AbbSum1030, AbbSum2295, AbbSum2121
	common /ddabbrev/ AbbSum2275, AbbSum2138, AbbSum1079, AbbSum1536
	common /ddabbrev/ AbbSum703, AbbSum531, AbbSum702, AbbSum530
	common /ddabbrev/ AbbSum1710, AbbSum1670, AbbSum1711, AbbSum1671
	common /ddabbrev/ AbbSum2009, AbbSum724, AbbSum505, AbbSum1983
	common /ddabbrev/ AbbSum585, AbbSum656, AbbSum2096, AbbSum2310
	common /ddabbrev/ AbbSum1987, AbbSum509, AbbSum728, AbbSum2013
	common /ddabbrev/ AbbSum143, AbbSum103, AbbSum1652, AbbSum1168
	common /ddabbrev/ AbbSum2143, AbbSum2252, AbbSum1458, AbbSum1728
	common /ddabbrev/ AbbSum1653, AbbSum2072, AbbSum2197, AbbSum1031
	common /ddabbrev/ AbbSum949, AbbSum1614, AbbSum1700, AbbSum1686
	common /ddabbrev/ AbbSum1960, AbbSum1851, AbbSum1432, AbbSum1208
	common /ddabbrev/ AbbSum1967, AbbSum2019, AbbSum751, AbbSum451
	common /ddabbrev/ AbbSum675, AbbSum566, AbbSum94, AbbSum161
	common /ddabbrev/ AbbSum32, AbbSum192, AbbSum641, AbbSum617
	common /ddabbrev/ AbbSum791, AbbSum238, AbbSum421, AbbSum894
	common /ddabbrev/ AbbSum1301, AbbSum2294, AbbSum2120, AbbSum1029
	common /ddabbrev/ AbbSum1566, AbbSum2051, AbbSum2211, AbbSum2263
	common /ddabbrev/ AbbSum1179, AbbSum1469, AbbSum2154, AbbSum2070
	common /ddabbrev/ AbbSum2196, AbbSum932, AbbSum1264, AbbSum1745
	common /ddabbrev/ AbbSum1634, AbbSum1946, AbbSum1865, AbbSum343
	common /ddabbrev/ AbbSum356, AbbSum1344, AbbSum846, AbbSum2219
	common /ddabbrev/ AbbSum2032, AbbSum2144, AbbSum1169, AbbSum2253
	common /ddabbrev/ AbbSum1459, AbbSum68, AbbSum166, AbbSum1989
	common /ddabbrev/ AbbSum2177, AbbSum2311, AbbSum1995, AbbSum70
	common /ddabbrev/ AbbSum1991, AbbSum2179, AbbSum2313, AbbSum1996
	common /ddabbrev/ AbbSum1990, AbbSum69, AbbSum2178, AbbSum2312
	common /ddabbrev/ AbbSum167, AbbSum1170, AbbSum2145, AbbSum2254
	common /ddabbrev/ AbbSum1460, AbbSum1174, AbbSum2149, AbbSum2258
	common /ddabbrev/ AbbSum1464, AbbSum642, AbbSum618, AbbSum1172
	common /ddabbrev/ AbbSum2147, AbbSum2256, AbbSum1462, AbbSum619
	common /ddabbrev/ AbbSum620, AbbSum1171, AbbSum2146, AbbSum2255
	common /ddabbrev/ AbbSum1461, AbbSum2207, AbbSum2047, AbbSum896
	common /ddabbrev/ AbbSum1303, AbbSum895, AbbSum1302, AbbSum1175
	common /ddabbrev/ AbbSum1465, AbbSum2259, AbbSum2150, AbbSum897
	common /ddabbrev/ AbbSum1304, AbbSum2257, AbbSum1173, AbbSum2148
	common /ddabbrev/ AbbSum1463, AbbSum1346, AbbSum848, AbbSum2221
	common /ddabbrev/ AbbSum2034, AbbSum1347, AbbSum849, AbbSum2222
	common /ddabbrev/ AbbSum2035, AbbSum847, AbbSum1345, AbbSum2033
	common /ddabbrev/ AbbSum2220, AbbSum71, AbbSum168, AbbSum1992
	common /ddabbrev/ AbbSum2180, AbbSum2314, AbbSum1997, AbbSum2260
	common /ddabbrev/ AbbSum1176, AbbSum1466, AbbSum2151, AbbSum592
	common /ddabbrev/ AbbSum593, AbbSum1138, AbbSum282, AbbSum245
	common /ddabbrev/ AbbSum442, AbbSum2101, AbbSum301, AbbSum298
	common /ddabbrev/ AbbSum994, AbbSum476, AbbSum912, AbbSum1230
	common /ddabbrev/ AbbSum477, AbbSum197, AbbSum313, AbbSum200
	common /ddabbrev/ AbbSum479, AbbSum1003, AbbSum903, AbbSum2074
	common /ddabbrev/ AbbSum1842, AbbSum328, AbbSum330, AbbSum610
	common /ddabbrev/ AbbSum2077, AbbSum213, AbbSum1976, AbbSum2041
	common /ddabbrev/ AbbSum1845, AbbSum45, AbbSum2039, AbbSum2076
	common /ddabbrev/ AbbSum1659, AbbSum96, AbbSum108, AbbSum46
	common /ddabbrev/ AbbSum2040, AbbSum49, AbbSum1011, AbbSum112
	common /ddabbrev/ AbbSum53, AbbSum1016, AbbSum841, AbbSum608
	common /ddabbrev/ AbbSum606, AbbSum607, AbbSum1059, AbbSum2045
	common /ddabbrev/ AbbSum1783, AbbSum1061, AbbSum1062, AbbSum1012
	common /ddabbrev/ AbbSum1063, AbbSum116, AbbSum494, AbbSum2129
	common /ddabbrev/ AbbSum1789, AbbSum1974, AbbSum1975, AbbSum983
	common /ddabbrev/ AbbSum255, AbbSum488, AbbSum52, AbbSum256
	common /ddabbrev/ AbbSum293, AbbSum292, AbbSum1767, AbbSum1906
	common /ddabbrev/ AbbSum1793, AbbSum1262, AbbSum1529, AbbSum1069
	common /ddabbrev/ AbbSum930, AbbSum2194, AbbSum2067, AbbSum117
	common /ddabbrev/ AbbSum865, AbbSum1319, AbbSum134, AbbSum1070
	common /ddabbrev/ AbbSum1794, AbbSum1907, AbbSum1122, AbbSum1493
	common /ddabbrev/ AbbSum808, AbbSum1359, AbbSum389, AbbSum295
	common /ddabbrev/ AbbSum258, AbbSum407, AbbSum366, AbbSum324
	common /ddabbrev/ AbbSum984, AbbSum1596, AbbSum2007, AbbSum1942
	common /ddabbrev/ AbbSum2006, AbbSum1978, AbbSum1756, AbbSum1977
	common /ddabbrev/ AbbSum1709, AbbSum1448, AbbSum1669, AbbSum1157
	common /ddabbrev/ AbbSum1298, AbbSum889, AbbSum1340, AbbSum2218
	common /ddabbrev/ AbbSum842, AbbSum2031, AbbSum55, AbbSum2027
	common /ddabbrev/ AbbSum1947, AbbSum1878, AbbSum496, AbbSum1017
	common /ddabbrev/ AbbSum1852, AbbSum497, AbbSum2206, AbbSum176
	common /ddabbrev/ AbbSum2228, AbbSum1861, AbbSum1838, AbbSum719
	common /ddabbrev/ AbbSum720, AbbSum1955, AbbSum1559, AbbSum2046
	common /ddabbrev/ AbbSum1905, AbbSum1792, AbbSum1528, AbbSum1068
	common /ddabbrev/ AbbSum1702, AbbSum1676, AbbSum750, AbbSum449
	common /ddabbrev/ AbbSum614, AbbSum214, AbbSum432, AbbSum1241
	common /ddabbrev/ AbbSum1410, AbbSum1072, AbbSum2172, AbbSum2244
	common /ddabbrev/ AbbSum1530, AbbSum1071, AbbSum638, AbbSum613
	common /ddabbrev/ AbbSum1795, AbbSum611, AbbSum1908, AbbSum636
	common /ddabbrev/ AbbSum2272, AbbSum2130, AbbSum612, AbbSum1791
	common /ddabbrev/ AbbSum637, AbbSum1904, AbbSum396, AbbSum681
	common /ddabbrev/ AbbSum1602, AbbSum81, AbbSum82, AbbSum1193
	common /ddabbrev/ AbbSum1580, AbbSum597, AbbSum1579, AbbSum372
	common /ddabbrev/ AbbSum1918, AbbSum80, AbbSum1577, AbbSum1576
	common /ddabbrev/ AbbSum153, AbbSum959, AbbSum958, AbbSum304
	common /ddabbrev/ AbbSum1775, AbbSum543, AbbSum937, AbbSum154
	common /ddabbrev/ AbbSum155, AbbSum1418, AbbSum962, AbbSum445
	common /ddabbrev/ AbbSum961, AbbSum270, AbbSum37, AbbSum1221
	common /ddabbrev/ AbbSum1307, AbbSum1475, AbbSum1324, AbbSum819
	common /ddabbrev/ AbbSum124, AbbSum1394, AbbSum854, AbbSum1101
	common /ddabbrev/ AbbSum221, AbbSum1420, AbbSum1584, AbbSum1582
	common /ddabbrev/ AbbSum1585, AbbSum967, AbbSum964, AbbSum966
	common /ddabbrev/ AbbSum326, AbbSum1195, AbbSum1520, AbbSum683
	common /ddabbrev/ AbbSum757, AbbSum21, AbbSum902, AbbSum2163
	common /ddabbrev/ AbbSum1937, AbbSum558, AbbSum1251, AbbSum911
	common /ddabbrev/ AbbSum1056, AbbSum461, AbbSum1381, AbbSum85
	common /ddabbrev/ AbbSum1271, AbbSum2235, AbbSum1749, AbbSum666
	common /ddabbrev/ AbbSum22, AbbSum1732, AbbSum758, AbbSum1382
	common /ddabbrev/ AbbSum86, AbbSum1617, AbbSum1485, AbbSum189
	common /ddabbrev/ AbbSum2282, AbbSum464, AbbSum560, AbbSum1005
	common /ddabbrev/ AbbSum600, AbbSum1405, AbbSum188, AbbSum190
	common /ddabbrev/ AbbSum2192, AbbSum599, AbbSum1294, AbbSum2065
	common /ddabbrev/ AbbSum626, AbbSum835, AbbSum1113, AbbSum24
	common /ddabbrev/ AbbSum2105, AbbSum736, AbbSum668, AbbSum1550
	common /ddabbrev/ AbbSum627, AbbSum1235, AbbSum23, AbbSum25
	common /ddabbrev/ AbbSum411, AbbSum141, AbbSum1421, AbbSum1004
	common /ddabbrev/ AbbSum273, AbbSum970, AbbSum223, AbbSum99
	common /ddabbrev/ AbbSum1196, AbbSum1588, AbbSum1549, AbbSum374
	common /ddabbrev/ AbbSum344, AbbSum816, AbbSum76, AbbSum539
	common /ddabbrev/ AbbSum1189, AbbSum1416, AbbSum218, AbbSum1280
	common /ddabbrev/ AbbSum152, AbbSum676, AbbSum1527, AbbSum928
	common /ddabbrev/ AbbSum133, AbbSum1428, AbbSum864, AbbSum944
	common /ddabbrev/ AbbSum1628, AbbSum1297, AbbSum1060, AbbSum1296
	common /ddabbrev/ AbbSum131, AbbSum863, AbbSum2288, AbbSum1595
	common /ddabbrev/ AbbSum132, AbbSum862, AbbSum2287, AbbSum114
	common /ddabbrev/ AbbSum1316, AbbSum2110, AbbSum982, AbbSum1260
	common /ddabbrev/ AbbSum115, AbbSum1203, AbbSum1318, AbbSum1609
	common /ddabbrev/ AbbSum1740, AbbSum888, AbbSum1523, AbbSum887
	common /ddabbrev/ AbbSum113, AbbSum1317, AbbSum2111, AbbSum1067
	common /ddabbrev/ AbbSum1558, AbbSum406, AbbSum294, AbbSum839
	common /ddabbrev/ AbbSum840, AbbSum175, AbbSum173, AbbSum1557
	common /ddabbrev/ AbbSum2109, AbbSum1556, AbbSum174, AbbSum1358
	common /ddabbrev/ AbbSum2268, AbbSum1013, AbbSum51, AbbSum807
	common /ddabbrev/ AbbSum2125, AbbSum2286, AbbSum1014, AbbSum1338
	common /ddabbrev/ AbbSum1339, AbbSum54, AbbSum50, AbbSum1015
	common /ddabbrev/ AbbSum257, AbbSum388, AbbSum714, AbbSum1971
	common /ddabbrev/ AbbSum1972, AbbSum1926, AbbSum1927, AbbSum2030
	common /ddabbrev/ AbbSum2205, AbbSum1925, AbbSum1492, AbbSum1385
	common /ddabbrev/ AbbSum605, AbbSum1409, AbbSum1969, AbbSum1968
	common /ddabbrev/ AbbSum1970, AbbSum467, AbbSum604, AbbSum431
	common /ddabbrev/ AbbSum430, AbbSum2203, AbbSum1275, AbbSum929
	common /ddabbrev/ AbbSum1766, AbbSum564, AbbSum29, AbbSum718
	common /ddabbrev/ AbbSum749, AbbSum716, AbbSum717, AbbSum1755
	common /ddabbrev/ AbbSum1155, AbbSum365, AbbSum1707, AbbSum1156
	common /ddabbrev/ AbbSum1708, AbbSum713, AbbSum2187, AbbSum1787
	common /ddabbrev/ AbbSum603, AbbSum1526, AbbSum1525, AbbSum2269
	common /ddabbrev/ AbbSum712, AbbSum1786, AbbSum602, AbbSum1408
	common /ddabbrev/ AbbSum429, AbbSum1120, AbbSum2270, AbbSum1754
	common /ddabbrev/ AbbSum1877, AbbSum364, AbbSum2171, AbbSum1784
	common /ddabbrev/ AbbSum1785, AbbSum1788, AbbSum1790, AbbSum1973
	common /ddabbrev/ AbbSum609, AbbSum715, AbbSum1524, AbbSum2271
	common /ddabbrev/ AbbSum2128, AbbSum1064, AbbSum2005, AbbSum635
	common /ddabbrev/ AbbSum491, AbbSum1239, AbbSum210, AbbSum1903
	common /ddabbrev/ AbbSum1902, AbbSum1491, AbbSum2127, AbbSum1940
	common /ddabbrev/ AbbSum1837, AbbSum322, AbbSum2243, AbbSum1898
	common /ddabbrev/ AbbSum1899, AbbSum631, AbbSum1066, AbbSum487
	common /ddabbrev/ AbbSum1900, AbbSum1065, AbbSum2126, AbbSum632
	common /ddabbrev/ AbbSum1901, AbbSum1121, AbbSum761, AbbSum634
	common /ddabbrev/ AbbSum1240, AbbSum2001, AbbSum2000, AbbSum2002
	common /ddabbrev/ AbbSum687, AbbSum633, AbbSum212, AbbSum211
	common /ddabbrev/ AbbSum2059, AbbSum907, AbbSum1261, AbbSum672
	common /ddabbrev/ AbbSum92, AbbSum495, AbbSum448, AbbSum492
	common /ddabbrev/ AbbSum493, AbbSum1941, AbbSum1446, AbbSum323
	common /ddabbrev/ AbbSum1667, AbbSum1447, AbbSum1668, AbbSum489
	common /ddabbrev/ AbbSum1764, AbbSum2003, AbbSum2004, AbbSum2217
	common /ddabbrev/ AbbSum2044, AbbSum1762, AbbSum1763, AbbSum490
	common /ddabbrev/ AbbSum518, AbbSum890, AbbSum700, AbbSum1021
	common /ddabbrev/ AbbSum375, AbbSum376, AbbSum30, AbbSum674
	common /ddabbrev/ AbbSum1560, AbbSum1778, AbbSum1612, AbbSum1779
	common /ddabbrev/ AbbSum519, AbbSum1562, AbbSum1611, AbbSum1863
	common /ddabbrev/ AbbSum946, AbbSum1944, AbbSum1022, AbbSum469
	common /ddabbrev/ AbbSum1561, AbbSum307, AbbSum309, AbbSum93
	common /ddabbrev/ AbbSum450, AbbSum1019, AbbSum1920, AbbSum947
	common /ddabbrev/ AbbSum1921, AbbSum695, AbbSum694, AbbSum1299
	common /ddabbrev/ AbbSum1282, AbbSum286, AbbSum1103, AbbSum299
	common /ddabbrev/ AbbSum1859, AbbSum870, AbbSum247, AbbSum1478
	common /ddabbrev/ AbbSum369, AbbSum1953, AbbSum393, AbbSum63
	common /ddabbrev/ AbbSum1773, AbbSum125, AbbSum1858, AbbSum1281
	common /ddabbrev/ AbbSum1951, AbbSum1950, AbbSum511, AbbSum1949
	common /ddabbrev/ AbbSum732, AbbSum1350, AbbSum798, AbbSum1855
	common /ddabbrev/ AbbSum459, AbbSum266, AbbSum162, AbbSum1916
	common /ddabbrev/ AbbSum690, AbbSum1856, AbbSum1857, AbbSum38
	common /ddabbrev/ AbbSum1952, AbbSum869, AbbSum1519, AbbSum127
	common /ddabbrev/ AbbSum1586, AbbSum1002, AbbSum478, AbbSum64
	common /ddabbrev/ AbbSum2164, AbbSum387, AbbSum1873, AbbSum361
	common /ddabbrev/ AbbSum201, AbbSum2304, AbbSum1936, AbbSum965
	common /ddabbrev/ AbbSum2162, AbbSum1934, AbbSum199, AbbSum1402
	common /ddabbrev/ AbbSum1228, AbbSum198, AbbSum1935, AbbSum1831
	common /ddabbrev/ AbbSum312, AbbSum1830, AbbSum314, AbbSum1229
	common /ddabbrev/ AbbSum2302, AbbSum2303, AbbSum2161, AbbSum2233
	common /ddabbrev/ AbbSum2086, AbbSum2087, AbbSum1747, AbbSum1401
	common /ddabbrev/ AbbSum1872, AbbSum359, AbbSum1871, AbbSum360
	common /ddabbrev/ AbbSum1583, AbbSum1746, AbbSum423, AbbSum1231
	common /ddabbrev/ AbbSum1400, AbbSum422, AbbSum2234, AbbSum1748
	common /ddabbrev/ AbbSum1055, AbbSum41, AbbSum968, AbbSum1548
	common /ddabbrev/ AbbSum163, AbbSum709, AbbSum2088, AbbSum2236
	common /ddabbrev/ AbbSum253, AbbSum424, AbbSum1832, AbbSum315
	common /ddabbrev/ AbbSum1618, AbbSum203, AbbSum225, AbbSum1885
	common /ddabbrev/ AbbSum2089, AbbSum1619, AbbSum645, AbbSum128
	common /ddabbrev/ AbbSum684, AbbSum222, AbbSum1807, AbbSum969
	common /ddabbrev/ AbbSum347, AbbSum462, AbbSum1423, AbbSum224
	common /ddabbrev/ AbbSum1422, AbbSum1272, AbbSum1365, AbbSum327
	common /ddabbrev/ AbbSum1605, AbbSum331, AbbSum226, AbbSum1620
	common /ddabbrev/ AbbSum1604, AbbSum1735, AbbSum939, AbbSum1197
	common /ddabbrev/ AbbSum904, AbbSum771, AbbSum414, AbbSum346
	common /ddabbrev/ AbbSum940, AbbSum348, AbbSum1198, AbbSum412
	common /ddabbrev/ AbbSum317, AbbSum413, AbbSum1808, AbbSum2305
	common /ddabbrev/ AbbSum1733, AbbSum1809, AbbSum42, AbbSum545
	common /ddabbrev/ AbbSum410, AbbSum1884, AbbSum734, AbbSum329
	common /ddabbrev/ AbbSum1587, AbbSum567, AbbSum1734, AbbSum920
	common /ddabbrev/ AbbSum773, AbbSum772, AbbSum775, AbbSum129
	common /ddabbrev/ AbbSum921, AbbSum801, AbbSum569, AbbSum1150
	common /ddabbrev/ AbbSum1834, AbbSum1835, AbbSum1200, AbbSum1607
	common /ddabbrev/ AbbSum334, AbbSum774, AbbSum227, AbbSum349
	common /ddabbrev/ AbbSum1848, AbbSum941, AbbSum362, AbbSum1886
	common /ddabbrev/ AbbSum2090, AbbSum2091, AbbSum710, AbbSum1146
	common /ddabbrev/ AbbSum1145, AbbSum1636, AbbSum1370, AbbSum1199
	common /ddabbrev/ AbbSum228, AbbSum1849, AbbSum2093, AbbSum2308
	common /ddabbrev/ AbbSum1424, AbbSum415, AbbSum1958, AbbSum1717
	common /ddabbrev/ AbbSum777, AbbSum481, AbbSum1440, AbbSum1439
	common /ddabbrev/ AbbSum1369, AbbSum43, AbbSum1256, AbbSum1310
	common /ddabbrev/ AbbSum646, AbbSum1442, AbbSum1874, AbbSum1875
	common /ddabbrev/ AbbSum1425, AbbSum942, AbbSum350, AbbSum1368
	common /ddabbrev/ AbbSum318, AbbSum333, AbbSum1957, AbbSum1606
	common /ddabbrev/ AbbSum319, AbbSum2306, AbbSum2307, AbbSum1810
	common /ddabbrev/ AbbSum1255, AbbSum1367, AbbSum1366, AbbSum1314
	common /ddabbrev/ AbbSum838, AbbSum804, AbbSum1312, AbbSum806
	common /ddabbrev/ AbbSum861, AbbSum803, AbbSum1315, AbbSum1354
	common /ddabbrev/ AbbSum1357, AbbSum886, AbbSum1355, AbbSum858
	common /ddabbrev/ AbbSum860, AbbSum149, AbbSum1571, AbbSum2053
	common /ddabbrev/ AbbSum1268, AbbSum150, AbbSum540, AbbSum662
	common /ddabbrev/ AbbSum1868, AbbSum357, AbbSum151, AbbSum1574
	common /ddabbrev/ AbbSum148, AbbSum7, AbbSum955, AbbSum2199
	common /ddabbrev/ AbbSum752, AbbSum8, AbbSum677, AbbSum554
	common /ddabbrev/ AbbSum1827, AbbSum194, AbbSum75, AbbSum9
	common /ddabbrev/ AbbSum953, AbbSum1487, AbbSum209, AbbSum630
	common /ddabbrev/ AbbSum337, AbbSum485, AbbSum1739, AbbSum1202
	common /ddabbrev/ AbbSum1554, AbbSum1116, AbbSum1488, AbbSum1009
	common /ddabbrev/ AbbSum321, AbbSum486, AbbSum231, AbbSum629
	common /ddabbrev/ AbbSum1626, AbbSum1427, AbbSum1115, AbbSum291
	common /ddabbrev/ AbbSum1490, AbbSum1311, AbbSum1356, AbbSum274
	common /ddabbrev/ AbbSum2202, AbbSum802, AbbSum1313, AbbSum686
	common /ddabbrev/ AbbSum2285, AbbSum1489, AbbSum2024, AbbSum130
	common /ddabbrev/ AbbSum44, AbbSum1008, AbbSum2025, AbbSum47
	common /ddabbrev/ AbbSum1010, AbbSum1007, AbbSum2284, AbbSum48
	common /ddabbrev/ AbbSum2026, AbbSum2283, AbbSum2106, AbbSum111
	common /ddabbrev/ AbbSum2227, AbbSum2107, AbbSum1552, AbbSum2226
	common /ddabbrev/ AbbSum110, AbbSum1555, AbbSum1553, AbbSum254
	common /ddabbrev/ AbbSum1119, AbbSum857, AbbSum805, AbbSum306
	common /ddabbrev/ AbbSum2058, AbbSum1353, AbbSum859, AbbSum1117
	common /ddabbrev/ AbbSum2225, AbbSum109, AbbSum172, AbbSum2108
	common /ddabbrev/ AbbSum466, AbbSum1818, AbbSum1375, AbbSum179
	common /ddabbrev/ AbbSum1362, AbbSum1450, AbbSum1452, AbbSum178
	common /ddabbrev/ AbbSum1361, AbbSum1454, AbbSum639, AbbSum1242
	common /ddabbrev/ AbbSum1449, AbbSum1451, AbbSum1816, AbbSum985
	common /ddabbrev/ AbbSum1374, AbbSum2095, AbbSum1682, AbbSum1430
	common /ddabbrev/ AbbSum1817, AbbSum1815, AbbSum576, AbbSum578
	common /ddabbrev/ AbbSum1360, AbbSum177, AbbSum577, AbbSum649
	common /ddabbrev/ AbbSum648, AbbSum650, AbbSum809, AbbSum56
	common /ddabbrev/ AbbSum1887, AbbSum1889, AbbSum1698, AbbSum1205
	common /ddabbrev/ AbbSum788, AbbSum59, AbbSum811, AbbSum1159
	common /ddabbrev/ AbbSum1161, AbbSum58, AbbSum810, AbbSum1163
	common /ddabbrev/ AbbSum502, AbbSum1411, AbbSum1158, AbbSum1160
	common /ddabbrev/ AbbSum1597, AbbSum787, AbbSum2309, AbbSum1888
	common /ddabbrev/ AbbSum1890, AbbSum2249, AbbSum1633, AbbSum215
	common /ddabbrev/ AbbSum408, AbbSum2250, AbbSum72, AbbSum1768
	common /ddabbrev/ AbbSum1650, AbbSum1798, AbbSum164, AbbSum67
	common /ddabbrev/ AbbSum2251, AbbSum1797, AbbSum1649, AbbSum1909
	common /ddabbrev/ AbbSum1725, AbbSum169, AbbSum1929, AbbSum1726
	common /ddabbrev/ AbbSum1910, AbbSum66, AbbSum165, AbbSum2140
	common /ddabbrev/ AbbSum2141, AbbSum1685, AbbSum342, AbbSum237
	common /ddabbrev/ AbbSum2142, AbbSum561, AbbSum1938, AbbSum2267
	common /ddabbrev/ AbbSum427, AbbSum562, AbbSum905, AbbSum1622
	common /ddabbrev/ AbbSum1591, AbbSum204, AbbSum1680, AbbSum205
	common /ddabbrev/ AbbSum1678, AbbSum2201, AbbSum1737, AbbSum971
	common /ddabbrev/ AbbSum972, AbbSum1720, AbbSum463, AbbSum685
	common /ddabbrev/ AbbSum1057, AbbSum559, AbbSum446, AbbSum87
	common /ddabbrev/ AbbSum90, AbbSum2241, AbbSum89, AbbSum1719
	common /ddabbrev/ AbbSum2237, AbbSum1383, AbbSum1254, AbbSum1147
	common /ddabbrev/ AbbSum2239, AbbSum1718, AbbSum1253, AbbSum2189
	common /ddabbrev/ AbbSum2191, AbbSum1679, AbbSum2215, AbbSum1236
	common /ddabbrev/ AbbSum1234, AbbSum2240, AbbSum2238, AbbSum1233
	common /ddabbrev/ AbbSum1252, AbbSum2190, AbbSum2063, AbbSum1403
	common /ddabbrev/ AbbSum917, AbbSum1404, AbbSum1406, AbbSum2168
	common /ddabbrev/ AbbSum2166, AbbSum2064, AbbSum1696, AbbSum2028
	common /ddabbrev/ AbbSum759, AbbSum919, AbbSum1637, AbbSum918
	common /ddabbrev/ AbbSum2062, AbbSum2167, AbbSum1441, AbbSum2165
	common /ddabbrev/ AbbSum207, AbbSum670, AbbSum1273, AbbSum1736
	common /ddabbrev/ AbbSum974, AbbSum425, AbbSum1624, AbbSum426
	common /ddabbrev/ AbbSum1695, AbbSum2057, AbbSum1623, AbbSum1589
	common /ddabbrev/ AbbSum1590, AbbSum1642, AbbSum735, AbbSum546
	common /ddabbrev/ AbbSum1521, AbbSum1638, AbbSum667, AbbSum747
	common /ddabbrev/ AbbSum157, AbbSum159, AbbSum2169, AbbSum158
	common /ddabbrev/ AbbSum2124, AbbSum1751, AbbSum669, AbbSum1125
	common /ddabbrev/ AbbSum398, AbbSum137, AbbSum844, AbbSum136
	common /ddabbrev/ AbbSum1531, AbbSum1126, AbbSum843, AbbSum1127
	common /ddabbrev/ AbbSum763, AbbSum468, AbbSum2289, AbbSum418
	common /ddabbrev/ AbbSum2290, AbbSum688, AbbSum908, AbbSum419
	common /ddabbrev/ AbbSum353, AbbSum2060, AbbSum1630, AbbSum2114
	common /ddabbrev/ AbbSum354, AbbSum1631, AbbSum1124, AbbSum866
	common /ddabbrev/ AbbSum135, AbbSum696, AbbSum352, AbbSum1629
	common /ddabbrev/ AbbSum417, AbbSum1162, AbbSum1499, AbbSum697
	common /ddabbrev/ AbbSum693, AbbSum515, AbbSum1128, AbbSum523
	common /ddabbrev/ AbbSum1320, AbbSum118, AbbSum522, AbbSum1453
	common /ddabbrev/ AbbSum232, AbbSum338, AbbSum1741, AbbSum1495
	common /ddabbrev/ AbbSum121, AbbSum1342, AbbSum119, AbbSum1074
	common /ddabbrev/ AbbSum1497, AbbSum1341, AbbSum1498, AbbSum1387
	common /ddabbrev/ AbbSum340, AbbSum1743, AbbSum2112, AbbSum233
	common /ddabbrev/ AbbSum2113, AbbSum547, AbbSum1276, AbbSum234
	common /ddabbrev/ AbbSum339, AbbSum2204, AbbSum1742, AbbSum2291
	common /ddabbrev/ AbbSum738, AbbSum277, AbbSum1496, AbbSum552
	common /ddabbrev/ AbbSum731, AbbSum368, AbbSum182, AbbSum553
	common /ddabbrev/ AbbSum591, AbbSum742, AbbSum1573, AbbSum379
	common /ddabbrev/ AbbSum743, AbbSum11, AbbSum730, AbbSum952
	common /ddabbrev/ AbbSum2017, AbbSum1250, AbbSum740, AbbSum663
	common /ddabbrev/ AbbSum741, AbbSum1390, AbbSum2016, AbbSum1249
	common /ddabbrev/ AbbSum1389, AbbSum1279, AbbSum2015, AbbSum1248
	common /ddabbrev/ AbbSum1570, AbbSum1391, AbbSum1388, AbbSum1214
	common /ddabbrev/ AbbSum1217, AbbSum951, AbbSum814, AbbSum1961
	common /ddabbrev/ AbbSum765, AbbSum1215, AbbSum1216, AbbSum555
	common /ddabbrev/ AbbSum437, AbbSum1962, AbbSum767, AbbSum436
	common /ddabbrev/ AbbSum456, AbbSum261, AbbSum10, AbbSum661
	common /ddabbrev/ AbbSum472, AbbSum438, AbbSum240, AbbSum439
	common /ddabbrev/ AbbSum183, AbbSum454, AbbSum1572, AbbSum1964
	common /ddabbrev/ AbbSum768, AbbSum954, AbbSum660, AbbSum297
	common /ddabbrev/ AbbSum817, AbbSum1506, AbbSum402, AbbSum140
	common /ddabbrev/ AbbSum1476, AbbSum594, AbbSum284, AbbSum382
	common /ddabbrev/ AbbSum263, AbbSum1323, AbbSum1034, AbbSum243
	common /ddabbrev/ AbbSum623, AbbSum98, AbbSum1102, AbbSum1393
	common /ddabbrev/ AbbSum458, AbbSum1507, AbbSum381, AbbSum358
	common /ddabbrev/ AbbSum403, AbbSum1508, AbbSum13, AbbSum1037
	common /ddabbrev/ AbbSum744, AbbSum380, AbbSum1664, AbbSum1140
	common /ddabbrev/ AbbSum401, AbbSum1510, AbbSum1704, AbbSum1139
	common /ddabbrev/ AbbSum1392, AbbSum1706, AbbSum1137, AbbSum1219
	common /ddabbrev/ AbbSum1433, AbbSum1665, AbbSum242, AbbSum1434
	common /ddabbrev/ AbbSum1663, AbbSum1038, AbbSum1220, AbbSum678
	common /ddabbrev/ AbbSum1035, AbbSum283, AbbSum195, AbbSum78
	common /ddabbrev/ AbbSum1509, AbbSum441, AbbSum1435, AbbSum281
	common /ddabbrev/ AbbSum1705, AbbSum1036, AbbSum244, AbbSum541
	common /ddabbrev/ AbbSum1774, AbbSum268, AbbSum14, AbbSum1105
	common /ddabbrev/ AbbSum265, AbbSum1104, AbbSum1269, AbbSum1480
	common /ddabbrev/ AbbSum1192, AbbSum219, AbbSum267, AbbSum300
	common /ddabbrev/ AbbSum302, AbbSum1543, AbbSum1575, AbbSum1481
	common /ddabbrev/ AbbSum2200, AbbSum2279, AbbSum2277, AbbSum2278
	common /ddabbrev/ AbbSum1544, AbbSum1772, AbbSum264, AbbSum1542
	common /ddabbrev/ AbbSum391, AbbSum993, AbbSum1915, AbbSum957
	common /ddabbrev/ AbbSum2054, AbbSum1107, AbbSum2100, AbbSum2098
	common /ddabbrev/ AbbSum2099, AbbSum996, AbbSum79, AbbSum1479
	common /ddabbrev/ AbbSum392, AbbSum753, AbbSum1106, AbbSum1477
	common /ddabbrev/ AbbSum995, AbbSum1417, AbbSum345, AbbSum394
	common /ddabbrev/ AbbSum370, AbbSum371, AbbSum395, AbbSum1917
	common /ddabbrev/ AbbSum679, AbbSum1222, AbbSum1882, AbbSum385
	common /ddabbrev/ AbbSum754, AbbSum1109, AbbSum999, AbbSum1483
	common /ddabbrev/ AbbSum460, AbbSum680, AbbSum2079, AbbSum2078
	common /ddabbrev/ AbbSum998, AbbSum1545, AbbSum2298, AbbSum2297
	common /ddabbrev/ AbbSum1395, AbbSum1805, AbbSum289, AbbSum1270
	common /ddabbrev/ AbbSum1482, AbbSum1546, AbbSum1111, AbbSum733
	common /ddabbrev/ AbbSum542, AbbSum1283, AbbSum249, AbbSum1284
	common /ddabbrev/ AbbSum876, AbbSum39, AbbSum799, AbbSum828
	common /ddabbrev/ AbbSum691, AbbSum1326, AbbSum1142, AbbSum1325
	common /ddabbrev/ AbbSum1143, AbbSum1223, AbbSum383, AbbSum2103
	common /ddabbrev/ AbbSum746, AbbSum875, AbbSum877, AbbSum1828
	common /ddabbrev/ AbbSum186, AbbSum185, AbbSum706, AbbSum707
	common /ddabbrev/ AbbSum624, AbbSum1803, AbbSum288, AbbSum598
	common /ddabbrev/ AbbSum825, AbbSum1327, AbbSum1308, AbbSum827
	common /ddabbrev/ AbbSum826, AbbSum2080, AbbSum960, AbbSum824
	common /ddabbrev/ AbbSum745, AbbSum2104, AbbSum1288, AbbSum664
	common /ddabbrev/ AbbSum184, AbbSum1804, AbbSum1512, AbbSum1513
	common /ddabbrev/ AbbSum1511, AbbSum1040, AbbSum1041, AbbSum1044
	common /ddabbrev/ AbbSum1328, AbbSum1881, AbbSum879, AbbSum556
	common /ddabbrev/ AbbSum15, AbbSum2281, AbbSum443, AbbSum822
	common /ddabbrev/ AbbSum855, AbbSum1331, AbbSum1330, AbbSum1578
	common /ddabbrev/ AbbSum2299, AbbSum1329, AbbSum126, AbbSum1351
	common /ddabbrev/ AbbSum1332, AbbSum513, AbbSum821, AbbSum1436
	common /ddabbrev/ AbbSum820, AbbSum1437, AbbSum1396, AbbSum1880
	common /ddabbrev/ AbbSum384, AbbSum625, AbbSum2280, AbbSum444
	common /ddabbrev/ AbbSum596, AbbSum2183, AbbSum1285, AbbSum1287
	common /ddabbrev/ AbbSum1869, AbbSum19, AbbSum16, AbbSum473
	common /ddabbrev/ AbbSum474, AbbSum287, AbbSum1286, AbbSum404
	common /ddabbrev/ AbbSum872, AbbSum2184, AbbSum871, AbbSum1094
	common /ddabbrev/ AbbSum35, AbbSum1095, AbbSum2097, AbbSum36
	common /ddabbrev/ AbbSum796, AbbSum797, AbbSum262, AbbSum1096
	common /ddabbrev/ AbbSum1097, AbbSum1098, AbbSum1190, AbbSum1099
	common /ddabbrev/ AbbSum1191, AbbSum241, AbbSum852, AbbSum853
	common /ddabbrev/ AbbSum1615, AbbSum40, AbbSum1001, AbbSum1616
	common /ddabbrev/ AbbSum1, AbbSum220, AbbSum900, AbbSum901
	common /ddabbrev/ AbbSum2055, AbbSum2056, AbbSum755, AbbSum756
	common /ddabbrev/ AbbSum196, AbbSum1825, AbbSum1209, AbbSum5
	common /ddabbrev/ AbbSum1210, AbbSum434, AbbSum435, AbbSum239
	common /ddabbrev/ AbbSum1186, AbbSum1211, AbbSum33, AbbSum1187
	common /ddabbrev/ AbbSum1212, AbbSum1213, AbbSum1188, AbbSum193
	common /ddabbrev/ AbbSum550, AbbSum551, AbbSum950, AbbSum6
	common /ddabbrev/ AbbSum1840, AbbSum260, AbbSum536, AbbSum537
	common /ddabbrev/ AbbSum34, AbbSum452, AbbSum453, AbbSum217

	double complex cint1, cint2, cint3, cint4, cint5, cint6, cint7
	double complex cint8, cint9, cint10, cint11, cint12, cint13
	double complex cint14, cint15, cint16, cint17, cint18, cint19
	double complex cint20, cint21, cint22, cint23, cint24, cint25
	double complex cint26, cint27, cint28, cint29, cint30, cint31
	double complex cint32, cint33, cint34, cint35, cint36, cint37
	double complex cint38, cint39, cint40, cint41, cint42, cint43
	double complex cint44, cint45, cint46, cint47, cint48, cint49
	double complex cint50, cint51, cint52, cint53, cint54, cint55
	double complex cint56, cint57, cint58, cint59, cint60
	common /ddloopint/ cint1, cint2, cint3, cint4, cint5, cint6
	common /ddloopint/ cint7, cint8, cint9, cint10, cint11, cint12
	common /ddloopint/ cint13, cint14, cint15, cint16, cint17
	common /ddloopint/ cint18, cint19, cint20, cint21, cint22
	common /ddloopint/ cint23, cint24, cint25, cint26, cint27
	common /ddloopint/ cint28, cint29, cint30, cint31, cint32
	common /ddloopint/ cint33, cint34, cint35, cint36, cint37
	common /ddloopint/ cint38, cint39, cint40, cint41, cint42
	common /ddloopint/ cint43, cint44, cint45, cint46, cint47
	common /ddloopint/ cint48, cint49, cint50, cint51, cint52
	common /ddloopint/ cint53, cint54, cint55, cint56, cint57
	common /ddloopint/ cint58, cint59, cint60

	integer*8 iint1, iint2, iint3, iint4, iint5, iint6, iint7, iint8
	integer*8 iint9, iint10, iint11, iint12, iint13, iint14, iint15
	integer*8 iint16, iint17, iint18, iint19, iint20, iint21, iint22
	common /ddloopint/ iint1, iint2, iint3, iint4, iint5, iint6
	common /ddloopint/ iint7, iint8, iint9, iint10, iint11, iint12
	common /ddloopint/ iint13, iint14, iint15, iint16, iint17
	common /ddloopint/ iint18, iint19, iint20, iint21, iint22

	double complex Cloop(48), Ctree(4), MatF(48,4)
	common /ddcoeff/ Cloop, Ctree, MatF
