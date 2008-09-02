
/*  copyright 1990  Richard H. Krukar all rights reserved

	Permission granted to buy, sell, or steal this software is granted.
    The author retains the right to distribute this software freely, but
    is not responsible for it's quality or maintainance. */

void idint(float *x, int length, float *wtab)
{
	register float tempr, tempi, temp2;
	switch (length) {
	    case 128 :
		tempr = x[64] - x[192];
		tempi = x[65] - x[193];

		x[64] += x[192];
		x[65] += x[193];

		x[192] = -tempi;  x[193] = tempr;
		tempr = x[0] - x[128];
		tempi = x[1] - x[129];

		x[0] += x[128];
		x[1] += x[129];

		x[128] = tempr + x[192];
		x[129] = tempi + x[193];

		temp2 = x[128]*wtab[0] - x[129]*wtab[1024]; x[129] = x[128]*wtab[1024] + x[129]*wtab[0]; x[128] = temp2;
		x[192] = tempr - x[192];
		x[193] = tempi - x[193];

		temp2 = x[192]*wtab[0] - x[193]*wtab[1024]; x[193] = x[192]*wtab[1024] + x[193]*wtab[0]; x[192] = temp2;
		tempr = x[66] - x[194];
		tempi = x[67] - x[195];

		x[66] += x[194];
		x[67] += x[195];

		x[194] = -tempi;  x[195] = tempr;
		tempr = x[2] - x[130];
		tempi = x[3] - x[131];

		x[2] += x[130];
		x[3] += x[131];

		x[130] = tempr + x[194];
		x[131] = tempi + x[195];

		temp2 = x[130]*wtab[32] - x[131]*wtab[992]; x[131] = x[130]*wtab[992] + x[131]*wtab[32]; x[130] = temp2;
		x[194] = tempr - x[194];
		x[195] = tempi - x[195];

		temp2 = x[194]*wtab[96] - x[195]*wtab[928]; x[195] = x[194]*wtab[928] + x[195]*wtab[96]; x[194] = temp2;
		tempr = x[68] - x[196];
		tempi = x[69] - x[197];

		x[68] += x[196];
		x[69] += x[197];

		x[196] = -tempi;  x[197] = tempr;
		tempr = x[4] - x[132];
		tempi = x[5] - x[133];

		x[4] += x[132];
		x[5] += x[133];

		x[132] = tempr + x[196];
		x[133] = tempi + x[197];

		temp2 = x[132]*wtab[64] - x[133]*wtab[960]; x[133] = x[132]*wtab[960] + x[133]*wtab[64]; x[132] = temp2;
		x[196] = tempr - x[196];
		x[197] = tempi - x[197];

		temp2 = x[196]*wtab[192] - x[197]*wtab[832]; x[197] = x[196]*wtab[832] + x[197]*wtab[192]; x[196] = temp2;
		tempr = x[70] - x[198];
		tempi = x[71] - x[199];

		x[70] += x[198];
		x[71] += x[199];

		x[198] = -tempi;  x[199] = tempr;
		tempr = x[6] - x[134];
		tempi = x[7] - x[135];

		x[6] += x[134];
		x[7] += x[135];

		x[134] = tempr + x[198];
		x[135] = tempi + x[199];

		temp2 = x[134]*wtab[96] - x[135]*wtab[928]; x[135] = x[134]*wtab[928] + x[135]*wtab[96]; x[134] = temp2;
		x[198] = tempr - x[198];
		x[199] = tempi - x[199];

		temp2 = x[198]*wtab[288] - x[199]*wtab[736]; x[199] = x[198]*wtab[736] + x[199]*wtab[288]; x[198] = temp2;
		tempr = x[72] - x[200];
		tempi = x[73] - x[201];

		x[72] += x[200];
		x[73] += x[201];

		x[200] = -tempi;  x[201] = tempr;
		tempr = x[8] - x[136];
		tempi = x[9] - x[137];

		x[8] += x[136];
		x[9] += x[137];

		x[136] = tempr + x[200];
		x[137] = tempi + x[201];

		temp2 = x[136]*wtab[128] - x[137]*wtab[896]; x[137] = x[136]*wtab[896] + x[137]*wtab[128]; x[136] = temp2;
		x[200] = tempr - x[200];
		x[201] = tempi - x[201];

		temp2 = x[200]*wtab[384] - x[201]*wtab[640]; x[201] = x[200]*wtab[640] + x[201]*wtab[384]; x[200] = temp2;
		tempr = x[74] - x[202];
		tempi = x[75] - x[203];

		x[74] += x[202];
		x[75] += x[203];

		x[202] = -tempi;  x[203] = tempr;
		tempr = x[10] - x[138];
		tempi = x[11] - x[139];

		x[10] += x[138];
		x[11] += x[139];

		x[138] = tempr + x[202];
		x[139] = tempi + x[203];

		temp2 = x[138]*wtab[160] - x[139]*wtab[864]; x[139] = x[138]*wtab[864] + x[139]*wtab[160]; x[138] = temp2;
		x[202] = tempr - x[202];
		x[203] = tempi - x[203];

		temp2 = x[202]*wtab[480] - x[203]*wtab[544]; x[203] = x[202]*wtab[544] + x[203]*wtab[480]; x[202] = temp2;
		tempr = x[76] - x[204];
		tempi = x[77] - x[205];

		x[76] += x[204];
		x[77] += x[205];

		x[204] = -tempi;  x[205] = tempr;
		tempr = x[12] - x[140];
		tempi = x[13] - x[141];

		x[12] += x[140];
		x[13] += x[141];

		x[140] = tempr + x[204];
		x[141] = tempi + x[205];

		temp2 = x[140]*wtab[192] - x[141]*wtab[832]; x[141] = x[140]*wtab[832] + x[141]*wtab[192]; x[140] = temp2;
		x[204] = tempr - x[204];
		x[205] = tempi - x[205];

		temp2 = x[204]*wtab[576] - x[205]*wtab[448]; x[205] = x[204]*wtab[448] + x[205]*wtab[576]; x[204] = temp2;
		tempr = x[78] - x[206];
		tempi = x[79] - x[207];

		x[78] += x[206];
		x[79] += x[207];

		x[206] = -tempi;  x[207] = tempr;
		tempr = x[14] - x[142];
		tempi = x[15] - x[143];

		x[14] += x[142];
		x[15] += x[143];

		x[142] = tempr + x[206];
		x[143] = tempi + x[207];

		temp2 = x[142]*wtab[224] - x[143]*wtab[800]; x[143] = x[142]*wtab[800] + x[143]*wtab[224]; x[142] = temp2;
		x[206] = tempr - x[206];
		x[207] = tempi - x[207];

		temp2 = x[206]*wtab[672] - x[207]*wtab[352]; x[207] = x[206]*wtab[352] + x[207]*wtab[672]; x[206] = temp2;
		tempr = x[80] - x[208];
		tempi = x[81] - x[209];

		x[80] += x[208];
		x[81] += x[209];

		x[208] = -tempi;  x[209] = tempr;
		tempr = x[16] - x[144];
		tempi = x[17] - x[145];

		x[16] += x[144];
		x[17] += x[145];

		x[144] = tempr + x[208];
		x[145] = tempi + x[209];

		temp2 = x[144]*wtab[256] - x[145]*wtab[768]; x[145] = x[144]*wtab[768] + x[145]*wtab[256]; x[144] = temp2;
		x[208] = tempr - x[208];
		x[209] = tempi - x[209];

		temp2 = x[208]*wtab[768] - x[209]*wtab[256]; x[209] = x[208]*wtab[256] + x[209]*wtab[768]; x[208] = temp2;
		tempr = x[82] - x[210];
		tempi = x[83] - x[211];

		x[82] += x[210];
		x[83] += x[211];

		x[210] = -tempi;  x[211] = tempr;
		tempr = x[18] - x[146];
		tempi = x[19] - x[147];

		x[18] += x[146];
		x[19] += x[147];

		x[146] = tempr + x[210];
		x[147] = tempi + x[211];

		temp2 = x[146]*wtab[288] - x[147]*wtab[736]; x[147] = x[146]*wtab[736] + x[147]*wtab[288]; x[146] = temp2;
		x[210] = tempr - x[210];
		x[211] = tempi - x[211];

		temp2 = x[210]*wtab[864] - x[211]*wtab[160]; x[211] = x[210]*wtab[160] + x[211]*wtab[864]; x[210] = temp2;
		tempr = x[84] - x[212];
		tempi = x[85] - x[213];

		x[84] += x[212];
		x[85] += x[213];

		x[212] = -tempi;  x[213] = tempr;
		tempr = x[20] - x[148];
		tempi = x[21] - x[149];

		x[20] += x[148];
		x[21] += x[149];

		x[148] = tempr + x[212];
		x[149] = tempi + x[213];

		temp2 = x[148]*wtab[320] - x[149]*wtab[704]; x[149] = x[148]*wtab[704] + x[149]*wtab[320]; x[148] = temp2;
		x[212] = tempr - x[212];
		x[213] = tempi - x[213];

		temp2 = x[212]*wtab[960] - x[213]*wtab[64]; x[213] = x[212]*wtab[64] + x[213]*wtab[960]; x[212] = temp2;
		tempr = x[86] - x[214];
		tempi = x[87] - x[215];

		x[86] += x[214];
		x[87] += x[215];

		x[214] = -tempi;  x[215] = tempr;
		tempr = x[22] - x[150];
		tempi = x[23] - x[151];

		x[22] += x[150];
		x[23] += x[151];

		x[150] = tempr + x[214];
		x[151] = tempi + x[215];

		temp2 = x[150]*wtab[352] - x[151]*wtab[672]; x[151] = x[150]*wtab[672] + x[151]*wtab[352]; x[150] = temp2;
		x[214] = tempr - x[214];
		x[215] = tempi - x[215];

		temp2 = x[214]*wtab[1056] - x[215]*wtab[32]; x[215] = x[214]*wtab[32] + x[215]*wtab[1056]; x[214] = temp2;
		tempr = x[88] - x[216];
		tempi = x[89] - x[217];

		x[88] += x[216];
		x[89] += x[217];

		x[216] = -tempi;  x[217] = tempr;
		tempr = x[24] - x[152];
		tempi = x[25] - x[153];

		x[24] += x[152];
		x[25] += x[153];

		x[152] = tempr + x[216];
		x[153] = tempi + x[217];

		temp2 = x[152]*wtab[384] - x[153]*wtab[640]; x[153] = x[152]*wtab[640] + x[153]*wtab[384]; x[152] = temp2;
		x[216] = tempr - x[216];
		x[217] = tempi - x[217];

		temp2 = x[216]*wtab[1152] - x[217]*wtab[128]; x[217] = x[216]*wtab[128] + x[217]*wtab[1152]; x[216] = temp2;
		tempr = x[90] - x[218];
		tempi = x[91] - x[219];

		x[90] += x[218];
		x[91] += x[219];

		x[218] = -tempi;  x[219] = tempr;
		tempr = x[26] - x[154];
		tempi = x[27] - x[155];

		x[26] += x[154];
		x[27] += x[155];

		x[154] = tempr + x[218];
		x[155] = tempi + x[219];

		temp2 = x[154]*wtab[416] - x[155]*wtab[608]; x[155] = x[154]*wtab[608] + x[155]*wtab[416]; x[154] = temp2;
		x[218] = tempr - x[218];
		x[219] = tempi - x[219];

		temp2 = x[218]*wtab[1248] - x[219]*wtab[224]; x[219] = x[218]*wtab[224] + x[219]*wtab[1248]; x[218] = temp2;
		tempr = x[92] - x[220];
		tempi = x[93] - x[221];

		x[92] += x[220];
		x[93] += x[221];

		x[220] = -tempi;  x[221] = tempr;
		tempr = x[28] - x[156];
		tempi = x[29] - x[157];

		x[28] += x[156];
		x[29] += x[157];

		x[156] = tempr + x[220];
		x[157] = tempi + x[221];

		temp2 = x[156]*wtab[448] - x[157]*wtab[576]; x[157] = x[156]*wtab[576] + x[157]*wtab[448]; x[156] = temp2;
		x[220] = tempr - x[220];
		x[221] = tempi - x[221];

		temp2 = x[220]*wtab[1344] - x[221]*wtab[320]; x[221] = x[220]*wtab[320] + x[221]*wtab[1344]; x[220] = temp2;
		tempr = x[94] - x[222];
		tempi = x[95] - x[223];

		x[94] += x[222];
		x[95] += x[223];

		x[222] = -tempi;  x[223] = tempr;
		tempr = x[30] - x[158];
		tempi = x[31] - x[159];

		x[30] += x[158];
		x[31] += x[159];

		x[158] = tempr + x[222];
		x[159] = tempi + x[223];

		temp2 = x[158]*wtab[480] - x[159]*wtab[544]; x[159] = x[158]*wtab[544] + x[159]*wtab[480]; x[158] = temp2;
		x[222] = tempr - x[222];
		x[223] = tempi - x[223];

		temp2 = x[222]*wtab[1440] - x[223]*wtab[416]; x[223] = x[222]*wtab[416] + x[223]*wtab[1440]; x[222] = temp2;
		tempr = x[96] - x[224];
		tempi = x[97] - x[225];

		x[96] += x[224];
		x[97] += x[225];

		x[224] = -tempi;  x[225] = tempr;
		tempr = x[32] - x[160];
		tempi = x[33] - x[161];

		x[32] += x[160];
		x[33] += x[161];

		x[160] = tempr + x[224];
		x[161] = tempi + x[225];

		temp2 = x[160]*wtab[512] - x[161]*wtab[512]; x[161] = x[160]*wtab[512] + x[161]*wtab[512]; x[160] = temp2;
		x[224] = tempr - x[224];
		x[225] = tempi - x[225];

		temp2 = x[224]*wtab[1536] - x[225]*wtab[512]; x[225] = x[224]*wtab[512] + x[225]*wtab[1536]; x[224] = temp2;
		tempr = x[98] - x[226];
		tempi = x[99] - x[227];

		x[98] += x[226];
		x[99] += x[227];

		x[226] = -tempi;  x[227] = tempr;
		tempr = x[34] - x[162];
		tempi = x[35] - x[163];

		x[34] += x[162];
		x[35] += x[163];

		x[162] = tempr + x[226];
		x[163] = tempi + x[227];

		temp2 = x[162]*wtab[544] - x[163]*wtab[480]; x[163] = x[162]*wtab[480] + x[163]*wtab[544]; x[162] = temp2;
		x[226] = tempr - x[226];
		x[227] = tempi - x[227];

		temp2 = x[226]*wtab[1632] - x[227]*wtab[608]; x[227] = x[226]*wtab[608] + x[227]*wtab[1632]; x[226] = temp2;
		tempr = x[100] - x[228];
		tempi = x[101] - x[229];

		x[100] += x[228];
		x[101] += x[229];

		x[228] = -tempi;  x[229] = tempr;
		tempr = x[36] - x[164];
		tempi = x[37] - x[165];

		x[36] += x[164];
		x[37] += x[165];

		x[164] = tempr + x[228];
		x[165] = tempi + x[229];

		temp2 = x[164]*wtab[576] - x[165]*wtab[448]; x[165] = x[164]*wtab[448] + x[165]*wtab[576]; x[164] = temp2;
		x[228] = tempr - x[228];
		x[229] = tempi - x[229];

		temp2 = x[228]*wtab[1728] - x[229]*wtab[704]; x[229] = x[228]*wtab[704] + x[229]*wtab[1728]; x[228] = temp2;
		tempr = x[102] - x[230];
		tempi = x[103] - x[231];

		x[102] += x[230];
		x[103] += x[231];

		x[230] = -tempi;  x[231] = tempr;
		tempr = x[38] - x[166];
		tempi = x[39] - x[167];

		x[38] += x[166];
		x[39] += x[167];

		x[166] = tempr + x[230];
		x[167] = tempi + x[231];

		temp2 = x[166]*wtab[608] - x[167]*wtab[416]; x[167] = x[166]*wtab[416] + x[167]*wtab[608]; x[166] = temp2;
		x[230] = tempr - x[230];
		x[231] = tempi - x[231];

		temp2 = x[230]*wtab[1824] - x[231]*wtab[800]; x[231] = x[230]*wtab[800] + x[231]*wtab[1824]; x[230] = temp2;
		tempr = x[104] - x[232];
		tempi = x[105] - x[233];

		x[104] += x[232];
		x[105] += x[233];

		x[232] = -tempi;  x[233] = tempr;
		tempr = x[40] - x[168];
		tempi = x[41] - x[169];

		x[40] += x[168];
		x[41] += x[169];

		x[168] = tempr + x[232];
		x[169] = tempi + x[233];

		temp2 = x[168]*wtab[640] - x[169]*wtab[384]; x[169] = x[168]*wtab[384] + x[169]*wtab[640]; x[168] = temp2;
		x[232] = tempr - x[232];
		x[233] = tempi - x[233];

		temp2 = x[232]*wtab[1920] - x[233]*wtab[896]; x[233] = x[232]*wtab[896] + x[233]*wtab[1920]; x[232] = temp2;
		tempr = x[106] - x[234];
		tempi = x[107] - x[235];

		x[106] += x[234];
		x[107] += x[235];

		x[234] = -tempi;  x[235] = tempr;
		tempr = x[42] - x[170];
		tempi = x[43] - x[171];

		x[42] += x[170];
		x[43] += x[171];

		x[170] = tempr + x[234];
		x[171] = tempi + x[235];

		temp2 = x[170]*wtab[672] - x[171]*wtab[352]; x[171] = x[170]*wtab[352] + x[171]*wtab[672]; x[170] = temp2;
		x[234] = tempr - x[234];
		x[235] = tempi - x[235];

		temp2 = x[234]*wtab[2016] - x[235]*wtab[992]; x[235] = x[234]*wtab[992] + x[235]*wtab[2016]; x[234] = temp2;
		tempr = x[108] - x[236];
		tempi = x[109] - x[237];

		x[108] += x[236];
		x[109] += x[237];

		x[236] = -tempi;  x[237] = tempr;
		tempr = x[44] - x[172];
		tempi = x[45] - x[173];

		x[44] += x[172];
		x[45] += x[173];

		x[172] = tempr + x[236];
		x[173] = tempi + x[237];

		temp2 = x[172]*wtab[704] - x[173]*wtab[320]; x[173] = x[172]*wtab[320] + x[173]*wtab[704]; x[172] = temp2;
		x[236] = tempr - x[236];
		x[237] = tempi - x[237];

		temp2 = x[236]*wtab[1984] - x[237]*wtab[1088]; x[237] = x[236]*wtab[1088] + x[237]*wtab[1984]; x[236] = temp2;
		tempr = x[110] - x[238];
		tempi = x[111] - x[239];

		x[110] += x[238];
		x[111] += x[239];

		x[238] = -tempi;  x[239] = tempr;
		tempr = x[46] - x[174];
		tempi = x[47] - x[175];

		x[46] += x[174];
		x[47] += x[175];

		x[174] = tempr + x[238];
		x[175] = tempi + x[239];

		temp2 = x[174]*wtab[736] - x[175]*wtab[288]; x[175] = x[174]*wtab[288] + x[175]*wtab[736]; x[174] = temp2;
		x[238] = tempr - x[238];
		x[239] = tempi - x[239];

		temp2 = x[238]*wtab[1888] - x[239]*wtab[1184]; x[239] = x[238]*wtab[1184] + x[239]*wtab[1888]; x[238] = temp2;
		tempr = x[112] - x[240];
		tempi = x[113] - x[241];

		x[112] += x[240];
		x[113] += x[241];

		x[240] = -tempi;  x[241] = tempr;
		tempr = x[48] - x[176];
		tempi = x[49] - x[177];

		x[48] += x[176];
		x[49] += x[177];

		x[176] = tempr + x[240];
		x[177] = tempi + x[241];

		temp2 = x[176]*wtab[768] - x[177]*wtab[256]; x[177] = x[176]*wtab[256] + x[177]*wtab[768]; x[176] = temp2;
		x[240] = tempr - x[240];
		x[241] = tempi - x[241];

		temp2 = x[240]*wtab[1792] - x[241]*wtab[1280]; x[241] = x[240]*wtab[1280] + x[241]*wtab[1792]; x[240] = temp2;
		tempr = x[114] - x[242];
		tempi = x[115] - x[243];

		x[114] += x[242];
		x[115] += x[243];

		x[242] = -tempi;  x[243] = tempr;
		tempr = x[50] - x[178];
		tempi = x[51] - x[179];

		x[50] += x[178];
		x[51] += x[179];

		x[178] = tempr + x[242];
		x[179] = tempi + x[243];

		temp2 = x[178]*wtab[800] - x[179]*wtab[224]; x[179] = x[178]*wtab[224] + x[179]*wtab[800]; x[178] = temp2;
		x[242] = tempr - x[242];
		x[243] = tempi - x[243];

		temp2 = x[242]*wtab[1696] - x[243]*wtab[1376]; x[243] = x[242]*wtab[1376] + x[243]*wtab[1696]; x[242] = temp2;
		tempr = x[116] - x[244];
		tempi = x[117] - x[245];

		x[116] += x[244];
		x[117] += x[245];

		x[244] = -tempi;  x[245] = tempr;
		tempr = x[52] - x[180];
		tempi = x[53] - x[181];

		x[52] += x[180];
		x[53] += x[181];

		x[180] = tempr + x[244];
		x[181] = tempi + x[245];

		temp2 = x[180]*wtab[832] - x[181]*wtab[192]; x[181] = x[180]*wtab[192] + x[181]*wtab[832]; x[180] = temp2;
		x[244] = tempr - x[244];
		x[245] = tempi - x[245];

		temp2 = x[244]*wtab[1600] - x[245]*wtab[1472]; x[245] = x[244]*wtab[1472] + x[245]*wtab[1600]; x[244] = temp2;
		tempr = x[118] - x[246];
		tempi = x[119] - x[247];

		x[118] += x[246];
		x[119] += x[247];

		x[246] = -tempi;  x[247] = tempr;
		tempr = x[54] - x[182];
		tempi = x[55] - x[183];

		x[54] += x[182];
		x[55] += x[183];

		x[182] = tempr + x[246];
		x[183] = tempi + x[247];

		temp2 = x[182]*wtab[864] - x[183]*wtab[160]; x[183] = x[182]*wtab[160] + x[183]*wtab[864]; x[182] = temp2;
		x[246] = tempr - x[246];
		x[247] = tempi - x[247];

		temp2 = x[246]*wtab[1504] - x[247]*wtab[1568]; x[247] = x[246]*wtab[1568] + x[247]*wtab[1504]; x[246] = temp2;
		tempr = x[120] - x[248];
		tempi = x[121] - x[249];

		x[120] += x[248];
		x[121] += x[249];

		x[248] = -tempi;  x[249] = tempr;
		tempr = x[56] - x[184];
		tempi = x[57] - x[185];

		x[56] += x[184];
		x[57] += x[185];

		x[184] = tempr + x[248];
		x[185] = tempi + x[249];

		temp2 = x[184]*wtab[896] - x[185]*wtab[128]; x[185] = x[184]*wtab[128] + x[185]*wtab[896]; x[184] = temp2;
		x[248] = tempr - x[248];
		x[249] = tempi - x[249];

		temp2 = x[248]*wtab[1408] - x[249]*wtab[1664]; x[249] = x[248]*wtab[1664] + x[249]*wtab[1408]; x[248] = temp2;
		tempr = x[122] - x[250];
		tempi = x[123] - x[251];

		x[122] += x[250];
		x[123] += x[251];

		x[250] = -tempi;  x[251] = tempr;
		tempr = x[58] - x[186];
		tempi = x[59] - x[187];

		x[58] += x[186];
		x[59] += x[187];

		x[186] = tempr + x[250];
		x[187] = tempi + x[251];

		temp2 = x[186]*wtab[928] - x[187]*wtab[96]; x[187] = x[186]*wtab[96] + x[187]*wtab[928]; x[186] = temp2;
		x[250] = tempr - x[250];
		x[251] = tempi - x[251];

		temp2 = x[250]*wtab[1312] - x[251]*wtab[1760]; x[251] = x[250]*wtab[1760] + x[251]*wtab[1312]; x[250] = temp2;
		tempr = x[124] - x[252];
		tempi = x[125] - x[253];

		x[124] += x[252];
		x[125] += x[253];

		x[252] = -tempi;  x[253] = tempr;
		tempr = x[60] - x[188];
		tempi = x[61] - x[189];

		x[60] += x[188];
		x[61] += x[189];

		x[188] = tempr + x[252];
		x[189] = tempi + x[253];

		temp2 = x[188]*wtab[960] - x[189]*wtab[64]; x[189] = x[188]*wtab[64] + x[189]*wtab[960]; x[188] = temp2;
		x[252] = tempr - x[252];
		x[253] = tempi - x[253];

		temp2 = x[252]*wtab[1216] - x[253]*wtab[1856]; x[253] = x[252]*wtab[1856] + x[253]*wtab[1216]; x[252] = temp2;
		tempr = x[126] - x[254];
		tempi = x[127] - x[255];

		x[126] += x[254];
		x[127] += x[255];

		x[254] = -tempi;  x[255] = tempr;
		tempr = x[62] - x[190];
		tempi = x[63] - x[191];

		x[62] += x[190];
		x[63] += x[191];

		x[190] = tempr + x[254];
		x[191] = tempi + x[255];

		temp2 = x[190]*wtab[992] - x[191]*wtab[32]; x[191] = x[190]*wtab[32] + x[191]*wtab[992]; x[190] = temp2;
		x[254] = tempr - x[254];
		x[255] = tempi - x[255];

		temp2 = x[254]*wtab[1120] - x[255]*wtab[1952]; x[255] = x[254]*wtab[1952] + x[255]*wtab[1120]; x[254] = temp2;
		tempr = x[144] - x[176];
		tempi = x[145] - x[177];

		x[144] += x[176];
		x[145] += x[177];

		x[176] = -tempi;  x[177] = tempr;
		tempr = x[128] - x[160];
		tempi = x[129] - x[161];

		x[128] += x[160];
		x[129] += x[161];

		x[160] = tempr + x[176];
		x[161] = tempi + x[177];

		temp2 = x[160]*wtab[0] - x[161]*wtab[1024]; x[161] = x[160]*wtab[1024] + x[161]*wtab[0]; x[160] = temp2;
		x[176] = tempr - x[176];
		x[177] = tempi - x[177];

		temp2 = x[176]*wtab[0] - x[177]*wtab[1024]; x[177] = x[176]*wtab[1024] + x[177]*wtab[0]; x[176] = temp2;
		tempr = x[146] - x[178];
		tempi = x[147] - x[179];

		x[146] += x[178];
		x[147] += x[179];

		x[178] = -tempi;  x[179] = tempr;
		tempr = x[130] - x[162];
		tempi = x[131] - x[163];

		x[130] += x[162];
		x[131] += x[163];

		x[162] = tempr + x[178];
		x[163] = tempi + x[179];

		temp2 = x[162]*wtab[128] - x[163]*wtab[896]; x[163] = x[162]*wtab[896] + x[163]*wtab[128]; x[162] = temp2;
		x[178] = tempr - x[178];
		x[179] = tempi - x[179];

		temp2 = x[178]*wtab[384] - x[179]*wtab[640]; x[179] = x[178]*wtab[640] + x[179]*wtab[384]; x[178] = temp2;
		tempr = x[148] - x[180];
		tempi = x[149] - x[181];

		x[148] += x[180];
		x[149] += x[181];

		x[180] = -tempi;  x[181] = tempr;
		tempr = x[132] - x[164];
		tempi = x[133] - x[165];

		x[132] += x[164];
		x[133] += x[165];

		x[164] = tempr + x[180];
		x[165] = tempi + x[181];

		temp2 = x[164]*wtab[256] - x[165]*wtab[768]; x[165] = x[164]*wtab[768] + x[165]*wtab[256]; x[164] = temp2;
		x[180] = tempr - x[180];
		x[181] = tempi - x[181];

		temp2 = x[180]*wtab[768] - x[181]*wtab[256]; x[181] = x[180]*wtab[256] + x[181]*wtab[768]; x[180] = temp2;
		tempr = x[150] - x[182];
		tempi = x[151] - x[183];

		x[150] += x[182];
		x[151] += x[183];

		x[182] = -tempi;  x[183] = tempr;
		tempr = x[134] - x[166];
		tempi = x[135] - x[167];

		x[134] += x[166];
		x[135] += x[167];

		x[166] = tempr + x[182];
		x[167] = tempi + x[183];

		temp2 = x[166]*wtab[384] - x[167]*wtab[640]; x[167] = x[166]*wtab[640] + x[167]*wtab[384]; x[166] = temp2;
		x[182] = tempr - x[182];
		x[183] = tempi - x[183];

		temp2 = x[182]*wtab[1152] - x[183]*wtab[128]; x[183] = x[182]*wtab[128] + x[183]*wtab[1152]; x[182] = temp2;
		tempr = x[152] - x[184];
		tempi = x[153] - x[185];

		x[152] += x[184];
		x[153] += x[185];

		x[184] = -tempi;  x[185] = tempr;
		tempr = x[136] - x[168];
		tempi = x[137] - x[169];

		x[136] += x[168];
		x[137] += x[169];

		x[168] = tempr + x[184];
		x[169] = tempi + x[185];

		temp2 = x[168]*wtab[512] - x[169]*wtab[512]; x[169] = x[168]*wtab[512] + x[169]*wtab[512]; x[168] = temp2;
		x[184] = tempr - x[184];
		x[185] = tempi - x[185];

		temp2 = x[184]*wtab[1536] - x[185]*wtab[512]; x[185] = x[184]*wtab[512] + x[185]*wtab[1536]; x[184] = temp2;
		tempr = x[154] - x[186];
		tempi = x[155] - x[187];

		x[154] += x[186];
		x[155] += x[187];

		x[186] = -tempi;  x[187] = tempr;
		tempr = x[138] - x[170];
		tempi = x[139] - x[171];

		x[138] += x[170];
		x[139] += x[171];

		x[170] = tempr + x[186];
		x[171] = tempi + x[187];

		temp2 = x[170]*wtab[640] - x[171]*wtab[384]; x[171] = x[170]*wtab[384] + x[171]*wtab[640]; x[170] = temp2;
		x[186] = tempr - x[186];
		x[187] = tempi - x[187];

		temp2 = x[186]*wtab[1920] - x[187]*wtab[896]; x[187] = x[186]*wtab[896] + x[187]*wtab[1920]; x[186] = temp2;
		tempr = x[156] - x[188];
		tempi = x[157] - x[189];

		x[156] += x[188];
		x[157] += x[189];

		x[188] = -tempi;  x[189] = tempr;
		tempr = x[140] - x[172];
		tempi = x[141] - x[173];

		x[140] += x[172];
		x[141] += x[173];

		x[172] = tempr + x[188];
		x[173] = tempi + x[189];

		temp2 = x[172]*wtab[768] - x[173]*wtab[256]; x[173] = x[172]*wtab[256] + x[173]*wtab[768]; x[172] = temp2;
		x[188] = tempr - x[188];
		x[189] = tempi - x[189];

		temp2 = x[188]*wtab[1792] - x[189]*wtab[1280]; x[189] = x[188]*wtab[1280] + x[189]*wtab[1792]; x[188] = temp2;
		tempr = x[158] - x[190];
		tempi = x[159] - x[191];

		x[158] += x[190];
		x[159] += x[191];

		x[190] = -tempi;  x[191] = tempr;
		tempr = x[142] - x[174];
		tempi = x[143] - x[175];

		x[142] += x[174];
		x[143] += x[175];

		x[174] = tempr + x[190];
		x[175] = tempi + x[191];

		temp2 = x[174]*wtab[896] - x[175]*wtab[128]; x[175] = x[174]*wtab[128] + x[175]*wtab[896]; x[174] = temp2;
		x[190] = tempr - x[190];
		x[191] = tempi - x[191];

		temp2 = x[190]*wtab[1408] - x[191]*wtab[1664]; x[191] = x[190]*wtab[1664] + x[191]*wtab[1408]; x[190] = temp2;
		tempr = x[164] - x[172];
		tempi = x[165] - x[173];

		x[164] += x[172];
		x[165] += x[173];

		x[172] = -tempi;  x[173] = tempr;
		tempr = x[160] - x[168];
		tempi = x[161] - x[169];

		x[160] += x[168];
		x[161] += x[169];

		x[168] = tempr + x[172];
		x[169] = tempi + x[173];

		temp2 = x[168]*wtab[0] - x[169]*wtab[1024]; x[169] = x[168]*wtab[1024] + x[169]*wtab[0]; x[168] = temp2;
		x[172] = tempr - x[172];
		x[173] = tempi - x[173];

		temp2 = x[172]*wtab[0] - x[173]*wtab[1024]; x[173] = x[172]*wtab[1024] + x[173]*wtab[0]; x[172] = temp2;
		tempr = x[166] - x[174];
		tempi = x[167] - x[175];

		x[166] += x[174];
		x[167] += x[175];

		x[174] = -tempi;  x[175] = tempr;
		tempr = x[162] - x[170];
		tempi = x[163] - x[171];

		x[162] += x[170];
		x[163] += x[171];

		x[170] = tempr + x[174];
		x[171] = tempi + x[175];

		temp2 = x[170]*wtab[512] - x[171]*wtab[512]; x[171] = x[170]*wtab[512] + x[171]*wtab[512]; x[170] = temp2;
		x[174] = tempr - x[174];
		x[175] = tempi - x[175];

		temp2 = x[174]*wtab[1536] - x[175]*wtab[512]; x[175] = x[174]*wtab[512] + x[175]*wtab[1536]; x[174] = temp2;
		tempr = x[168] - x[170];
		tempi = x[169] - x[171];

		x[168] += x[170];
		x[169] += x[171];

		x[170] = tempr;
		x[171] = tempi;

		tempr = x[172] - x[174];
		tempi = x[173] - x[175];

		x[172] += x[174];
		x[173] += x[175];

		x[174] = tempr;
		x[175] = tempi;

		tempr = x[160] - x[164];
		tempi = x[161] - x[165];

		x[160] += x[164];
		x[161] += x[165];

		x[164] = tempr;
		x[165] = tempi;

		tempr = x[162] - x[166];
		tempi = x[163] - x[167];

		x[162] += x[166];
		x[163] += x[167];

		temp2 = -tempi; tempi = tempr; tempr=temp2;
		x[166] = x[164] - tempr;
		x[167] = x[165] - tempi;

		x[164] += tempr;
		x[165] += tempi;

		tempr = x[160] - x[162];
		tempi = x[161] - x[163];

		x[160] += x[162];
		x[161] += x[163];

		x[162] = tempr;
		x[163] = tempi;

		tempr = x[180] - x[188];
		tempi = x[181] - x[189];

		x[180] += x[188];
		x[181] += x[189];

		x[188] = -tempi;  x[189] = tempr;
		tempr = x[176] - x[184];
		tempi = x[177] - x[185];

		x[176] += x[184];
		x[177] += x[185];

		x[184] = tempr + x[188];
		x[185] = tempi + x[189];

		temp2 = x[184]*wtab[0] - x[185]*wtab[1024]; x[185] = x[184]*wtab[1024] + x[185]*wtab[0]; x[184] = temp2;
		x[188] = tempr - x[188];
		x[189] = tempi - x[189];

		temp2 = x[188]*wtab[0] - x[189]*wtab[1024]; x[189] = x[188]*wtab[1024] + x[189]*wtab[0]; x[188] = temp2;
		tempr = x[182] - x[190];
		tempi = x[183] - x[191];

		x[182] += x[190];
		x[183] += x[191];

		x[190] = -tempi;  x[191] = tempr;
		tempr = x[178] - x[186];
		tempi = x[179] - x[187];

		x[178] += x[186];
		x[179] += x[187];

		x[186] = tempr + x[190];
		x[187] = tempi + x[191];

		temp2 = x[186]*wtab[512] - x[187]*wtab[512]; x[187] = x[186]*wtab[512] + x[187]*wtab[512]; x[186] = temp2;
		x[190] = tempr - x[190];
		x[191] = tempi - x[191];

		temp2 = x[190]*wtab[1536] - x[191]*wtab[512]; x[191] = x[190]*wtab[512] + x[191]*wtab[1536]; x[190] = temp2;
		tempr = x[184] - x[186];
		tempi = x[185] - x[187];

		x[184] += x[186];
		x[185] += x[187];

		x[186] = tempr;
		x[187] = tempi;

		tempr = x[188] - x[190];
		tempi = x[189] - x[191];

		x[188] += x[190];
		x[189] += x[191];

		x[190] = tempr;
		x[191] = tempi;

		tempr = x[176] - x[180];
		tempi = x[177] - x[181];

		x[176] += x[180];
		x[177] += x[181];

		x[180] = tempr;
		x[181] = tempi;

		tempr = x[178] - x[182];
		tempi = x[179] - x[183];

		x[178] += x[182];
		x[179] += x[183];

		temp2 = -tempi; tempi = tempr; tempr=temp2;
		x[182] = x[180] - tempr;
		x[183] = x[181] - tempi;

		x[180] += tempr;
		x[181] += tempi;

		tempr = x[176] - x[178];
		tempi = x[177] - x[179];

		x[176] += x[178];
		x[177] += x[179];

		x[178] = tempr;
		x[179] = tempi;

		tempr = x[136] - x[152];
		tempi = x[137] - x[153];

		x[136] += x[152];
		x[137] += x[153];

		x[152] = -tempi;  x[153] = tempr;
		tempr = x[128] - x[144];
		tempi = x[129] - x[145];

		x[128] += x[144];
		x[129] += x[145];

		x[144] = tempr + x[152];
		x[145] = tempi + x[153];

		temp2 = x[144]*wtab[0] - x[145]*wtab[1024]; x[145] = x[144]*wtab[1024] + x[145]*wtab[0]; x[144] = temp2;
		x[152] = tempr - x[152];
		x[153] = tempi - x[153];

		temp2 = x[152]*wtab[0] - x[153]*wtab[1024]; x[153] = x[152]*wtab[1024] + x[153]*wtab[0]; x[152] = temp2;
		tempr = x[138] - x[154];
		tempi = x[139] - x[155];

		x[138] += x[154];
		x[139] += x[155];

		x[154] = -tempi;  x[155] = tempr;
		tempr = x[130] - x[146];
		tempi = x[131] - x[147];

		x[130] += x[146];
		x[131] += x[147];

		x[146] = tempr + x[154];
		x[147] = tempi + x[155];

		temp2 = x[146]*wtab[256] - x[147]*wtab[768]; x[147] = x[146]*wtab[768] + x[147]*wtab[256]; x[146] = temp2;
		x[154] = tempr - x[154];
		x[155] = tempi - x[155];

		temp2 = x[154]*wtab[768] - x[155]*wtab[256]; x[155] = x[154]*wtab[256] + x[155]*wtab[768]; x[154] = temp2;
		tempr = x[140] - x[156];
		tempi = x[141] - x[157];

		x[140] += x[156];
		x[141] += x[157];

		x[156] = -tempi;  x[157] = tempr;
		tempr = x[132] - x[148];
		tempi = x[133] - x[149];

		x[132] += x[148];
		x[133] += x[149];

		x[148] = tempr + x[156];
		x[149] = tempi + x[157];

		temp2 = x[148]*wtab[512] - x[149]*wtab[512]; x[149] = x[148]*wtab[512] + x[149]*wtab[512]; x[148] = temp2;
		x[156] = tempr - x[156];
		x[157] = tempi - x[157];

		temp2 = x[156]*wtab[1536] - x[157]*wtab[512]; x[157] = x[156]*wtab[512] + x[157]*wtab[1536]; x[156] = temp2;
		tempr = x[142] - x[158];
		tempi = x[143] - x[159];

		x[142] += x[158];
		x[143] += x[159];

		x[158] = -tempi;  x[159] = tempr;
		tempr = x[134] - x[150];
		tempi = x[135] - x[151];

		x[134] += x[150];
		x[135] += x[151];

		x[150] = tempr + x[158];
		x[151] = tempi + x[159];

		temp2 = x[150]*wtab[768] - x[151]*wtab[256]; x[151] = x[150]*wtab[256] + x[151]*wtab[768]; x[150] = temp2;
		x[158] = tempr - x[158];
		x[159] = tempi - x[159];

		temp2 = x[158]*wtab[1792] - x[159]*wtab[1280]; x[159] = x[158]*wtab[1280] + x[159]*wtab[1792]; x[158] = temp2;
		tempr = x[144] - x[148];
		tempi = x[145] - x[149];

		x[144] += x[148];
		x[145] += x[149];

		x[148] = tempr;
		x[149] = tempi;

		tempr = x[146] - x[150];
		tempi = x[147] - x[151];

		x[146] += x[150];
		x[147] += x[151];

		temp2 = -tempi; tempi = tempr; tempr=temp2;
		x[150] = x[148] - tempr;
		x[151] = x[149] - tempi;

		x[148] += tempr;
		x[149] += tempi;

		tempr = x[144] - x[146];
		tempi = x[145] - x[147];

		x[144] += x[146];
		x[145] += x[147];

		x[146] = tempr;
		x[147] = tempi;

		tempr = x[152] - x[156];
		tempi = x[153] - x[157];

		x[152] += x[156];
		x[153] += x[157];

		x[156] = tempr;
		x[157] = tempi;

		tempr = x[154] - x[158];
		tempi = x[155] - x[159];

		x[154] += x[158];
		x[155] += x[159];

		temp2 = -tempi; tempi = tempr; tempr=temp2;
		x[158] = x[156] - tempr;
		x[159] = x[157] - tempi;

		x[156] += tempr;
		x[157] += tempi;

		tempr = x[152] - x[154];
		tempi = x[153] - x[155];

		x[152] += x[154];
		x[153] += x[155];

		x[154] = tempr;
		x[155] = tempi;

		tempr = x[132] - x[140];
		tempi = x[133] - x[141];

		x[132] += x[140];
		x[133] += x[141];

		x[140] = -tempi;  x[141] = tempr;
		tempr = x[128] - x[136];
		tempi = x[129] - x[137];

		x[128] += x[136];
		x[129] += x[137];

		x[136] = tempr + x[140];
		x[137] = tempi + x[141];

		temp2 = x[136]*wtab[0] - x[137]*wtab[1024]; x[137] = x[136]*wtab[1024] + x[137]*wtab[0]; x[136] = temp2;
		x[140] = tempr - x[140];
		x[141] = tempi - x[141];

		temp2 = x[140]*wtab[0] - x[141]*wtab[1024]; x[141] = x[140]*wtab[1024] + x[141]*wtab[0]; x[140] = temp2;
		tempr = x[134] - x[142];
		tempi = x[135] - x[143];

		x[134] += x[142];
		x[135] += x[143];

		x[142] = -tempi;  x[143] = tempr;
		tempr = x[130] - x[138];
		tempi = x[131] - x[139];

		x[130] += x[138];
		x[131] += x[139];

		x[138] = tempr + x[142];
		x[139] = tempi + x[143];

		temp2 = x[138]*wtab[512] - x[139]*wtab[512]; x[139] = x[138]*wtab[512] + x[139]*wtab[512]; x[138] = temp2;
		x[142] = tempr - x[142];
		x[143] = tempi - x[143];

		temp2 = x[142]*wtab[1536] - x[143]*wtab[512]; x[143] = x[142]*wtab[512] + x[143]*wtab[1536]; x[142] = temp2;
		tempr = x[136] - x[138];
		tempi = x[137] - x[139];

		x[136] += x[138];
		x[137] += x[139];

		x[138] = tempr;
		x[139] = tempi;

		tempr = x[140] - x[142];
		tempi = x[141] - x[143];

		x[140] += x[142];
		x[141] += x[143];

		x[142] = tempr;
		x[143] = tempi;

		tempr = x[128] - x[132];
		tempi = x[129] - x[133];

		x[128] += x[132];
		x[129] += x[133];

		x[132] = tempr;
		x[133] = tempi;

		tempr = x[130] - x[134];
		tempi = x[131] - x[135];

		x[130] += x[134];
		x[131] += x[135];

		temp2 = -tempi; tempi = tempr; tempr=temp2;
		x[134] = x[132] - tempr;
		x[135] = x[133] - tempi;

		x[132] += tempr;
		x[133] += tempi;

		tempr = x[128] - x[130];
		tempi = x[129] - x[131];

		x[128] += x[130];
		x[129] += x[131];

		x[130] = tempr;
		x[131] = tempi;

		tempr = x[208] - x[240];
		tempi = x[209] - x[241];

		x[208] += x[240];
		x[209] += x[241];

		x[240] = -tempi;  x[241] = tempr;
		tempr = x[192] - x[224];
		tempi = x[193] - x[225];

		x[192] += x[224];
		x[193] += x[225];

		x[224] = tempr + x[240];
		x[225] = tempi + x[241];

		temp2 = x[224]*wtab[0] - x[225]*wtab[1024]; x[225] = x[224]*wtab[1024] + x[225]*wtab[0]; x[224] = temp2;
		x[240] = tempr - x[240];
		x[241] = tempi - x[241];

		temp2 = x[240]*wtab[0] - x[241]*wtab[1024]; x[241] = x[240]*wtab[1024] + x[241]*wtab[0]; x[240] = temp2;
		tempr = x[210] - x[242];
		tempi = x[211] - x[243];

		x[210] += x[242];
		x[211] += x[243];

		x[242] = -tempi;  x[243] = tempr;
		tempr = x[194] - x[226];
		tempi = x[195] - x[227];

		x[194] += x[226];
		x[195] += x[227];

		x[226] = tempr + x[242];
		x[227] = tempi + x[243];

		temp2 = x[226]*wtab[128] - x[227]*wtab[896]; x[227] = x[226]*wtab[896] + x[227]*wtab[128]; x[226] = temp2;
		x[242] = tempr - x[242];
		x[243] = tempi - x[243];

		temp2 = x[242]*wtab[384] - x[243]*wtab[640]; x[243] = x[242]*wtab[640] + x[243]*wtab[384]; x[242] = temp2;
		tempr = x[212] - x[244];
		tempi = x[213] - x[245];

		x[212] += x[244];
		x[213] += x[245];

		x[244] = -tempi;  x[245] = tempr;
		tempr = x[196] - x[228];
		tempi = x[197] - x[229];

		x[196] += x[228];
		x[197] += x[229];

		x[228] = tempr + x[244];
		x[229] = tempi + x[245];

		temp2 = x[228]*wtab[256] - x[229]*wtab[768]; x[229] = x[228]*wtab[768] + x[229]*wtab[256]; x[228] = temp2;
		x[244] = tempr - x[244];
		x[245] = tempi - x[245];

		temp2 = x[244]*wtab[768] - x[245]*wtab[256]; x[245] = x[244]*wtab[256] + x[245]*wtab[768]; x[244] = temp2;
		tempr = x[214] - x[246];
		tempi = x[215] - x[247];

		x[214] += x[246];
		x[215] += x[247];

		x[246] = -tempi;  x[247] = tempr;
		tempr = x[198] - x[230];
		tempi = x[199] - x[231];

		x[198] += x[230];
		x[199] += x[231];

		x[230] = tempr + x[246];
		x[231] = tempi + x[247];

		temp2 = x[230]*wtab[384] - x[231]*wtab[640]; x[231] = x[230]*wtab[640] + x[231]*wtab[384]; x[230] = temp2;
		x[246] = tempr - x[246];
		x[247] = tempi - x[247];

		temp2 = x[246]*wtab[1152] - x[247]*wtab[128]; x[247] = x[246]*wtab[128] + x[247]*wtab[1152]; x[246] = temp2;
		tempr = x[216] - x[248];
		tempi = x[217] - x[249];

		x[216] += x[248];
		x[217] += x[249];

		x[248] = -tempi;  x[249] = tempr;
		tempr = x[200] - x[232];
		tempi = x[201] - x[233];

		x[200] += x[232];
		x[201] += x[233];

		x[232] = tempr + x[248];
		x[233] = tempi + x[249];

		temp2 = x[232]*wtab[512] - x[233]*wtab[512]; x[233] = x[232]*wtab[512] + x[233]*wtab[512]; x[232] = temp2;
		x[248] = tempr - x[248];
		x[249] = tempi - x[249];

		temp2 = x[248]*wtab[1536] - x[249]*wtab[512]; x[249] = x[248]*wtab[512] + x[249]*wtab[1536]; x[248] = temp2;
		tempr = x[218] - x[250];
		tempi = x[219] - x[251];

		x[218] += x[250];
		x[219] += x[251];

		x[250] = -tempi;  x[251] = tempr;
		tempr = x[202] - x[234];
		tempi = x[203] - x[235];

		x[202] += x[234];
		x[203] += x[235];

		x[234] = tempr + x[250];
		x[235] = tempi + x[251];

		temp2 = x[234]*wtab[640] - x[235]*wtab[384]; x[235] = x[234]*wtab[384] + x[235]*wtab[640]; x[234] = temp2;
		x[250] = tempr - x[250];
		x[251] = tempi - x[251];

		temp2 = x[250]*wtab[1920] - x[251]*wtab[896]; x[251] = x[250]*wtab[896] + x[251]*wtab[1920]; x[250] = temp2;
		tempr = x[220] - x[252];
		tempi = x[221] - x[253];

		x[220] += x[252];
		x[221] += x[253];

		x[252] = -tempi;  x[253] = tempr;
		tempr = x[204] - x[236];
		tempi = x[205] - x[237];

		x[204] += x[236];
		x[205] += x[237];

		x[236] = tempr + x[252];
		x[237] = tempi + x[253];

		temp2 = x[236]*wtab[768] - x[237]*wtab[256]; x[237] = x[236]*wtab[256] + x[237]*wtab[768]; x[236] = temp2;
		x[252] = tempr - x[252];
		x[253] = tempi - x[253];

		temp2 = x[252]*wtab[1792] - x[253]*wtab[1280]; x[253] = x[252]*wtab[1280] + x[253]*wtab[1792]; x[252] = temp2;
		tempr = x[222] - x[254];
		tempi = x[223] - x[255];

		x[222] += x[254];
		x[223] += x[255];

		x[254] = -tempi;  x[255] = tempr;
		tempr = x[206] - x[238];
		tempi = x[207] - x[239];

		x[206] += x[238];
		x[207] += x[239];

		x[238] = tempr + x[254];
		x[239] = tempi + x[255];

		temp2 = x[238]*wtab[896] - x[239]*wtab[128]; x[239] = x[238]*wtab[128] + x[239]*wtab[896]; x[238] = temp2;
		x[254] = tempr - x[254];
		x[255] = tempi - x[255];

		temp2 = x[254]*wtab[1408] - x[255]*wtab[1664]; x[255] = x[254]*wtab[1664] + x[255]*wtab[1408]; x[254] = temp2;
		tempr = x[228] - x[236];
		tempi = x[229] - x[237];

		x[228] += x[236];
		x[229] += x[237];

		x[236] = -tempi;  x[237] = tempr;
		tempr = x[224] - x[232];
		tempi = x[225] - x[233];

		x[224] += x[232];
		x[225] += x[233];

		x[232] = tempr + x[236];
		x[233] = tempi + x[237];

		temp2 = x[232]*wtab[0] - x[233]*wtab[1024]; x[233] = x[232]*wtab[1024] + x[233]*wtab[0]; x[232] = temp2;
		x[236] = tempr - x[236];
		x[237] = tempi - x[237];

		temp2 = x[236]*wtab[0] - x[237]*wtab[1024]; x[237] = x[236]*wtab[1024] + x[237]*wtab[0]; x[236] = temp2;
		tempr = x[230] - x[238];
		tempi = x[231] - x[239];

		x[230] += x[238];
		x[231] += x[239];

		x[238] = -tempi;  x[239] = tempr;
		tempr = x[226] - x[234];
		tempi = x[227] - x[235];

		x[226] += x[234];
		x[227] += x[235];

		x[234] = tempr + x[238];
		x[235] = tempi + x[239];

		temp2 = x[234]*wtab[512] - x[235]*wtab[512]; x[235] = x[234]*wtab[512] + x[235]*wtab[512]; x[234] = temp2;
		x[238] = tempr - x[238];
		x[239] = tempi - x[239];

		temp2 = x[238]*wtab[1536] - x[239]*wtab[512]; x[239] = x[238]*wtab[512] + x[239]*wtab[1536]; x[238] = temp2;
		tempr = x[232] - x[234];
		tempi = x[233] - x[235];

		x[232] += x[234];
		x[233] += x[235];

		x[234] = tempr;
		x[235] = tempi;

		tempr = x[236] - x[238];
		tempi = x[237] - x[239];

		x[236] += x[238];
		x[237] += x[239];

		x[238] = tempr;
		x[239] = tempi;

		tempr = x[224] - x[228];
		tempi = x[225] - x[229];

		x[224] += x[228];
		x[225] += x[229];

		x[228] = tempr;
		x[229] = tempi;

		tempr = x[226] - x[230];
		tempi = x[227] - x[231];

		x[226] += x[230];
		x[227] += x[231];

		temp2 = -tempi; tempi = tempr; tempr=temp2;
		x[230] = x[228] - tempr;
		x[231] = x[229] - tempi;

		x[228] += tempr;
		x[229] += tempi;

		tempr = x[224] - x[226];
		tempi = x[225] - x[227];

		x[224] += x[226];
		x[225] += x[227];

		x[226] = tempr;
		x[227] = tempi;

		tempr = x[244] - x[252];
		tempi = x[245] - x[253];

		x[244] += x[252];
		x[245] += x[253];

		x[252] = -tempi;  x[253] = tempr;
		tempr = x[240] - x[248];
		tempi = x[241] - x[249];

		x[240] += x[248];
		x[241] += x[249];

		x[248] = tempr + x[252];
		x[249] = tempi + x[253];

		temp2 = x[248]*wtab[0] - x[249]*wtab[1024]; x[249] = x[248]*wtab[1024] + x[249]*wtab[0]; x[248] = temp2;
		x[252] = tempr - x[252];
		x[253] = tempi - x[253];

		temp2 = x[252]*wtab[0] - x[253]*wtab[1024]; x[253] = x[252]*wtab[1024] + x[253]*wtab[0]; x[252] = temp2;
		tempr = x[246] - x[254];
		tempi = x[247] - x[255];

		x[246] += x[254];
		x[247] += x[255];

		x[254] = -tempi;  x[255] = tempr;
		tempr = x[242] - x[250];
		tempi = x[243] - x[251];

		x[242] += x[250];
		x[243] += x[251];

		x[250] = tempr + x[254];
		x[251] = tempi + x[255];

		temp2 = x[250]*wtab[512] - x[251]*wtab[512]; x[251] = x[250]*wtab[512] + x[251]*wtab[512]; x[250] = temp2;
		x[254] = tempr - x[254];
		x[255] = tempi - x[255];

		temp2 = x[254]*wtab[1536] - x[255]*wtab[512]; x[255] = x[254]*wtab[512] + x[255]*wtab[1536]; x[254] = temp2;
		tempr = x[248] - x[250];
		tempi = x[249] - x[251];

		x[248] += x[250];
		x[249] += x[251];

		x[250] = tempr;
		x[251] = tempi;

		tempr = x[252] - x[254];
		tempi = x[253] - x[255];

		x[252] += x[254];
		x[253] += x[255];

		x[254] = tempr;
		x[255] = tempi;

		tempr = x[240] - x[244];
		tempi = x[241] - x[245];

		x[240] += x[244];
		x[241] += x[245];

		x[244] = tempr;
		x[245] = tempi;

		tempr = x[242] - x[246];
		tempi = x[243] - x[247];

		x[242] += x[246];
		x[243] += x[247];

		temp2 = -tempi; tempi = tempr; tempr=temp2;
		x[246] = x[244] - tempr;
		x[247] = x[245] - tempi;

		x[244] += tempr;
		x[245] += tempi;

		tempr = x[240] - x[242];
		tempi = x[241] - x[243];

		x[240] += x[242];
		x[241] += x[243];

		x[242] = tempr;
		x[243] = tempi;

		tempr = x[200] - x[216];
		tempi = x[201] - x[217];

		x[200] += x[216];
		x[201] += x[217];

		x[216] = -tempi;  x[217] = tempr;
		tempr = x[192] - x[208];
		tempi = x[193] - x[209];

		x[192] += x[208];
		x[193] += x[209];

		x[208] = tempr + x[216];
		x[209] = tempi + x[217];

		temp2 = x[208]*wtab[0] - x[209]*wtab[1024]; x[209] = x[208]*wtab[1024] + x[209]*wtab[0]; x[208] = temp2;
		x[216] = tempr - x[216];
		x[217] = tempi - x[217];

		temp2 = x[216]*wtab[0] - x[217]*wtab[1024]; x[217] = x[216]*wtab[1024] + x[217]*wtab[0]; x[216] = temp2;
		tempr = x[202] - x[218];
		tempi = x[203] - x[219];

		x[202] += x[218];
		x[203] += x[219];

		x[218] = -tempi;  x[219] = tempr;
		tempr = x[194] - x[210];
		tempi = x[195] - x[211];

		x[194] += x[210];
		x[195] += x[211];

		x[210] = tempr + x[218];
		x[211] = tempi + x[219];

		temp2 = x[210]*wtab[256] - x[211]*wtab[768]; x[211] = x[210]*wtab[768] + x[211]*wtab[256]; x[210] = temp2;
		x[218] = tempr - x[218];
		x[219] = tempi - x[219];

		temp2 = x[218]*wtab[768] - x[219]*wtab[256]; x[219] = x[218]*wtab[256] + x[219]*wtab[768]; x[218] = temp2;
		tempr = x[204] - x[220];
		tempi = x[205] - x[221];

		x[204] += x[220];
		x[205] += x[221];

		x[220] = -tempi;  x[221] = tempr;
		tempr = x[196] - x[212];
		tempi = x[197] - x[213];

		x[196] += x[212];
		x[197] += x[213];

		x[212] = tempr + x[220];
		x[213] = tempi + x[221];

		temp2 = x[212]*wtab[512] - x[213]*wtab[512]; x[213] = x[212]*wtab[512] + x[213]*wtab[512]; x[212] = temp2;
		x[220] = tempr - x[220];
		x[221] = tempi - x[221];

		temp2 = x[220]*wtab[1536] - x[221]*wtab[512]; x[221] = x[220]*wtab[512] + x[221]*wtab[1536]; x[220] = temp2;
		tempr = x[206] - x[222];
		tempi = x[207] - x[223];

		x[206] += x[222];
		x[207] += x[223];

		x[222] = -tempi;  x[223] = tempr;
		tempr = x[198] - x[214];
		tempi = x[199] - x[215];

		x[198] += x[214];
		x[199] += x[215];

		x[214] = tempr + x[222];
		x[215] = tempi + x[223];

		temp2 = x[214]*wtab[768] - x[215]*wtab[256]; x[215] = x[214]*wtab[256] + x[215]*wtab[768]; x[214] = temp2;
		x[222] = tempr - x[222];
		x[223] = tempi - x[223];

		temp2 = x[222]*wtab[1792] - x[223]*wtab[1280]; x[223] = x[222]*wtab[1280] + x[223]*wtab[1792]; x[222] = temp2;
		tempr = x[208] - x[212];
		tempi = x[209] - x[213];

		x[208] += x[212];
		x[209] += x[213];

		x[212] = tempr;
		x[213] = tempi;

		tempr = x[210] - x[214];
		tempi = x[211] - x[215];

		x[210] += x[214];
		x[211] += x[215];

		temp2 = -tempi; tempi = tempr; tempr=temp2;
		x[214] = x[212] - tempr;
		x[215] = x[213] - tempi;

		x[212] += tempr;
		x[213] += tempi;

		tempr = x[208] - x[210];
		tempi = x[209] - x[211];

		x[208] += x[210];
		x[209] += x[211];

		x[210] = tempr;
		x[211] = tempi;

		tempr = x[216] - x[220];
		tempi = x[217] - x[221];

		x[216] += x[220];
		x[217] += x[221];

		x[220] = tempr;
		x[221] = tempi;

		tempr = x[218] - x[222];
		tempi = x[219] - x[223];

		x[218] += x[222];
		x[219] += x[223];

		temp2 = -tempi; tempi = tempr; tempr=temp2;
		x[222] = x[220] - tempr;
		x[223] = x[221] - tempi;

		x[220] += tempr;
		x[221] += tempi;

		tempr = x[216] - x[218];
		tempi = x[217] - x[219];

		x[216] += x[218];
		x[217] += x[219];

		x[218] = tempr;
		x[219] = tempi;

		tempr = x[196] - x[204];
		tempi = x[197] - x[205];

		x[196] += x[204];
		x[197] += x[205];

		x[204] = -tempi;  x[205] = tempr;
		tempr = x[192] - x[200];
		tempi = x[193] - x[201];

		x[192] += x[200];
		x[193] += x[201];

		x[200] = tempr + x[204];
		x[201] = tempi + x[205];

		temp2 = x[200]*wtab[0] - x[201]*wtab[1024]; x[201] = x[200]*wtab[1024] + x[201]*wtab[0]; x[200] = temp2;
		x[204] = tempr - x[204];
		x[205] = tempi - x[205];

		temp2 = x[204]*wtab[0] - x[205]*wtab[1024]; x[205] = x[204]*wtab[1024] + x[205]*wtab[0]; x[204] = temp2;
		tempr = x[198] - x[206];
		tempi = x[199] - x[207];

		x[198] += x[206];
		x[199] += x[207];

		x[206] = -tempi;  x[207] = tempr;
		tempr = x[194] - x[202];
		tempi = x[195] - x[203];

		x[194] += x[202];
		x[195] += x[203];

		x[202] = tempr + x[206];
		x[203] = tempi + x[207];

		temp2 = x[202]*wtab[512] - x[203]*wtab[512]; x[203] = x[202]*wtab[512] + x[203]*wtab[512]; x[202] = temp2;
		x[206] = tempr - x[206];
		x[207] = tempi - x[207];

		temp2 = x[206]*wtab[1536] - x[207]*wtab[512]; x[207] = x[206]*wtab[512] + x[207]*wtab[1536]; x[206] = temp2;
		tempr = x[200] - x[202];
		tempi = x[201] - x[203];

		x[200] += x[202];
		x[201] += x[203];

		x[202] = tempr;
		x[203] = tempi;

		tempr = x[204] - x[206];
		tempi = x[205] - x[207];

		x[204] += x[206];
		x[205] += x[207];

		x[206] = tempr;
		x[207] = tempi;

		tempr = x[192] - x[196];
		tempi = x[193] - x[197];

		x[192] += x[196];
		x[193] += x[197];

		x[196] = tempr;
		x[197] = tempi;

		tempr = x[194] - x[198];
		tempi = x[195] - x[199];

		x[194] += x[198];
		x[195] += x[199];

		temp2 = -tempi; tempi = tempr; tempr=temp2;
		x[198] = x[196] - tempr;
		x[199] = x[197] - tempi;

		x[196] += tempr;
		x[197] += tempi;

		tempr = x[192] - x[194];
		tempi = x[193] - x[195];

		x[192] += x[194];
		x[193] += x[195];

		x[194] = tempr;
		x[195] = tempi;

	    case 64 :
		tempr = x[32] - x[96];
		tempi = x[33] - x[97];

		x[32] += x[96];
		x[33] += x[97];

		x[96] = -tempi;  x[97] = tempr;
		tempr = x[0] - x[64];
		tempi = x[1] - x[65];

		x[0] += x[64];
		x[1] += x[65];

		x[64] = tempr + x[96];
		x[65] = tempi + x[97];

		temp2 = x[64]*wtab[0] - x[65]*wtab[1024]; x[65] = x[64]*wtab[1024] + x[65]*wtab[0]; x[64] = temp2;
		x[96] = tempr - x[96];
		x[97] = tempi - x[97];

		temp2 = x[96]*wtab[0] - x[97]*wtab[1024]; x[97] = x[96]*wtab[1024] + x[97]*wtab[0]; x[96] = temp2;
		tempr = x[34] - x[98];
		tempi = x[35] - x[99];

		x[34] += x[98];
		x[35] += x[99];

		x[98] = -tempi;  x[99] = tempr;
		tempr = x[2] - x[66];
		tempi = x[3] - x[67];

		x[2] += x[66];
		x[3] += x[67];

		x[66] = tempr + x[98];
		x[67] = tempi + x[99];

		temp2 = x[66]*wtab[64] - x[67]*wtab[960]; x[67] = x[66]*wtab[960] + x[67]*wtab[64]; x[66] = temp2;
		x[98] = tempr - x[98];
		x[99] = tempi - x[99];

		temp2 = x[98]*wtab[192] - x[99]*wtab[832]; x[99] = x[98]*wtab[832] + x[99]*wtab[192]; x[98] = temp2;
		tempr = x[36] - x[100];
		tempi = x[37] - x[101];

		x[36] += x[100];
		x[37] += x[101];

		x[100] = -tempi;  x[101] = tempr;
		tempr = x[4] - x[68];
		tempi = x[5] - x[69];

		x[4] += x[68];
		x[5] += x[69];

		x[68] = tempr + x[100];
		x[69] = tempi + x[101];

		temp2 = x[68]*wtab[128] - x[69]*wtab[896]; x[69] = x[68]*wtab[896] + x[69]*wtab[128]; x[68] = temp2;
		x[100] = tempr - x[100];
		x[101] = tempi - x[101];

		temp2 = x[100]*wtab[384] - x[101]*wtab[640]; x[101] = x[100]*wtab[640] + x[101]*wtab[384]; x[100] = temp2;
		tempr = x[38] - x[102];
		tempi = x[39] - x[103];

		x[38] += x[102];
		x[39] += x[103];

		x[102] = -tempi;  x[103] = tempr;
		tempr = x[6] - x[70];
		tempi = x[7] - x[71];

		x[6] += x[70];
		x[7] += x[71];

		x[70] = tempr + x[102];
		x[71] = tempi + x[103];

		temp2 = x[70]*wtab[192] - x[71]*wtab[832]; x[71] = x[70]*wtab[832] + x[71]*wtab[192]; x[70] = temp2;
		x[102] = tempr - x[102];
		x[103] = tempi - x[103];

		temp2 = x[102]*wtab[576] - x[103]*wtab[448]; x[103] = x[102]*wtab[448] + x[103]*wtab[576]; x[102] = temp2;
		tempr = x[40] - x[104];
		tempi = x[41] - x[105];

		x[40] += x[104];
		x[41] += x[105];

		x[104] = -tempi;  x[105] = tempr;
		tempr = x[8] - x[72];
		tempi = x[9] - x[73];

		x[8] += x[72];
		x[9] += x[73];

		x[72] = tempr + x[104];
		x[73] = tempi + x[105];

		temp2 = x[72]*wtab[256] - x[73]*wtab[768]; x[73] = x[72]*wtab[768] + x[73]*wtab[256]; x[72] = temp2;
		x[104] = tempr - x[104];
		x[105] = tempi - x[105];

		temp2 = x[104]*wtab[768] - x[105]*wtab[256]; x[105] = x[104]*wtab[256] + x[105]*wtab[768]; x[104] = temp2;
		tempr = x[42] - x[106];
		tempi = x[43] - x[107];

		x[42] += x[106];
		x[43] += x[107];

		x[106] = -tempi;  x[107] = tempr;
		tempr = x[10] - x[74];
		tempi = x[11] - x[75];

		x[10] += x[74];
		x[11] += x[75];

		x[74] = tempr + x[106];
		x[75] = tempi + x[107];

		temp2 = x[74]*wtab[320] - x[75]*wtab[704]; x[75] = x[74]*wtab[704] + x[75]*wtab[320]; x[74] = temp2;
		x[106] = tempr - x[106];
		x[107] = tempi - x[107];

		temp2 = x[106]*wtab[960] - x[107]*wtab[64]; x[107] = x[106]*wtab[64] + x[107]*wtab[960]; x[106] = temp2;
		tempr = x[44] - x[108];
		tempi = x[45] - x[109];

		x[44] += x[108];
		x[45] += x[109];

		x[108] = -tempi;  x[109] = tempr;
		tempr = x[12] - x[76];
		tempi = x[13] - x[77];

		x[12] += x[76];
		x[13] += x[77];

		x[76] = tempr + x[108];
		x[77] = tempi + x[109];

		temp2 = x[76]*wtab[384] - x[77]*wtab[640]; x[77] = x[76]*wtab[640] + x[77]*wtab[384]; x[76] = temp2;
		x[108] = tempr - x[108];
		x[109] = tempi - x[109];

		temp2 = x[108]*wtab[1152] - x[109]*wtab[128]; x[109] = x[108]*wtab[128] + x[109]*wtab[1152]; x[108] = temp2;
		tempr = x[46] - x[110];
		tempi = x[47] - x[111];

		x[46] += x[110];
		x[47] += x[111];

		x[110] = -tempi;  x[111] = tempr;
		tempr = x[14] - x[78];
		tempi = x[15] - x[79];

		x[14] += x[78];
		x[15] += x[79];

		x[78] = tempr + x[110];
		x[79] = tempi + x[111];

		temp2 = x[78]*wtab[448] - x[79]*wtab[576]; x[79] = x[78]*wtab[576] + x[79]*wtab[448]; x[78] = temp2;
		x[110] = tempr - x[110];
		x[111] = tempi - x[111];

		temp2 = x[110]*wtab[1344] - x[111]*wtab[320]; x[111] = x[110]*wtab[320] + x[111]*wtab[1344]; x[110] = temp2;
		tempr = x[48] - x[112];
		tempi = x[49] - x[113];

		x[48] += x[112];
		x[49] += x[113];

		x[112] = -tempi;  x[113] = tempr;
		tempr = x[16] - x[80];
		tempi = x[17] - x[81];

		x[16] += x[80];
		x[17] += x[81];

		x[80] = tempr + x[112];
		x[81] = tempi + x[113];

		temp2 = x[80]*wtab[512] - x[81]*wtab[512]; x[81] = x[80]*wtab[512] + x[81]*wtab[512]; x[80] = temp2;
		x[112] = tempr - x[112];
		x[113] = tempi - x[113];

		temp2 = x[112]*wtab[1536] - x[113]*wtab[512]; x[113] = x[112]*wtab[512] + x[113]*wtab[1536]; x[112] = temp2;
		tempr = x[50] - x[114];
		tempi = x[51] - x[115];

		x[50] += x[114];
		x[51] += x[115];

		x[114] = -tempi;  x[115] = tempr;
		tempr = x[18] - x[82];
		tempi = x[19] - x[83];

		x[18] += x[82];
		x[19] += x[83];

		x[82] = tempr + x[114];
		x[83] = tempi + x[115];

		temp2 = x[82]*wtab[576] - x[83]*wtab[448]; x[83] = x[82]*wtab[448] + x[83]*wtab[576]; x[82] = temp2;
		x[114] = tempr - x[114];
		x[115] = tempi - x[115];

		temp2 = x[114]*wtab[1728] - x[115]*wtab[704]; x[115] = x[114]*wtab[704] + x[115]*wtab[1728]; x[114] = temp2;
		tempr = x[52] - x[116];
		tempi = x[53] - x[117];

		x[52] += x[116];
		x[53] += x[117];

		x[116] = -tempi;  x[117] = tempr;
		tempr = x[20] - x[84];
		tempi = x[21] - x[85];

		x[20] += x[84];
		x[21] += x[85];

		x[84] = tempr + x[116];
		x[85] = tempi + x[117];

		temp2 = x[84]*wtab[640] - x[85]*wtab[384]; x[85] = x[84]*wtab[384] + x[85]*wtab[640]; x[84] = temp2;
		x[116] = tempr - x[116];
		x[117] = tempi - x[117];

		temp2 = x[116]*wtab[1920] - x[117]*wtab[896]; x[117] = x[116]*wtab[896] + x[117]*wtab[1920]; x[116] = temp2;
		tempr = x[54] - x[118];
		tempi = x[55] - x[119];

		x[54] += x[118];
		x[55] += x[119];

		x[118] = -tempi;  x[119] = tempr;
		tempr = x[22] - x[86];
		tempi = x[23] - x[87];

		x[22] += x[86];
		x[23] += x[87];

		x[86] = tempr + x[118];
		x[87] = tempi + x[119];

		temp2 = x[86]*wtab[704] - x[87]*wtab[320]; x[87] = x[86]*wtab[320] + x[87]*wtab[704]; x[86] = temp2;
		x[118] = tempr - x[118];
		x[119] = tempi - x[119];

		temp2 = x[118]*wtab[1984] - x[119]*wtab[1088]; x[119] = x[118]*wtab[1088] + x[119]*wtab[1984]; x[118] = temp2;
		tempr = x[56] - x[120];
		tempi = x[57] - x[121];

		x[56] += x[120];
		x[57] += x[121];

		x[120] = -tempi;  x[121] = tempr;
		tempr = x[24] - x[88];
		tempi = x[25] - x[89];

		x[24] += x[88];
		x[25] += x[89];

		x[88] = tempr + x[120];
		x[89] = tempi + x[121];

		temp2 = x[88]*wtab[768] - x[89]*wtab[256]; x[89] = x[88]*wtab[256] + x[89]*wtab[768]; x[88] = temp2;
		x[120] = tempr - x[120];
		x[121] = tempi - x[121];

		temp2 = x[120]*wtab[1792] - x[121]*wtab[1280]; x[121] = x[120]*wtab[1280] + x[121]*wtab[1792]; x[120] = temp2;
		tempr = x[58] - x[122];
		tempi = x[59] - x[123];

		x[58] += x[122];
		x[59] += x[123];

		x[122] = -tempi;  x[123] = tempr;
		tempr = x[26] - x[90];
		tempi = x[27] - x[91];

		x[26] += x[90];
		x[27] += x[91];

		x[90] = tempr + x[122];
		x[91] = tempi + x[123];

		temp2 = x[90]*wtab[832] - x[91]*wtab[192]; x[91] = x[90]*wtab[192] + x[91]*wtab[832]; x[90] = temp2;
		x[122] = tempr - x[122];
		x[123] = tempi - x[123];

		temp2 = x[122]*wtab[1600] - x[123]*wtab[1472]; x[123] = x[122]*wtab[1472] + x[123]*wtab[1600]; x[122] = temp2;
		tempr = x[60] - x[124];
		tempi = x[61] - x[125];

		x[60] += x[124];
		x[61] += x[125];

		x[124] = -tempi;  x[125] = tempr;
		tempr = x[28] - x[92];
		tempi = x[29] - x[93];

		x[28] += x[92];
		x[29] += x[93];

		x[92] = tempr + x[124];
		x[93] = tempi + x[125];

		temp2 = x[92]*wtab[896] - x[93]*wtab[128]; x[93] = x[92]*wtab[128] + x[93]*wtab[896]; x[92] = temp2;
		x[124] = tempr - x[124];
		x[125] = tempi - x[125];

		temp2 = x[124]*wtab[1408] - x[125]*wtab[1664]; x[125] = x[124]*wtab[1664] + x[125]*wtab[1408]; x[124] = temp2;
		tempr = x[62] - x[126];
		tempi = x[63] - x[127];

		x[62] += x[126];
		x[63] += x[127];

		x[126] = -tempi;  x[127] = tempr;
		tempr = x[30] - x[94];
		tempi = x[31] - x[95];

		x[30] += x[94];
		x[31] += x[95];

		x[94] = tempr + x[126];
		x[95] = tempi + x[127];

		temp2 = x[94]*wtab[960] - x[95]*wtab[64]; x[95] = x[94]*wtab[64] + x[95]*wtab[960]; x[94] = temp2;
		x[126] = tempr - x[126];
		x[127] = tempi - x[127];

		temp2 = x[126]*wtab[1216] - x[127]*wtab[1856]; x[127] = x[126]*wtab[1856] + x[127]*wtab[1216]; x[126] = temp2;
		tempr = x[72] - x[88];
		tempi = x[73] - x[89];

		x[72] += x[88];
		x[73] += x[89];

		x[88] = -tempi;  x[89] = tempr;
		tempr = x[64] - x[80];
		tempi = x[65] - x[81];

		x[64] += x[80];
		x[65] += x[81];

		x[80] = tempr + x[88];
		x[81] = tempi + x[89];

		temp2 = x[80]*wtab[0] - x[81]*wtab[1024]; x[81] = x[80]*wtab[1024] + x[81]*wtab[0]; x[80] = temp2;
		x[88] = tempr - x[88];
		x[89] = tempi - x[89];

		temp2 = x[88]*wtab[0] - x[89]*wtab[1024]; x[89] = x[88]*wtab[1024] + x[89]*wtab[0]; x[88] = temp2;
		tempr = x[74] - x[90];
		tempi = x[75] - x[91];

		x[74] += x[90];
		x[75] += x[91];

		x[90] = -tempi;  x[91] = tempr;
		tempr = x[66] - x[82];
		tempi = x[67] - x[83];

		x[66] += x[82];
		x[67] += x[83];

		x[82] = tempr + x[90];
		x[83] = tempi + x[91];

		temp2 = x[82]*wtab[256] - x[83]*wtab[768]; x[83] = x[82]*wtab[768] + x[83]*wtab[256]; x[82] = temp2;
		x[90] = tempr - x[90];
		x[91] = tempi - x[91];

		temp2 = x[90]*wtab[768] - x[91]*wtab[256]; x[91] = x[90]*wtab[256] + x[91]*wtab[768]; x[90] = temp2;
		tempr = x[76] - x[92];
		tempi = x[77] - x[93];

		x[76] += x[92];
		x[77] += x[93];

		x[92] = -tempi;  x[93] = tempr;
		tempr = x[68] - x[84];
		tempi = x[69] - x[85];

		x[68] += x[84];
		x[69] += x[85];

		x[84] = tempr + x[92];
		x[85] = tempi + x[93];

		temp2 = x[84]*wtab[512] - x[85]*wtab[512]; x[85] = x[84]*wtab[512] + x[85]*wtab[512]; x[84] = temp2;
		x[92] = tempr - x[92];
		x[93] = tempi - x[93];

		temp2 = x[92]*wtab[1536] - x[93]*wtab[512]; x[93] = x[92]*wtab[512] + x[93]*wtab[1536]; x[92] = temp2;
		tempr = x[78] - x[94];
		tempi = x[79] - x[95];

		x[78] += x[94];
		x[79] += x[95];

		x[94] = -tempi;  x[95] = tempr;
		tempr = x[70] - x[86];
		tempi = x[71] - x[87];

		x[70] += x[86];
		x[71] += x[87];

		x[86] = tempr + x[94];
		x[87] = tempi + x[95];

		temp2 = x[86]*wtab[768] - x[87]*wtab[256]; x[87] = x[86]*wtab[256] + x[87]*wtab[768]; x[86] = temp2;
		x[94] = tempr - x[94];
		x[95] = tempi - x[95];

		temp2 = x[94]*wtab[1792] - x[95]*wtab[1280]; x[95] = x[94]*wtab[1280] + x[95]*wtab[1792]; x[94] = temp2;
		tempr = x[80] - x[84];
		tempi = x[81] - x[85];

		x[80] += x[84];
		x[81] += x[85];

		x[84] = tempr;
		x[85] = tempi;

		tempr = x[82] - x[86];
		tempi = x[83] - x[87];

		x[82] += x[86];
		x[83] += x[87];

		temp2 = -tempi; tempi = tempr; tempr=temp2;
		x[86] = x[84] - tempr;
		x[87] = x[85] - tempi;

		x[84] += tempr;
		x[85] += tempi;

		tempr = x[80] - x[82];
		tempi = x[81] - x[83];

		x[80] += x[82];
		x[81] += x[83];

		x[82] = tempr;
		x[83] = tempi;

		tempr = x[88] - x[92];
		tempi = x[89] - x[93];

		x[88] += x[92];
		x[89] += x[93];

		x[92] = tempr;
		x[93] = tempi;

		tempr = x[90] - x[94];
		tempi = x[91] - x[95];

		x[90] += x[94];
		x[91] += x[95];

		temp2 = -tempi; tempi = tempr; tempr=temp2;
		x[94] = x[92] - tempr;
		x[95] = x[93] - tempi;

		x[92] += tempr;
		x[93] += tempi;

		tempr = x[88] - x[90];
		tempi = x[89] - x[91];

		x[88] += x[90];
		x[89] += x[91];

		x[90] = tempr;
		x[91] = tempi;

		tempr = x[68] - x[76];
		tempi = x[69] - x[77];

		x[68] += x[76];
		x[69] += x[77];

		x[76] = -tempi;  x[77] = tempr;
		tempr = x[64] - x[72];
		tempi = x[65] - x[73];

		x[64] += x[72];
		x[65] += x[73];

		x[72] = tempr + x[76];
		x[73] = tempi + x[77];

		temp2 = x[72]*wtab[0] - x[73]*wtab[1024]; x[73] = x[72]*wtab[1024] + x[73]*wtab[0]; x[72] = temp2;
		x[76] = tempr - x[76];
		x[77] = tempi - x[77];

		temp2 = x[76]*wtab[0] - x[77]*wtab[1024]; x[77] = x[76]*wtab[1024] + x[77]*wtab[0]; x[76] = temp2;
		tempr = x[70] - x[78];
		tempi = x[71] - x[79];

		x[70] += x[78];
		x[71] += x[79];

		x[78] = -tempi;  x[79] = tempr;
		tempr = x[66] - x[74];
		tempi = x[67] - x[75];

		x[66] += x[74];
		x[67] += x[75];

		x[74] = tempr + x[78];
		x[75] = tempi + x[79];

		temp2 = x[74]*wtab[512] - x[75]*wtab[512]; x[75] = x[74]*wtab[512] + x[75]*wtab[512]; x[74] = temp2;
		x[78] = tempr - x[78];
		x[79] = tempi - x[79];

		temp2 = x[78]*wtab[1536] - x[79]*wtab[512]; x[79] = x[78]*wtab[512] + x[79]*wtab[1536]; x[78] = temp2;
		tempr = x[72] - x[74];
		tempi = x[73] - x[75];

		x[72] += x[74];
		x[73] += x[75];

		x[74] = tempr;
		x[75] = tempi;

		tempr = x[76] - x[78];
		tempi = x[77] - x[79];

		x[76] += x[78];
		x[77] += x[79];

		x[78] = tempr;
		x[79] = tempi;

		tempr = x[64] - x[68];
		tempi = x[65] - x[69];

		x[64] += x[68];
		x[65] += x[69];

		x[68] = tempr;
		x[69] = tempi;

		tempr = x[66] - x[70];
		tempi = x[67] - x[71];

		x[66] += x[70];
		x[67] += x[71];

		temp2 = -tempi; tempi = tempr; tempr=temp2;
		x[70] = x[68] - tempr;
		x[71] = x[69] - tempi;

		x[68] += tempr;
		x[69] += tempi;

		tempr = x[64] - x[66];
		tempi = x[65] - x[67];

		x[64] += x[66];
		x[65] += x[67];

		x[66] = tempr;
		x[67] = tempi;

		tempr = x[104] - x[120];
		tempi = x[105] - x[121];

		x[104] += x[120];
		x[105] += x[121];

		x[120] = -tempi;  x[121] = tempr;
		tempr = x[96] - x[112];
		tempi = x[97] - x[113];

		x[96] += x[112];
		x[97] += x[113];

		x[112] = tempr + x[120];
		x[113] = tempi + x[121];

		temp2 = x[112]*wtab[0] - x[113]*wtab[1024]; x[113] = x[112]*wtab[1024] + x[113]*wtab[0]; x[112] = temp2;
		x[120] = tempr - x[120];
		x[121] = tempi - x[121];

		temp2 = x[120]*wtab[0] - x[121]*wtab[1024]; x[121] = x[120]*wtab[1024] + x[121]*wtab[0]; x[120] = temp2;
		tempr = x[106] - x[122];
		tempi = x[107] - x[123];

		x[106] += x[122];
		x[107] += x[123];

		x[122] = -tempi;  x[123] = tempr;
		tempr = x[98] - x[114];
		tempi = x[99] - x[115];

		x[98] += x[114];
		x[99] += x[115];

		x[114] = tempr + x[122];
		x[115] = tempi + x[123];

		temp2 = x[114]*wtab[256] - x[115]*wtab[768]; x[115] = x[114]*wtab[768] + x[115]*wtab[256]; x[114] = temp2;
		x[122] = tempr - x[122];
		x[123] = tempi - x[123];

		temp2 = x[122]*wtab[768] - x[123]*wtab[256]; x[123] = x[122]*wtab[256] + x[123]*wtab[768]; x[122] = temp2;
		tempr = x[108] - x[124];
		tempi = x[109] - x[125];

		x[108] += x[124];
		x[109] += x[125];

		x[124] = -tempi;  x[125] = tempr;
		tempr = x[100] - x[116];
		tempi = x[101] - x[117];

		x[100] += x[116];
		x[101] += x[117];

		x[116] = tempr + x[124];
		x[117] = tempi + x[125];

		temp2 = x[116]*wtab[512] - x[117]*wtab[512]; x[117] = x[116]*wtab[512] + x[117]*wtab[512]; x[116] = temp2;
		x[124] = tempr - x[124];
		x[125] = tempi - x[125];

		temp2 = x[124]*wtab[1536] - x[125]*wtab[512]; x[125] = x[124]*wtab[512] + x[125]*wtab[1536]; x[124] = temp2;
		tempr = x[110] - x[126];
		tempi = x[111] - x[127];

		x[110] += x[126];
		x[111] += x[127];

		x[126] = -tempi;  x[127] = tempr;
		tempr = x[102] - x[118];
		tempi = x[103] - x[119];

		x[102] += x[118];
		x[103] += x[119];

		x[118] = tempr + x[126];
		x[119] = tempi + x[127];

		temp2 = x[118]*wtab[768] - x[119]*wtab[256]; x[119] = x[118]*wtab[256] + x[119]*wtab[768]; x[118] = temp2;
		x[126] = tempr - x[126];
		x[127] = tempi - x[127];

		temp2 = x[126]*wtab[1792] - x[127]*wtab[1280]; x[127] = x[126]*wtab[1280] + x[127]*wtab[1792]; x[126] = temp2;
		tempr = x[112] - x[116];
		tempi = x[113] - x[117];

		x[112] += x[116];
		x[113] += x[117];

		x[116] = tempr;
		x[117] = tempi;

		tempr = x[114] - x[118];
		tempi = x[115] - x[119];

		x[114] += x[118];
		x[115] += x[119];

		temp2 = -tempi; tempi = tempr; tempr=temp2;
		x[118] = x[116] - tempr;
		x[119] = x[117] - tempi;

		x[116] += tempr;
		x[117] += tempi;

		tempr = x[112] - x[114];
		tempi = x[113] - x[115];

		x[112] += x[114];
		x[113] += x[115];

		x[114] = tempr;
		x[115] = tempi;

		tempr = x[120] - x[124];
		tempi = x[121] - x[125];

		x[120] += x[124];
		x[121] += x[125];

		x[124] = tempr;
		x[125] = tempi;

		tempr = x[122] - x[126];
		tempi = x[123] - x[127];

		x[122] += x[126];
		x[123] += x[127];

		temp2 = -tempi; tempi = tempr; tempr=temp2;
		x[126] = x[124] - tempr;
		x[127] = x[125] - tempi;

		x[124] += tempr;
		x[125] += tempi;

		tempr = x[120] - x[122];
		tempi = x[121] - x[123];

		x[120] += x[122];
		x[121] += x[123];

		x[122] = tempr;
		x[123] = tempi;

		tempr = x[100] - x[108];
		tempi = x[101] - x[109];

		x[100] += x[108];
		x[101] += x[109];

		x[108] = -tempi;  x[109] = tempr;
		tempr = x[96] - x[104];
		tempi = x[97] - x[105];

		x[96] += x[104];
		x[97] += x[105];

		x[104] = tempr + x[108];
		x[105] = tempi + x[109];

		temp2 = x[104]*wtab[0] - x[105]*wtab[1024]; x[105] = x[104]*wtab[1024] + x[105]*wtab[0]; x[104] = temp2;
		x[108] = tempr - x[108];
		x[109] = tempi - x[109];

		temp2 = x[108]*wtab[0] - x[109]*wtab[1024]; x[109] = x[108]*wtab[1024] + x[109]*wtab[0]; x[108] = temp2;
		tempr = x[102] - x[110];
		tempi = x[103] - x[111];

		x[102] += x[110];
		x[103] += x[111];

		x[110] = -tempi;  x[111] = tempr;
		tempr = x[98] - x[106];
		tempi = x[99] - x[107];

		x[98] += x[106];
		x[99] += x[107];

		x[106] = tempr + x[110];
		x[107] = tempi + x[111];

		temp2 = x[106]*wtab[512] - x[107]*wtab[512]; x[107] = x[106]*wtab[512] + x[107]*wtab[512]; x[106] = temp2;
		x[110] = tempr - x[110];
		x[111] = tempi - x[111];

		temp2 = x[110]*wtab[1536] - x[111]*wtab[512]; x[111] = x[110]*wtab[512] + x[111]*wtab[1536]; x[110] = temp2;
		tempr = x[104] - x[106];
		tempi = x[105] - x[107];

		x[104] += x[106];
		x[105] += x[107];

		x[106] = tempr;
		x[107] = tempi;

		tempr = x[108] - x[110];
		tempi = x[109] - x[111];

		x[108] += x[110];
		x[109] += x[111];

		x[110] = tempr;
		x[111] = tempi;

		tempr = x[96] - x[100];
		tempi = x[97] - x[101];

		x[96] += x[100];
		x[97] += x[101];

		x[100] = tempr;
		x[101] = tempi;

		tempr = x[98] - x[102];
		tempi = x[99] - x[103];

		x[98] += x[102];
		x[99] += x[103];

		temp2 = -tempi; tempi = tempr; tempr=temp2;
		x[102] = x[100] - tempr;
		x[103] = x[101] - tempi;

		x[100] += tempr;
		x[101] += tempi;

		tempr = x[96] - x[98];
		tempi = x[97] - x[99];

		x[96] += x[98];
		x[97] += x[99];

		x[98] = tempr;
		x[99] = tempi;

	    case 32 :
		tempr = x[16] - x[48];
		tempi = x[17] - x[49];

		x[16] += x[48];
		x[17] += x[49];

		x[48] = -tempi;  x[49] = tempr;
		tempr = x[0] - x[32];
		tempi = x[1] - x[33];

		x[0] += x[32];
		x[1] += x[33];

		x[32] = tempr + x[48];
		x[33] = tempi + x[49];

		temp2 = x[32]*wtab[0] - x[33]*wtab[1024]; x[33] = x[32]*wtab[1024] + x[33]*wtab[0]; x[32] = temp2;
		x[48] = tempr - x[48];
		x[49] = tempi - x[49];

		temp2 = x[48]*wtab[0] - x[49]*wtab[1024]; x[49] = x[48]*wtab[1024] + x[49]*wtab[0]; x[48] = temp2;
		tempr = x[18] - x[50];
		tempi = x[19] - x[51];

		x[18] += x[50];
		x[19] += x[51];

		x[50] = -tempi;  x[51] = tempr;
		tempr = x[2] - x[34];
		tempi = x[3] - x[35];

		x[2] += x[34];
		x[3] += x[35];

		x[34] = tempr + x[50];
		x[35] = tempi + x[51];

		temp2 = x[34]*wtab[128] - x[35]*wtab[896]; x[35] = x[34]*wtab[896] + x[35]*wtab[128]; x[34] = temp2;
		x[50] = tempr - x[50];
		x[51] = tempi - x[51];

		temp2 = x[50]*wtab[384] - x[51]*wtab[640]; x[51] = x[50]*wtab[640] + x[51]*wtab[384]; x[50] = temp2;
		tempr = x[20] - x[52];
		tempi = x[21] - x[53];

		x[20] += x[52];
		x[21] += x[53];

		x[52] = -tempi;  x[53] = tempr;
		tempr = x[4] - x[36];
		tempi = x[5] - x[37];

		x[4] += x[36];
		x[5] += x[37];

		x[36] = tempr + x[52];
		x[37] = tempi + x[53];

		temp2 = x[36]*wtab[256] - x[37]*wtab[768]; x[37] = x[36]*wtab[768] + x[37]*wtab[256]; x[36] = temp2;
		x[52] = tempr - x[52];
		x[53] = tempi - x[53];

		temp2 = x[52]*wtab[768] - x[53]*wtab[256]; x[53] = x[52]*wtab[256] + x[53]*wtab[768]; x[52] = temp2;
		tempr = x[22] - x[54];
		tempi = x[23] - x[55];

		x[22] += x[54];
		x[23] += x[55];

		x[54] = -tempi;  x[55] = tempr;
		tempr = x[6] - x[38];
		tempi = x[7] - x[39];

		x[6] += x[38];
		x[7] += x[39];

		x[38] = tempr + x[54];
		x[39] = tempi + x[55];

		temp2 = x[38]*wtab[384] - x[39]*wtab[640]; x[39] = x[38]*wtab[640] + x[39]*wtab[384]; x[38] = temp2;
		x[54] = tempr - x[54];
		x[55] = tempi - x[55];

		temp2 = x[54]*wtab[1152] - x[55]*wtab[128]; x[55] = x[54]*wtab[128] + x[55]*wtab[1152]; x[54] = temp2;
		tempr = x[24] - x[56];
		tempi = x[25] - x[57];

		x[24] += x[56];
		x[25] += x[57];

		x[56] = -tempi;  x[57] = tempr;
		tempr = x[8] - x[40];
		tempi = x[9] - x[41];

		x[8] += x[40];
		x[9] += x[41];

		x[40] = tempr + x[56];
		x[41] = tempi + x[57];

		temp2 = x[40]*wtab[512] - x[41]*wtab[512]; x[41] = x[40]*wtab[512] + x[41]*wtab[512]; x[40] = temp2;
		x[56] = tempr - x[56];
		x[57] = tempi - x[57];

		temp2 = x[56]*wtab[1536] - x[57]*wtab[512]; x[57] = x[56]*wtab[512] + x[57]*wtab[1536]; x[56] = temp2;
		tempr = x[26] - x[58];
		tempi = x[27] - x[59];

		x[26] += x[58];
		x[27] += x[59];

		x[58] = -tempi;  x[59] = tempr;
		tempr = x[10] - x[42];
		tempi = x[11] - x[43];

		x[10] += x[42];
		x[11] += x[43];

		x[42] = tempr + x[58];
		x[43] = tempi + x[59];

		temp2 = x[42]*wtab[640] - x[43]*wtab[384]; x[43] = x[42]*wtab[384] + x[43]*wtab[640]; x[42] = temp2;
		x[58] = tempr - x[58];
		x[59] = tempi - x[59];

		temp2 = x[58]*wtab[1920] - x[59]*wtab[896]; x[59] = x[58]*wtab[896] + x[59]*wtab[1920]; x[58] = temp2;
		tempr = x[28] - x[60];
		tempi = x[29] - x[61];

		x[28] += x[60];
		x[29] += x[61];

		x[60] = -tempi;  x[61] = tempr;
		tempr = x[12] - x[44];
		tempi = x[13] - x[45];

		x[12] += x[44];
		x[13] += x[45];

		x[44] = tempr + x[60];
		x[45] = tempi + x[61];

		temp2 = x[44]*wtab[768] - x[45]*wtab[256]; x[45] = x[44]*wtab[256] + x[45]*wtab[768]; x[44] = temp2;
		x[60] = tempr - x[60];
		x[61] = tempi - x[61];

		temp2 = x[60]*wtab[1792] - x[61]*wtab[1280]; x[61] = x[60]*wtab[1280] + x[61]*wtab[1792]; x[60] = temp2;
		tempr = x[30] - x[62];
		tempi = x[31] - x[63];

		x[30] += x[62];
		x[31] += x[63];

		x[62] = -tempi;  x[63] = tempr;
		tempr = x[14] - x[46];
		tempi = x[15] - x[47];

		x[14] += x[46];
		x[15] += x[47];

		x[46] = tempr + x[62];
		x[47] = tempi + x[63];

		temp2 = x[46]*wtab[896] - x[47]*wtab[128]; x[47] = x[46]*wtab[128] + x[47]*wtab[896]; x[46] = temp2;
		x[62] = tempr - x[62];
		x[63] = tempi - x[63];

		temp2 = x[62]*wtab[1408] - x[63]*wtab[1664]; x[63] = x[62]*wtab[1664] + x[63]*wtab[1408]; x[62] = temp2;
		tempr = x[36] - x[44];
		tempi = x[37] - x[45];

		x[36] += x[44];
		x[37] += x[45];

		x[44] = -tempi;  x[45] = tempr;
		tempr = x[32] - x[40];
		tempi = x[33] - x[41];

		x[32] += x[40];
		x[33] += x[41];

		x[40] = tempr + x[44];
		x[41] = tempi + x[45];

		temp2 = x[40]*wtab[0] - x[41]*wtab[1024]; x[41] = x[40]*wtab[1024] + x[41]*wtab[0]; x[40] = temp2;
		x[44] = tempr - x[44];
		x[45] = tempi - x[45];

		temp2 = x[44]*wtab[0] - x[45]*wtab[1024]; x[45] = x[44]*wtab[1024] + x[45]*wtab[0]; x[44] = temp2;
		tempr = x[38] - x[46];
		tempi = x[39] - x[47];

		x[38] += x[46];
		x[39] += x[47];

		x[46] = -tempi;  x[47] = tempr;
		tempr = x[34] - x[42];
		tempi = x[35] - x[43];

		x[34] += x[42];
		x[35] += x[43];

		x[42] = tempr + x[46];
		x[43] = tempi + x[47];

		temp2 = x[42]*wtab[512] - x[43]*wtab[512]; x[43] = x[42]*wtab[512] + x[43]*wtab[512]; x[42] = temp2;
		x[46] = tempr - x[46];
		x[47] = tempi - x[47];

		temp2 = x[46]*wtab[1536] - x[47]*wtab[512]; x[47] = x[46]*wtab[512] + x[47]*wtab[1536]; x[46] = temp2;
		tempr = x[40] - x[42];
		tempi = x[41] - x[43];

		x[40] += x[42];
		x[41] += x[43];

		x[42] = tempr;
		x[43] = tempi;

		tempr = x[44] - x[46];
		tempi = x[45] - x[47];

		x[44] += x[46];
		x[45] += x[47];

		x[46] = tempr;
		x[47] = tempi;

		tempr = x[32] - x[36];
		tempi = x[33] - x[37];

		x[32] += x[36];
		x[33] += x[37];

		x[36] = tempr;
		x[37] = tempi;

		tempr = x[34] - x[38];
		tempi = x[35] - x[39];

		x[34] += x[38];
		x[35] += x[39];

		temp2 = -tempi; tempi = tempr; tempr=temp2;
		x[38] = x[36] - tempr;
		x[39] = x[37] - tempi;

		x[36] += tempr;
		x[37] += tempi;

		tempr = x[32] - x[34];
		tempi = x[33] - x[35];

		x[32] += x[34];
		x[33] += x[35];

		x[34] = tempr;
		x[35] = tempi;

		tempr = x[52] - x[60];
		tempi = x[53] - x[61];

		x[52] += x[60];
		x[53] += x[61];

		x[60] = -tempi;  x[61] = tempr;
		tempr = x[48] - x[56];
		tempi = x[49] - x[57];

		x[48] += x[56];
		x[49] += x[57];

		x[56] = tempr + x[60];
		x[57] = tempi + x[61];

		temp2 = x[56]*wtab[0] - x[57]*wtab[1024]; x[57] = x[56]*wtab[1024] + x[57]*wtab[0]; x[56] = temp2;
		x[60] = tempr - x[60];
		x[61] = tempi - x[61];

		temp2 = x[60]*wtab[0] - x[61]*wtab[1024]; x[61] = x[60]*wtab[1024] + x[61]*wtab[0]; x[60] = temp2;
		tempr = x[54] - x[62];
		tempi = x[55] - x[63];

		x[54] += x[62];
		x[55] += x[63];

		x[62] = -tempi;  x[63] = tempr;
		tempr = x[50] - x[58];
		tempi = x[51] - x[59];

		x[50] += x[58];
		x[51] += x[59];

		x[58] = tempr + x[62];
		x[59] = tempi + x[63];

		temp2 = x[58]*wtab[512] - x[59]*wtab[512]; x[59] = x[58]*wtab[512] + x[59]*wtab[512]; x[58] = temp2;
		x[62] = tempr - x[62];
		x[63] = tempi - x[63];

		temp2 = x[62]*wtab[1536] - x[63]*wtab[512]; x[63] = x[62]*wtab[512] + x[63]*wtab[1536]; x[62] = temp2;
		tempr = x[56] - x[58];
		tempi = x[57] - x[59];

		x[56] += x[58];
		x[57] += x[59];

		x[58] = tempr;
		x[59] = tempi;

		tempr = x[60] - x[62];
		tempi = x[61] - x[63];

		x[60] += x[62];
		x[61] += x[63];

		x[62] = tempr;
		x[63] = tempi;

		tempr = x[48] - x[52];
		tempi = x[49] - x[53];

		x[48] += x[52];
		x[49] += x[53];

		x[52] = tempr;
		x[53] = tempi;

		tempr = x[50] - x[54];
		tempi = x[51] - x[55];

		x[50] += x[54];
		x[51] += x[55];

		temp2 = -tempi; tempi = tempr; tempr=temp2;
		x[54] = x[52] - tempr;
		x[55] = x[53] - tempi;

		x[52] += tempr;
		x[53] += tempi;

		tempr = x[48] - x[50];
		tempi = x[49] - x[51];

		x[48] += x[50];
		x[49] += x[51];

		x[50] = tempr;
		x[51] = tempi;

	    case 16 :
		tempr = x[8] - x[24];
		tempi = x[9] - x[25];

		x[8] += x[24];
		x[9] += x[25];

		x[24] = -tempi;  x[25] = tempr;
		tempr = x[0] - x[16];
		tempi = x[1] - x[17];

		x[0] += x[16];
		x[1] += x[17];

		x[16] = tempr + x[24];
		x[17] = tempi + x[25];

		temp2 = x[16]*wtab[0] - x[17]*wtab[1024]; x[17] = x[16]*wtab[1024] + x[17]*wtab[0]; x[16] = temp2;
		x[24] = tempr - x[24];
		x[25] = tempi - x[25];

		temp2 = x[24]*wtab[0] - x[25]*wtab[1024]; x[25] = x[24]*wtab[1024] + x[25]*wtab[0]; x[24] = temp2;
		tempr = x[10] - x[26];
		tempi = x[11] - x[27];

		x[10] += x[26];
		x[11] += x[27];

		x[26] = -tempi;  x[27] = tempr;
		tempr = x[2] - x[18];
		tempi = x[3] - x[19];

		x[2] += x[18];
		x[3] += x[19];

		x[18] = tempr + x[26];
		x[19] = tempi + x[27];

		temp2 = x[18]*wtab[256] - x[19]*wtab[768]; x[19] = x[18]*wtab[768] + x[19]*wtab[256]; x[18] = temp2;
		x[26] = tempr - x[26];
		x[27] = tempi - x[27];

		temp2 = x[26]*wtab[768] - x[27]*wtab[256]; x[27] = x[26]*wtab[256] + x[27]*wtab[768]; x[26] = temp2;
		tempr = x[12] - x[28];
		tempi = x[13] - x[29];

		x[12] += x[28];
		x[13] += x[29];

		x[28] = -tempi;  x[29] = tempr;
		tempr = x[4] - x[20];
		tempi = x[5] - x[21];

		x[4] += x[20];
		x[5] += x[21];

		x[20] = tempr + x[28];
		x[21] = tempi + x[29];

		temp2 = x[20]*wtab[512] - x[21]*wtab[512]; x[21] = x[20]*wtab[512] + x[21]*wtab[512]; x[20] = temp2;
		x[28] = tempr - x[28];
		x[29] = tempi - x[29];

		temp2 = x[28]*wtab[1536] - x[29]*wtab[512]; x[29] = x[28]*wtab[512] + x[29]*wtab[1536]; x[28] = temp2;
		tempr = x[14] - x[30];
		tempi = x[15] - x[31];

		x[14] += x[30];
		x[15] += x[31];

		x[30] = -tempi;  x[31] = tempr;
		tempr = x[6] - x[22];
		tempi = x[7] - x[23];

		x[6] += x[22];
		x[7] += x[23];

		x[22] = tempr + x[30];
		x[23] = tempi + x[31];

		temp2 = x[22]*wtab[768] - x[23]*wtab[256]; x[23] = x[22]*wtab[256] + x[23]*wtab[768]; x[22] = temp2;
		x[30] = tempr - x[30];
		x[31] = tempi - x[31];

		temp2 = x[30]*wtab[1792] - x[31]*wtab[1280]; x[31] = x[30]*wtab[1280] + x[31]*wtab[1792]; x[30] = temp2;
		tempr = x[16] - x[20];
		tempi = x[17] - x[21];

		x[16] += x[20];
		x[17] += x[21];

		x[20] = tempr;
		x[21] = tempi;

		tempr = x[18] - x[22];
		tempi = x[19] - x[23];

		x[18] += x[22];
		x[19] += x[23];

		temp2 = -tempi; tempi = tempr; tempr=temp2;
		x[22] = x[20] - tempr;
		x[23] = x[21] - tempi;

		x[20] += tempr;
		x[21] += tempi;

		tempr = x[16] - x[18];
		tempi = x[17] - x[19];

		x[16] += x[18];
		x[17] += x[19];

		x[18] = tempr;
		x[19] = tempi;

		tempr = x[24] - x[28];
		tempi = x[25] - x[29];

		x[24] += x[28];
		x[25] += x[29];

		x[28] = tempr;
		x[29] = tempi;

		tempr = x[26] - x[30];
		tempi = x[27] - x[31];

		x[26] += x[30];
		x[27] += x[31];

		temp2 = -tempi; tempi = tempr; tempr=temp2;
		x[30] = x[28] - tempr;
		x[31] = x[29] - tempi;

		x[28] += tempr;
		x[29] += tempi;

		tempr = x[24] - x[26];
		tempi = x[25] - x[27];

		x[24] += x[26];
		x[25] += x[27];

		x[26] = tempr;
		x[27] = tempi;

	    case 8 :
		tempr = x[4] - x[12];
		tempi = x[5] - x[13];

		x[4] += x[12];
		x[5] += x[13];

		x[12] = -tempi;  x[13] = tempr;
		tempr = x[0] - x[8];
		tempi = x[1] - x[9];

		x[0] += x[8];
		x[1] += x[9];

		x[8] = tempr + x[12];
		x[9] = tempi + x[13];

		temp2 = x[8]*wtab[0] - x[9]*wtab[1024]; x[9] = x[8]*wtab[1024] + x[9]*wtab[0]; x[8] = temp2;
		x[12] = tempr - x[12];
		x[13] = tempi - x[13];

		temp2 = x[12]*wtab[0] - x[13]*wtab[1024]; x[13] = x[12]*wtab[1024] + x[13]*wtab[0]; x[12] = temp2;
		tempr = x[6] - x[14];
		tempi = x[7] - x[15];

		x[6] += x[14];
		x[7] += x[15];

		x[14] = -tempi;  x[15] = tempr;
		tempr = x[2] - x[10];
		tempi = x[3] - x[11];

		x[2] += x[10];
		x[3] += x[11];

		x[10] = tempr + x[14];
		x[11] = tempi + x[15];

		temp2 = x[10]*wtab[512] - x[11]*wtab[512]; x[11] = x[10]*wtab[512] + x[11]*wtab[512]; x[10] = temp2;
		x[14] = tempr - x[14];
		x[15] = tempi - x[15];

		temp2 = x[14]*wtab[1536] - x[15]*wtab[512]; x[15] = x[14]*wtab[512] + x[15]*wtab[1536]; x[14] = temp2;
		tempr = x[8] - x[10];
		tempi = x[9] - x[11];

		x[8] += x[10];
		x[9] += x[11];

		x[10] = tempr;
		x[11] = tempi;

		tempr = x[12] - x[14];
		tempi = x[13] - x[15];

		x[12] += x[14];
		x[13] += x[15];

		x[14] = tempr;
		x[15] = tempi;

	    case 4 :
		tempr = x[0] - x[4];
		tempi = x[1] - x[5];

		x[0] += x[4];
		x[1] += x[5];

		x[4] = tempr;
		x[5] = tempi;

		tempr = x[2] - x[6];
		tempi = x[3] - x[7];

		x[2] += x[6];
		x[3] += x[7];

		temp2 = -tempi; tempi = tempr; tempr=temp2;
		x[6] = x[4] - tempr;
		x[7] = x[5] - tempi;

		x[4] += tempr;
		x[5] += tempi;

		tempr = x[0] - x[2];
		tempi = x[1] - x[3];

		x[0] += x[2];
		x[1] += x[3];

		x[2] = tempr;
		x[3] = tempi;

		return;
	    case 2 :
		tempr = x[0] - x[2];
		tempi = x[1] - x[3];

		x[0] += x[2];
		x[1] += x[3];

		x[2] = tempr;
		x[3] = tempi;
	}
}

