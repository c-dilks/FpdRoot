#include "FpdMap.h"
ClassImp(FpdMap)
enum setname {ENorth,ESouth,ETop,EBottom,EpsNorth,EpsSouth,
		WNorth,WSouth,WTop,WBottom,WpsNorth,WpsSouth};

Int_t fpdenmap_1[49]=    {39, 38, 37, 36, 35, 34, 33, 
			  7, 6, 5, 23, 22, 21, 55, 
			  4, 3, 2, 20, 19, 18, 54, 
			  1, 0, 15, 17, 16, 31, 53, 
			  14, 13, 12, 30, 29, 28, 52, 
			  11, 10, 9, 27, 26, 25, 51, 
			  32, 47, 46, 45, 44, 43, 42};

Int_t  fpdesmap_1[49]={103, 101, 100, 99, 98, 97, 96, 
			   71, 70, 69, 87, 86, 85, 48, 
			   68, 67, 66, 84, 83, 82, 63, 
			   65, 64, 79, 81, 80, 95, 61, 
			   78, 77, 76, 94, 93, 92, 60, 
			   75, 74, 73, 91, 90, 89, 59, 
			   111, 110, 109, 108, 107, 106, 105};

Int_t fpdetmap_1[25]=            {135, 134, 133, 132, 131, 
				  130, 129, 128, 143, 142, 
				  119, 118, 117, 116, 115, 
				  114, 113, 112, 127, 126, 
				  125, 124, 123, 122, 121};

Int_t fpdebmap_1[25]=            {151, 150, 149, 148, 147, 
                                  146, 145, 144, 159, 158, 
                                  157, 156, 155, 154, 153, 
                                  167, 166, 165, 164, 163, 
                                  162, 161, 160, 175, 174};
		
Int_t fpdwsmap_0[49]=           {7, 6, 25, 22, 39, 38, 48, 
  				 5, 4, 21, 20, 37, 36, 54, 
				 3, 2, 19, 18, 35, 34, 53, 
				 1, 0, 17, 16, 33, 32, 52, 
				 15, 14, 31, 30, 47, 46, 51, 
				 13, 12, 29, 28, 45, 44, 50, 
				 9, 10, 27, 26, 43, 42, 49};
                    
Int_t fpdwsmap_1[49]=           { 28, 29, 30, 31, 40, 41, 42,
                                  12, 13, 10, 11, 56, 33, 43,
                                  14, 15,  8, 9, 34, 35, 44,
                                  6,  7,  3,  0, 36, 37, 45,
                                  4,  5,  1,  2, 38, 39, 46,
		                  24, 25, 16, 17, 22, 23, 47,
		                  26, 27, 18, 19, 57, 49, 50};

Int_t fpdwbmap_1[25]=           {77, 70, 69, 68, 67, 
				 66, 65, 64, 79, 78, 
				 87, 86, 85, 84, 83, 
				 82, 81, 80, 95, 94, 
				 93, 76, 91, 90, 89};

Int_t fpd_prswn_map[7]=          {-1,-1,-1,-1,-1,-1,-1};
Int_t  fpd_prsws_map[7]=         {63,62,61 ,60,59,58,-1};
Int_t  fpd_prsws_map_dau[7]=     {62,63,64,65,66,67,68};
Int_t fpd_prsen_map[7]=          {50, 49, 141, 140, 139, 138, 137};
Int_t fpd_prses_map[7]=          {58, 57, 173, 172, 171, 170, 169};

/* 2005 entries */
    
Int_t  fpd_wn_map2005[49]=
     { 39, 38, 37, 36, 35, 34, 33,
       7,  6,  5, 23, 22, 21, 55,
       4,  3,  2, 20, 19, 18, 54,
       1,  0, 15, 17, 16, 31, 53,
      14, 13, 12, 30, 29, 28, 52,
      41, 10,  9, 27, 26, 25, 51,
	 32, 47, 46, 45, 44, 43, 42};

Int_t  fpd_wt_map2005[25]=
     {141, 134, 133, 132, 131, 130, 129, 140, 143, 142,
     119, 118, 117, 116, 115, 114, 113, 112,
	127, 120, 125, 124, 123, 122, 121};
  


Int_t fpd_wn_map2005pp[49]=
     { 39, 38, 37, 36, 35, 34, 33,
       7,  6,  5, 23, 22, 21, 55,
       4,  3,  2, 20, 19, 18, 54,
       1,  0, 15, 17, 16, 31, 53,
      14, 13, 12, 30, 29, 28, 52,
      41, 10,  9, 27, 26, 25, 51,
	 32, 47, 46, 45, 44, 43, -1};

Int_t fpd_ws_map2005pp[49]=
     {103,102,101,100,  99,  98, 88,
      71, 70, 69, 87,  86,  85, 48,
      68, 67, 66, 84,  83,  82, 63,
      65, 64, 97, 81,  80,  95, 62,
      78, 77 ,76, 94,  93,  92, 61,
     106, 74, 73, 91,  90,  89, 60,
      96,111,110,109,108, 107, 105};

Int_t  fpd_ws_map2005pp2[49]=
     {103,102,101,100,  99,  98, 88,
     107, 70, 69,105,  86,  85, 48,
      68, 67, 66, 84,  83,  82, 63,
      65, 64, 97, 81,  80,  95, 62,
      78, 77 ,76, 94,  93,  92, 61,
     106, 74, 73, 91,  90,  89, 60,
      96,111,110,109, 108,  40, 72};

// since run # 7094066 
Int_t  fpd_en_map_2006[49] =
  { 39, 38, 41, 36, 35, 34, 33, 
    7,  6,  5, 23, 22, 21, 55,
    4,  3,  2, 20, 19, 18, 54,
    1,  0, 15, 17, 16, 31, 53,
    14, 13, 12, 30, 29, 28, 52, 
    11, 10,  9, 24, 26, 25, 51,
    32, 47, 46, 45, 44, 43, 42};
Int_t fpd_es_map_2006[49]=
  {102,101,100, 99, 98, 97, 96, 
   71, 70, 69, 87, 86, 85, 48,
   68, 67, 66, 84, 83, 82, 63,
   65, 64, 79, 81, 80, 95, 61,
   78, 77 ,76, 94, 93, 92, 60,
   75, 74, 73, 91, 90, 89, 59,
   111,110,109,108,107,106,105
  };
Int_t fpd_et_map_2006[25]; // don't exist (-1)??
Int_t fpd_eb_map_2006[25]; // don't exist (-1)??

Int_t  fpd_wnsml_map_2006[36]=           // since run # 7094066 
  {  55,119,118,135,134, 53,
     117,116,115,133,132,131,
     114,113,112,130,129,128,
     127,126,125,143,142,141,
     124,123,122,140,139,138,
     52,121,120,137,136, 51
  };

Int_t fpd_wnsml_map2_2006[36]=          // since run # 7105001 
     {
       55, 23, 22, 39, 38, 53,
       21, 20, 19, 37, 36, 35,
       18, 17, 16, 34, 33, 32,
       31, 30, 29, 47, 46, 45,
       28, 27, 26, 44, 43, 42,
       52, 25, 24, 41, 40, 51
     };

Int_t fpd_wnsml_map3_2006[36]=          // since run # 7152000
     {
       55, 119,118,135,134,53,      // D00/H00/H01/I00/I01/D02
       117,116,115,133,132,131,      // H02/H03/H04/I02/I03/I04
       114,113,112,130,129,128,      // H05/H06/H07/I05/I06/I07
       127,126,125,143,142,141,      // H08/H09/H10/I08/I09/I10
       124,123,122,140,139,138,      // H11/H12/H13/I11/I12/I13
       52, 121,120,137,136,51      // D03/H14/H15/I14/I15/D01
     };

Int_t fpd_wssml_map_2006[36]=           // since run # 7094066 
     {
       63,151,150,167,166, 62,
       149,148,147,165,164,163,
       146,145,144,162,161,160,
       159,158,157,175,174,173,
       156,155,154,172,171,170,
       61,153,152,169,168, 60
     };
Int_t fpd_wssml_map2_2006[36]=          // since run # 7105001 
  {
    63, 87, 86,103,102, 62,
    85, 84, 83,101,100, 99,
    82, 81, 80, 98, 97, 96,
    95, 94, 93,111,110,109,
    92, 91, 90,108,107,106,
    61, 89, 88,105,104, 60
  };
Int_t fpd_wssml_map3_2006[36]=         // since run # 7152000
  {
    63, 151,150,167,166,62,      // D08/J00/J01/K00/K01/D09
    149,148,147,165,164,163,      // J02/J03/J04/K02/K03/K04
    146,145,144,162,161,160,      // J05/J06/J07/K05/K06/K07
    159,158,157,175,174,173,      // J08/J09/J10/K08/K09/K10
    156,155,154,172,171,170,      // J11/J12/J13/K11/K12/K13
    61, 153,152,169,168,60      // D10/J14/J15/K14/K15/D11
  };

Int_t fpd_prsen_map_2006[7]={50, 49,173,172,171,170,169}; // since run #7094066 
Int_t fpd_prses_map_2006[7]={58, 57,162,161,160,175,174}; // till run #7130004
Int_t fpd_wnlrg_map_2006[196]=          // since run # 7094066 
  {
    -1, -1,216,217,218,219,220,221,222,223,224,225, -1, -1,
    -1,248,249,250,251,252,253,254,255,256,257,258,259, -1,
    280,281,282,283,284,285,286,287,288,289,290,291,292,293,
    296,304,212,  7,  6, 23, 22, 21, 20, 39, 38,244,210,424,
    297,305,213,  5,  4, 19, 18, 17, 16, 37, 36,245,211,425,
    298,306,214,  3,  2, -1, -1, -1, -1, 35, 34,246,226,426,
    299,307,215,  1,  0, -1, -1, -1, -1, 33, 32,247,227,427,
    300,308,228, 15, 14, -1, -1, -1, -1, 47, 46,260,278,428,
    301,309,229, 13, 12, -1, -1, -1, -1, 45, 44,261,279,429,
    302,310,230, 11, 10, 31, 30, 29, 28, 43, 42,262,294,430,
    303,311,231,  9,  8, 27, 26, 25, 24, 41, 40,263,295,431,
    264,265,266,267,268,269,270,271,272,273,274,275,276,477,
    -1,232,233,234,235,236,237,238,239,240,241,242,243, -1,
    -1, -1,200,201,202,203,204,205,206,207,208,209, -1, -1
  };
Int_t fpd_wnlrg_map2_2006[196]=          // since run # 7102049 
  {
    -1, -1,216,217,218,219,220,221,222,223,224,225, -1, -1,
    -1,248,249,250,251,252,253,254,255,256,257,258,259, -1,
    280,281,282,283,284,285,286,287,288,289,290,291,292,293,
    296,304,212,472,473,489,488,490,491,504,505,244,210,424,
    297,305,213,474,475,492,493,494,495,506,507,245,211,425,
    298,306,214,476,477, -1, -1, -1, -1,508,509,246,226,426,
    299,307,215,478,479, -1, -1, -1, -1,510,511,247,227,427,
    300,308,228,480,481, -1, -1, -1, -1,512,513,260,278,428,
    301,309,229,482,483, -1, -1, -1, -1,514,515,261,279,429,
    302,310,230,484,485,496,497,498,499,516,517,262,294,430,
    303,311,231,486,487,500,501,502,503,518,519,263,295,431,
    264,265,266,267,268,269,270,271,272,273,274,275,276,277,
    -1,232,233,234,235,236,237,238,239,240,241,242,243, -1,
    -1, -1,200,201,202,203,204,205,206,207,208,209, -1, -1
  };

Int_t fpd_wnlrg_map3_2006[196]=          // since run # 7117038 
  {
    -1, -1,216,217,218,219,220,221,222,223,224,225, -1, -1,
    -1,248,249,250,251,252,253,254,255,256,257,258,259, -1,
    280,281,282,283,284,285,286,287,288,289,290,291,292,293,
    296,304,212,472,473,489,488,490,491,504,505,244,210,424,
    297,305,213,474,475,492,493,494,495,506,507,245,211,425,
    298,306,214,476,477, -1, -1, -1, -1,508,509,246,226,426,
    299,307,215,478,479, -1, -1, -1, -1,510,511,247,227,427,
    300,308,228,480,481, -1, -1, -1, -1,512,513,260,278,428,
    301,309,229,482,483, -1, -1, -1, -1,514,515,261,279,429,
    302,310,230,484,485,496,497,498,499,516,517,262,294,430,
    303,311,231,486,487,500,501,502,503,469,471,263,295,431,
    264,265,266,267,268,269,270,271,272,273,274,275,276,277,
    -1,232,233,234,235,236,237,238,239,240,241,242,243, -1,
    -1, -1,200,201,202,203,204,205,206,207,208,209, -1, -1
  };

Int_t fpd_wnlrg_map4_2006[196]=          // since run # 7126001
  {
    -1, -1,216,217,218,219,220,221,222,223,224,225, -1, -1,
    -1,248,249,250,251,252,253,254,255,256,257,258,259, -1,
    280,281,282,283,284,285,286,287,288,289,290,291,292,293,
    296,304,212,466,473,489,488,490,491,504,505,244,210,424,
    297,305,213,474,475,492,493,494,495,506,507,245,211,425,
    298,306,214,476,477, -1, -1, -1, -1,508,509,246,226,426,
    299,307,215,478,479, -1, -1, -1, -1,468,511,247,227,427,
    300,308,228,480,481, -1, -1, -1, -1,512,513,260,278,428,
    301,309,229,482,483, -1, -1, -1, -1,514,515,261,279,429,
    302,310,230,484,485,496,497,498,499,516,517,262,294,430,
    303,311,231,486,487,500,501,502,503, -1,471,263,295,431,
    264,265,266,267,268,269,270,271,272,273,274,275,276,277,
    -1,232,233,234,235,236,237,238,239,240,241,242,243, -1,
    -1, -1,200,201,202,203,204,205,206,207,208,209, -1, -1
  };

Int_t fpd_wslrg_map_2006[196]=           // since run # 7094066 
  {
    -1, -1,328,329,330,331,332,333,334,335,336,337, -1, -1,
    -1,360,361,362,363,364,365,366,367,368,369,370,371, -1,
    392,393,394,395,396,397,398,399,400,401,402,403,404,405,
    408,416,324, 71, 70, 87, 86, 85, 84,103,102,356,322,440,
    409,417,325, 69, 68, 83, 82, 81, 80,101,100,357,323,441,
    410,418,326, 67, 66, -1, -1, -1, -1, 99, 98,358,338,442,
    411,419,327, 65, 64, -1, -1, -1, -1, 97, 96,359,339,443,
    412,420,340, 79, 78, -1, -1, -1, -1,111,110,372,390,444,
    413,421,341, 77, 76, -1, -1, -1, -1,109,108,373,391,445,
    414,422,342, 75, 74, 95, 94, 93, 92,107,106,374,406,446,
    415,423,343, 73, 72, 91, 90, 89, 88,105,104,375,407,447,
    376,377,378,379,380,381,382,383,384,385,386,387,388,432,
    -1,344,345,346,347,348,349,350,351,352,353,354,355, -1,
    -1, -1,312,313,314,315,316,317,318,319,320,321, -1, -1
  };


Int_t fpd_wslrg_map2_2006[196]=          // since run # 7102049 
  {
    -1, -1,328,329,330,331,332,333,334,335,336,337, -1, -1,
    -1,360,361,362,363,364,365,366,367,368,369,370,371, -1,
    392,393,394,395,396,397,398,399,400,401,402,403,404,405,
    408,416,324,520,521,536,537,538,539,432,433,356,322,440,
    409,417,325,522,523,540,541,542,543,434,435,357,323,441,
    410,418,326,524,525, -1, -1, -1, -1,436,437,358,338,442,
    411,419,327,526,527, -1, -1, -1, -1,438,439,359,339,443,
    412,420,340,528,529, -1, -1, -1, -1,448,449,372,390,444,
    413,421,341,530,531, -1, -1, -1, -1,450,451,373,391,445,
    414,422,342,532,533,544,545,546,547,452,453,374,406,446,
    415,423,343,534,535,548,549,550,551,454,455,375,407,447,
    376,377,378,379,380,381,382,383,384,385,386,387,388, -1,
    -1,344,345,346,347,348,349,350,351,352,353,354,355, -1,
    -1, -1,312,313,314,315,316,317,318,319,320,321, -1, -1
  };


// since run # 9067058

Int_t  fpd_en_map_2008[49] =
  { 
     39, 38, 37, 36, 35, 34, 33, 
      7,  6,  5, 23, 22, 21, 55,
      4,  3,  2, 20, 19, 18, 50,
      1,  0, 15, 17, 16, 31, 53,
     14, 13, 12, 30, 29, 28, 52, 
     11, 10,  9, 27, 26, 25, 51,
     32, 47, 46, 45, 44, 43, 42
  };

Int_t fpd_es_map_2008[49]=
  {
     104,101,100, 99, 98, 97, 96, 
     71, 70, 69, 87, 86, 85, 48,
     68, 67, 66, 84, 83, 82, 62,
     65, 64, 79, 81, 80, 95, 61,
     78, 77 ,76, 94, 93, 92, 60,
     75, 74, 73, 91, 90, 89, 59,
    111,110,109,108,107,106,105
  };

//***********************************

Int_t  fpd_prswn_map2005[7]={50, 49, 70, 69, 68, 67, 66};

Int_t fpd_wt_map2005pp[25]=
     {135, 134, 133, 132, 120, 130, 129, 131, 143, 142,
     119, 118, 117, 116, 115, 114, 113, 112,
	127, 126, 125, 124, 123, 122, 121};

Int_t fpd_wb_map2005pp[25]=
     {151, 150, 149, 148, 147, 146, 145, 144,
     159, 158, 157, 156, 155, 154, 153,
	174, 166, 165, 164, 163, 162, 161, 160, 175, 168};

Int_t fpd_prswn_map2005pp[7]= {50, 49, 141, 140, 136, 138, 137};
Int_t fpd_prsws_map2005pp[7]= {56, 58, 173, 172, 171, 170, 169};

TMatrix MakeMatIII(Int_t nrows,Int_t ncols, Int_t* z)
{
  // put array z into nrows by ncols matrix

   TMatrix* zz=new TMatrix(nrows,ncols);
  for(int i=0;i<nrows;i++)
    {    
      for(int j=0;j<ncols;j++)
	{
	  (*zz)(i,j)=*z;
	  z++;
	};
    };
  return *zz;
};
TMatrix MakeMatII(Int_t n, Int_t* z)
{
  // put array z into an nxn matrix

  return MakeMatIII(n,n,z);
};

FpdMap::FpdMap(Int_t RNum)
{
  RunNum=RNum;
};
TMatrix FpdMap::GetMatrix(Int_t EW,Int_t NSTB)
{
  // EW=   1->east; 2->west
  // NSTB  1->north; 2->south; 3->top; 4->bottom; 5->Presh north; 6->Presh sout
  Int_t selectvar=(EW-1)*6+NSTB-1;
//from  GR20031112 add dau map for specific run range  
  if(RunNum>=4081126 && RunNum<=4082067)
    {
      switch (selectvar)
	{
	case ENorth: break;
	case ESouth: break;
	case ETop: break;
	case EBottom: break;
	case EpsNorth: break;
	case EpsSouth: break;
	case WNorth: break;
	case WSouth:        
	  return MakeMatII(7, fpdwsmap_1);
	  break;
	case WTop: break;
	case WBottom: break;
	case WpsNorth: break;
	case WpsSouth: 
	  return MakeMatIII(1,7,fpd_prsws_map_dau);
	}
    }
  else if ((RunNum>=6000000)&&(RunNum<=6100000))
    {
      switch (selectvar)
	{
	case ENorth:
	  return MakeMatII(7, fpdenmap_1);
	case ESouth:
	  return MakeMatII(7, fpdesmap_1);
	case ETop:
	  return MakeMatII(5, fpdetmap_1);
	case EBottom:
	  return MakeMatII(5, fpdebmap_1);
	case EpsNorth:
	  return MakeMatII(5, fpd_prsen_map);
	case EpsSouth: break;
	  return MakeMatII(5, fpd_prses_map);
	case WNorth: 
	  return MakeMatII(7, fpd_wn_map2005);
	case WSouth:        
	  break;
	case WTop:
	  return MakeMatII(5, fpd_wt_map2005);
	case WBottom:
	  break;
	case WpsNorth:
	  return MakeMatIII(1,7,fpd_prswn_map2005);
	case WpsSouth: 
	  break;
	};
    }
  else if((RunNum>6100000)&(RunNum<6115046))
    {
      switch (selectvar)
	{
	case ENorth:
	  return MakeMatII(7, fpdenmap_1);
	case ESouth:
	  return MakeMatII(7, fpdesmap_1);
	case ETop:
	  return MakeMatII(5, fpdetmap_1);
	case EBottom:
	  return MakeMatII(5, fpdebmap_1);
	case EpsNorth:
	  return MakeMatII(5, fpd_prsen_map);
	case EpsSouth: break;
	  return MakeMatII(5, fpd_prses_map);
	case WNorth: 
	  return MakeMatII(7, fpd_wn_map2005pp);
	case WSouth:        
	  return MakeMatII(7, fpd_ws_map2005pp);
	case WTop:
	  return MakeMatII(5, fpd_wt_map2005pp);
	case WBottom:
	  return MakeMatII(5, fpd_wb_map2005pp);
	case WpsNorth:
	  return MakeMatIII(1,7,fpd_prswn_map2005pp);
	case WpsSouth: 
	  return MakeMatIII(1,7,fpd_prsws_map2005pp);
	};
    }
  else if((RunNum>=6115046)&&(RunNum<7094066))
    {
      switch (selectvar)
	{
	case ENorth:
	  return MakeMatII(7, fpdenmap_1);
	case ESouth:
	  return MakeMatII(7, fpdesmap_1);
	case ETop:
	  return MakeMatII(5, fpdetmap_1);
	case EBottom:
	  return MakeMatII(5, fpdebmap_1);
	case EpsNorth:
	  return MakeMatII(5, fpd_prsen_map);
	case EpsSouth: break;
	  return MakeMatII(5, fpd_prses_map);
	case WNorth: 
	  return MakeMatII(7, fpd_wn_map2005pp);
	case WSouth:        
	  return MakeMatII(7, fpd_ws_map2005pp2);
	case WTop:
	  return MakeMatII(5, fpd_wt_map2005pp);
	case WBottom:
	  return MakeMatII(5, fpd_wb_map2005pp);
	case WpsNorth:
	  return MakeMatIII(1,7,fpd_prswn_map2005pp);
	case WpsSouth: 
	  return MakeMatIII(1,7,fpd_prsws_map2005pp);
	};      
    } 
  else if((RunNum>=7094066)&&(RunNum< 7102049))
    {
      switch (selectvar)
	{
	case ENorth:
	  return MakeMatII(7, fpd_en_map_2006);
	case ESouth:
	  return MakeMatII(7, fpd_es_map_2006);
	case ETop:
	  return TMatrix(5,5)=-1;
	case EBottom:
	  return TMatrix(5,5)=-1;;
	case EpsNorth:
	  return MakeMatII(5, fpd_prsen_map_2006);
	case EpsSouth: break;
	  return MakeMatII(5, fpd_prses_map_2006);
	case WNorth: 
	  return MakeMatII(14, fpd_wnlrg_map_2006);
	case WSouth:        
	  return MakeMatII(14, fpd_wslrg_map_2006);
	case WTop:
	  return MakeMatII(6, fpd_wnsml_map_2006);
	case WBottom:
	  return MakeMatII(6, fpd_wssml_map_2006);
	case WpsNorth:
	  return TMatrix(1,7)=-1;
	case WpsSouth: 
	  return TMatrix(1,7)=-1;
	};
    }
  else if((RunNum>=7102049)&&(RunNum< 7105001))
    {
      switch (selectvar)
	{
	case ENorth:
	  return MakeMatII(7, fpd_en_map_2006);
	case ESouth:
	  return MakeMatII(7, fpd_es_map_2006);
	case ETop:
	  return TMatrix(5,5)=-1;
	case EBottom:
	  return TMatrix(5,5)=-1;;
	case EpsNorth:
	  return MakeMatII(5, fpd_prsen_map_2006);
	case EpsSouth: break;
	  return MakeMatII(5, fpd_prses_map_2006);
	case WNorth: 
	  return MakeMatII(14, fpd_wnlrg_map2_2006);
	case WSouth:        
	  return MakeMatII(14, fpd_wslrg_map2_2006);
	case WTop:
	  return MakeMatII(6, fpd_wnsml_map_2006);
	case WBottom:
	  return MakeMatII(6, fpd_wssml_map_2006);
	case WpsNorth:
	  return TMatrix(1,7)=-1;
	case WpsSouth: 
	  return TMatrix(1,7)=-1;
	};
    }
  else if((RunNum>=7105001)&&(RunNum< 7117038))
    {
      switch (selectvar)
	{
	case ENorth:
	  return MakeMatII(7, fpd_en_map_2006);
	case ESouth:
	  return MakeMatII(7, fpd_es_map_2006);
	case ETop:
	  return TMatrix(5,5)=-1;
	case EBottom:
	  return TMatrix(5,5)=-1;;
	case EpsNorth:
	  return MakeMatII(5, fpd_prsen_map_2006);
	case EpsSouth: break;
	  return MakeMatII(5, fpd_prses_map_2006);
	case WNorth: 
	  return MakeMatII(14, fpd_wnlrg_map2_2006);
	case WSouth:        
	  return MakeMatII(14, fpd_wslrg_map2_2006);
	case WTop:
	  return MakeMatII(6, fpd_wnsml_map2_2006);
	case WBottom:
	  return MakeMatII(6, fpd_wssml_map2_2006);
	case WpsNorth:
	  return TMatrix(1,7)=-1;
	case WpsSouth: 
	  return TMatrix(1,7)=-1;
	};      
    } 
  else if((RunNum>=7117038)&&(RunNum< 7126001))
    {
      switch (selectvar)
	{
	case ENorth:
	  return MakeMatII(7, fpd_en_map_2006);
	case ESouth:
	  return MakeMatII(7, fpd_es_map_2006);
	case ETop:
	  return TMatrix(5,5)=-1;
	case EBottom:
	  return TMatrix(5,5)=-1;;
	case EpsNorth:
	  return MakeMatII(5, fpd_prsen_map_2006);
	case EpsSouth: break;
	  return MakeMatII(5, fpd_prses_map_2006);
	case WNorth: 
	  return MakeMatII(14, fpd_wnlrg_map3_2006);
	case WSouth:        
	  return MakeMatII(14, fpd_wslrg_map2_2006);
	case WTop:
	  return MakeMatII(6, fpd_wnsml_map2_2006);
	case WBottom:
	  return MakeMatII(6, fpd_wssml_map2_2006);
	case WpsNorth:
	  return TMatrix(1,7)=-1;
	case WpsSouth: 
	  return TMatrix(1,7)=-1;
	};      
    } 
  else if((RunNum>=7126001)&&(RunNum< 7130000))
    {
      switch (selectvar)
	{
	case ENorth:
	  return MakeMatII(7, fpd_en_map_2006);
	case ESouth:
	  return MakeMatII(7, fpd_es_map_2006);
	case ETop:
	  return TMatrix(5,5)=-1;
	case EBottom:
	  return TMatrix(5,5)=-1;;
	case EpsNorth:
	  return MakeMatII(5, fpd_prsen_map_2006);
	case EpsSouth: break;
	  return MakeMatII(5, fpd_prses_map_2006);
	case WNorth: 
	  return MakeMatII(14, fpd_wnlrg_map4_2006);
	case WSouth:        
	  return MakeMatII(14, fpd_wslrg_map2_2006);
	case WTop:
	  return MakeMatII(6, fpd_wnsml_map2_2006);
	case WBottom:
	  return MakeMatII(6, fpd_wssml_map2_2006);
	case WpsNorth:
	  return TMatrix(1,7)=-1;
	case WpsSouth: 
	  return TMatrix(1,7)=-1;
	};      
    }  
  else if((RunNum>=7130000)&&(RunNum< 7152000))
    {
      switch (selectvar)
	{
	case ENorth:
	  return MakeMatII(7, fpd_en_map_2006);
	case ESouth:
	  return MakeMatII(7, fpd_es_map_2006);
	case ETop:
	  return TMatrix(5,5)=-1;
	case EBottom:
	  return TMatrix(5,5)=-1;
	case EpsNorth:
	  return MakeMatII(5, fpd_prsen_map_2006);
	case EpsSouth: break;
	  return MakeMatII(5, fpd_prses_map_2006);
	case WNorth: 
	  return MakeMatII(14, fpd_wnlrg_map4_2006);
	case WSouth:        
	  return MakeMatII(14, fpd_wslrg_map2_2006);
	case WTop:
	  return MakeMatII(6, fpd_wnsml_map2_2006);
	case WBottom:
	  return MakeMatII(6, fpd_wssml_map2_2006);
	case WpsNorth:
	  return TMatrix(1,7)=-1;
	case WpsSouth: 
	  return TMatrix(1,7)=-1;
	};      
    }
  else if((RunNum>=7152000)&&(RunNum< 7158000))
    {
      switch (selectvar)
	{
	case ENorth:
	  return MakeMatII(7, fpd_en_map_2006);
	case ESouth:
	  return MakeMatII(7, fpd_es_map_2006);
	case ETop:
	  return TMatrix(5,5)=-1;
	case EBottom:
	  return TMatrix(5,5)=-1;
	case EpsNorth:
	  return MakeMatII(5, fpd_prsen_map_2006);
	case EpsSouth: break;
	  return MakeMatII(5, fpd_prses_map_2006);
	case WNorth: 
	  return MakeMatII(14, fpd_wnlrg_map4_2006);
	case WSouth:        
	  return MakeMatII(14, fpd_wslrg_map2_2006);
	case WTop:
	  return MakeMatII(6, fpd_wnsml_map3_2006);
	case WBottom:
	  return MakeMatII(6, fpd_wssml_map3_2006);
	case WpsNorth:
	  return TMatrix(1,7)=-1;
	case WpsSouth: 
	  return TMatrix(1,7)=-1;
	};      
    }
  else if((RunNum>=9067058)&&(RunNum< 9070008 ))
    {
      switch (selectvar)
	{
	case ENorth:
	  return MakeMatII(7, fpd_en_map_2008);
	case ESouth:
	  return MakeMatII(7, fpd_es_map_2008);
	case ETop:
	  return TMatrix(5,5)=-1;
	case EBottom:
	  return TMatrix(5,5)=-1;
	case EpsNorth:
	  return TMatrix(1,7)=-1;
	  return TMatrix(5,5)=-1;
	case EpsSouth: break;
	  return TMatrix(1,7)=-1;
	case WNorth: 
	  return TMatrix(5,5)=-1;
	case WSouth:        
	  return TMatrix(5,5)=-1;
	case WTop:
	  return TMatrix(5,5)=-1;
	case WBottom:
	  return TMatrix(5,5)=-1;
	case WpsNorth:
	  return TMatrix(1,7)=-1;
	case WpsSouth: 
	  return TMatrix(1,7)=-1;
	};      
    };


  return TMatrix();
};
